# Torch model training functions
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy

import numpy as np
import torch
from tqdm import tqdm

from core.inference.lzmlp import LZMLP

# Share one local ISO timestamp across surrogate trainers
from core.tune.runtime.process import training_timestamp


# Apply the requested batch reduction to one per-row loss
def _reduce_batch_loss(loss: torch.Tensor, reduction: str) -> torch.Tensor:
    if reduction == "mean":
        return torch.mean(loss)
    if reduction == "sum":
        return torch.sum(loss)
    if reduction == "none":
        return loss
    raise ValueError(f"Unknown loss reduction: {reduction}")


# Align scalar targets by event before applying componentwise regression losses
def _regression_values(predictions, targets, errs):
    predictions = predictions[:, None] if predictions.ndim == 1 else predictions
    targets = targets[:, None] if targets.ndim == 1 else targets
    if predictions.shape != targets.shape:
        raise ValueError("Neural predictions and targets must have matching event and component dimensions")
    if errs is not None and errs.ndim == 1:
        errs = errs[:, None]
    return predictions, targets, errs


# Compute the uncertainty-weighted mean-square error
def mse_error_loss(
    predictions: torch.Tensor, targets: torch.Tensor, errs: torch.Tensor = None, reduction="sum", EPS: float = 1e-12
):
    predictions, targets, errs = _regression_values(predictions, targets, errs)
    # Take mean over component dimensions
    loss = (predictions - targets) ** 2
    if errs is not None:
        loss = loss / (errs**2 + EPS)
    if loss.ndim > 1:
        loss = torch.mean(loss, dim=-1)

    return _reduce_batch_loss(loss, reduction)


# Compute the uncertainty-weighted mean-absolute error
def mae_error_loss(
    predictions: torch.Tensor, targets: torch.Tensor, errs: torch.Tensor = None, reduction="sum", EPS: float = 1e-12
):
    predictions, targets, errs = _regression_values(predictions, targets, errs)
    # Take mean over component dimensions
    loss = torch.abs(predictions - targets)
    if errs is not None:
        loss = loss / (errs + EPS)
    if loss.ndim > 1:
        loss = torch.mean(loss, dim=-1)

    return _reduce_batch_loss(loss, reduction)


# Compute the requested neural surrogate loss function
def neural_loss(name: str):
    losses = {"MAE": mae_error_loss, "MSE": mse_error_loss}
    try:
        return losses[name]
    except KeyError as exc:
        raise ValueError(f"Unknown loss function: {name}") from exc


# Build the common LZMLP architecture parameters
def lzmlp_parameters(in_dim: int, out_dim: int, hidden_dim: int, hidden_layers: int) -> dict:
    return {
        "in_dim": int(in_dim),
        "out_dim": int(out_dim),
        "mlp_dim": [int(hidden_dim)] * int(hidden_layers),
        "skip_connections": True,
        "activation": "silu",
        "layer_norm": True,
        "batch_norm": False,
    }


# Build one tensor dataset with optional target uncertainties
def _tensor_dataset(X: np.ndarray, y: np.ndarray, errs: np.ndarray | None, device: str, dtype: torch.dtype):
    tensors = [torch.tensor(values, dtype=dtype, device=device) for values in (X, y)]
    if errs is not None:
        tensors.append(torch.tensor(errs, dtype=dtype, device=device))
    return torch.utils.data.TensorDataset(*tensors)


# Unpack one supervised tensor batch with optional uncertainties
def _batch_values(batch) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor | None]:
    return batch[0], batch[1], batch[2] if len(batch) == 3 else None


# Evaluate a model on one DataLoader with the shared validation loss convention
def _evaluate_loader(model: torch.nn.Module, val_loader, lossfunc: callable) -> float:
    model.eval()
    eval_loss = 0.0
    num_samples = 0

    with torch.no_grad():
        for batch in val_loader:
            X_batch, y_batch, errs_batch = _batch_values(batch)

            pred = model(X_batch)
            loss = lossfunc(pred, y_batch, errs_batch, "sum")
            if not bool(torch.isfinite(loss)):
                raise FloatingPointError("Neural validation loss is non-finite")

            eval_loss += loss.item()
            num_samples += X_batch.shape[0]

    if num_samples == 0:
        raise Exception("No validation batches processed (batch_size > num of events ?)")

    return eval_loss / num_samples


# Train the Torch surrogate model with validation-based best snapshot selection
def train_torch(
    model: torch.nn.Module,
    optimizer,
    scheduler,
    lossfunc: callable,
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_val: np.ndarray,
    y_val: np.ndarray,
    errs_train: np.ndarray = None,
    errs_val: np.ndarray = None,
    num_epochs: int = 100,
    batch_size: int = 256,
    lipschitz: float = None,
    grad_clip: float = 10.0,
    device: str = "cpu",
    dtype: torch.dtype = torch.float32,
    patience: int = None,
    min_delta: float = 0.0,
):
    print(
        f"{training_timestamp()} {__name__}.train_torch: N_train = {X_train.shape[0]}, N_validation = {X_val.shape[0]}"
    )

    train_dataset = _tensor_dataset(X_train, y_train, errs_train, device, dtype)
    val_dataset = _tensor_dataset(X_val, y_val, errs_val, device, dtype)

    # Train
    train_loader = torch.utils.data.DataLoader(train_dataset, batch_size=batch_size, shuffle=True, drop_last=False)

    # Validation
    val_loader = torch.utils.data.DataLoader(val_dataset, batch_size=batch_size, shuffle=False, drop_last=False)

    best_model = None
    best_val_loss = float("inf")
    best_epoch = None
    stopped_epoch = None
    stale_epochs = 0

    stats = {
        "train_loss": [],
        "eval_loss": [],
        "initial_eval_loss": None,
        "best_val_loss": None,
        "best_epoch": None,
        "stopped_epoch": None,
    }

    initial_eval_loss = _evaluate_loader(model, val_loader, lossfunc)
    best_model = copy.deepcopy(model)
    best_val_loss = initial_eval_loss
    best_epoch = 0
    stats["initial_eval_loss"] = initial_eval_loss
    print(f"{training_timestamp()} Initial Val Loss: {initial_eval_loss:.4f}")

    for epoch in tqdm(range(num_epochs)):
        model.train()  #!
        train_loss = 0.0
        num_samples = 0

        ## --------------------------------------------------------
        ## Training
        ## --------------------------------------------------------

        for batch in train_loader:
            X_batch, y_batch, errs_batch = _batch_values(batch)

            optimizer.zero_grad()
            pred = model(X_batch)
            loss = lossfunc(pred, y_batch, errs_batch, "mean")

            # Lipschitz smoothness
            if lipschitz is not None and lipschitz > 0.0 and hasattr(model, "get_lipschitz_loss"):
                loss = loss + lipschitz * model.get_lipschitz_loss()
            if not bool(torch.isfinite(loss)):
                raise FloatingPointError("Neural training loss is non-finite")
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), grad_clip, error_if_nonfinite=True)
            optimizer.step()

            train_loss += loss.item() * X_batch.shape[0]
            num_samples += X_batch.shape[0]

        if num_samples == 0:
            raise Exception("No training batches processed (batch_size > num of events ?)")

        train_loss /= num_samples
        scheduler.step()

        stats["train_loss"].append(train_loss)

        ## --------------------------------------------------------
        ## Validation
        ## --------------------------------------------------------

        eval_loss = _evaluate_loader(model, val_loader, lossfunc)

        if eval_loss < best_val_loss - min_delta:
            best_model = copy.deepcopy(model)
            best_val_loss = eval_loss
            best_epoch = epoch + 1
            stale_epochs = 0

            print(
                f"{training_timestamp()} Epoch {epoch + 1}/{num_epochs} | "
                f"Train Loss: {train_loss:.4f}, Val Loss: {eval_loss:.4f} | "
                f"lr: {scheduler.get_last_lr()[0]:.2e}"
            )
        else:
            stale_epochs += 1

        stats["eval_loss"].append(eval_loss)

        if patience is not None and patience > 0 and stale_epochs >= patience:
            stopped_epoch = epoch + 1
            best_epoch_label = "initial" if best_epoch == 0 else str(best_epoch)
            print(
                f"{training_timestamp()} Early stopping neural training at "
                f"epoch {stopped_epoch}; best validation loss was "
                f"{best_val_loss:.4f} at epoch {best_epoch_label}"
            )
            break

    if best_model is None:
        best_model = copy.deepcopy(model)
        best_val_loss = stats["eval_loss"][-1] if stats["eval_loss"] else float("inf")
        best_epoch = len(stats["eval_loss"])

    stats["best_val_loss"] = best_val_loss
    stats["best_epoch"] = best_epoch
    stats["stopped_epoch"] = stopped_epoch
    return best_model, stats


# Refit a selected neural architecture for a fixed number of epochs
def train_torch_fixed_epochs(
    model: torch.nn.Module,
    lossfunc: callable,
    X: np.ndarray,
    y: np.ndarray,
    errs: np.ndarray = None,
    num_epochs: int = 1,
    batch_size: int = 256,
    lr: float = 1.0e-2,
    weight_decay: float = 1.0e-3,
    gamma: float = 1.0e-4,
    lipschitz: float = None,
    grad_clip: float = 10.0,
    device: str = "cpu",
    dtype: torch.dtype = torch.float32,
) -> torch.nn.Module:
    if num_epochs < 1:
        return model.eval().requires_grad_(False)
    dataset = _tensor_dataset(X, y, errs, device, dtype)
    loader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=True, drop_last=False)
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
    scheduler = torch.optim.lr_scheduler.ExponentialLR(optimizer, gamma=1.0 - gamma)
    log_interval = max(1, int(num_epochs) // 10)
    for epoch in range(int(num_epochs)):
        model.train()
        total_loss = 0.0
        total_rows = 0
        for batch in loader:
            X_batch, y_batch, errs_batch = _batch_values(batch)
            optimizer.zero_grad(set_to_none=True)
            loss = lossfunc(model(X_batch), y_batch, errs_batch, "mean")
            if lipschitz is not None and lipschitz > 0.0 and hasattr(model, "get_lipschitz_loss"):
                loss = loss + lipschitz * model.get_lipschitz_loss()
            if not bool(torch.isfinite(loss)):
                raise FloatingPointError("Neural training loss is non-finite")
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), grad_clip, error_if_nonfinite=True)
            optimizer.step()
            total_loss += float(loss.detach()) * len(X_batch)
            total_rows += len(X_batch)
        scheduler.step()
        if epoch % log_interval == 0 or epoch + 1 == int(num_epochs):
            mean_loss = total_loss / max(1, total_rows)
            print(f"{training_timestamp()} Neural full refit epoch {epoch + 1}/{num_epochs}: loss = {mean_loss:.6f}")
    model.eval()
    for parameter in model.parameters():
        parameter.requires_grad_(False)
    return model


# Select a neural model on validation data and refit it on all retained rows
def train_validated_lzmlp(
    *,
    model_param: dict,
    training: tuple[np.ndarray, np.ndarray, np.ndarray | None],
    validation: tuple[np.ndarray, np.ndarray, np.ndarray | None],
    full: tuple[np.ndarray, np.ndarray, np.ndarray | None],
    options: dict,
    use_errors: bool,
    print_model: bool = False,
) -> tuple[torch.nn.Module, torch.nn.Module, dict]:
    X_train, y_train, errs_train = training
    X_val, y_val, errs_val = validation
    X_full, y_full, errs_full = full
    if not use_errors:
        errs_train = errs_val = errs_full = None
    dtype = options.get("dtype", torch.float64)
    model = LZMLP(**model_param).to(device=options["device"], dtype=dtype)
    if print_model:
        print(model)
    initial_state = copy.deepcopy(model.state_dict())
    optimizer = torch.optim.AdamW(model.parameters(), lr=options["lr"], weight_decay=options["weight_decay"])
    scheduler = torch.optim.lr_scheduler.ExponentialLR(optimizer, gamma=1.0 - options["gamma"])
    lossfunc = neural_loss(options["loss_name"])
    selected_model, stats = train_torch(
        model=model,
        optimizer=optimizer,
        scheduler=scheduler,
        lossfunc=lossfunc,
        X_train=X_train,
        y_train=y_train,
        X_val=X_val,
        y_val=y_val,
        errs_train=errs_train,
        errs_val=errs_val,
        num_epochs=options["num_epochs"],
        batch_size=options["batch_size"],
        lipschitz=options["lipschitz"],
        device=options["device"],
        dtype=dtype,
        patience=options["patience"],
        min_delta=options["min_delta"],
    )
    full_model = LZMLP(**model_param).to(device=options["device"], dtype=dtype)
    full_model.load_state_dict(initial_state)
    full_model = train_torch_fixed_epochs(
        model=full_model,
        lossfunc=lossfunc,
        X=X_full,
        y=y_full,
        errs=errs_full,
        num_epochs=int(stats.get("best_epoch") or 0),
        batch_size=options["batch_size"],
        lr=options["lr"],
        weight_decay=options["weight_decay"],
        gamma=options["gamma"],
        lipschitz=options["lipschitz"],
        device=options["device"],
        dtype=dtype,
    )
    return selected_model, full_model, stats
