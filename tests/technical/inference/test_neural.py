# Numerical tests for neural surrogate training and Lipschitz layers
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import numpy as np
import pytest
import torch
from core.inference import neural
from core.inference.lzmlp import LZMLP, LipschitzLinear


# Check zero affine rows retain finite gradients and can learn nonzero weights
@pytest.mark.parametrize("dtype", [torch.float32, torch.float64])
def test_lipschitz_zero_rows_finite_gradients(dtype):
    layer = LipschitzLinear(3, 2).to(dtype=dtype)
    with torch.no_grad():
        layer.weight.zero_()
    inputs = torch.tensor([[1.0, 2.0, 3.0]], dtype=dtype, requires_grad=True)
    layer(inputs).sum().backward()
    for parameter in layer.parameters():
        assert torch.isfinite(parameter.grad).all()
    torch.testing.assert_close(layer.weight.grad, inputs.detach().expand(2, -1))
    assert torch.isfinite(inputs.grad).all()


# Check initialization preserves the affine map and bounds subsequent large weights
def test_lipschitz_affine_bound():
    torch.manual_seed(21)
    layer = LipschitzLinear(3, 2)
    torch.testing.assert_close(layer.get_lipschitz_constant().squeeze(), layer.weight.abs().sum(1).max())
    with torch.no_grad():
        layer.weight.mul_(100.0)
    inputs = torch.eye(3)
    effective_weight = (layer(inputs) - layer(torch.zeros_like(inputs))).T
    assert torch.all(effective_weight.abs().sum(1) <= layer.get_lipschitz_constant() * (1.0 + 1.0e-6))


# Check scalar and one-column regression losses have identical batch reductions
@pytest.mark.parametrize("loss", [neural.mse_error_loss, neural.mae_error_loss])
@pytest.mark.parametrize("reduction", ["sum", "mean", "none"])
def test_neural_scalar_targets_match_columns(loss, reduction):
    predicted = torch.tensor([0.0, 2.0, 5.0])
    target = torch.tensor([1.0, 1.0, 2.0])
    errors = torch.tensor([0.5, 1.0, 2.0])
    actual = loss(predicted, target, errors, reduction)
    expected = loss(predicted[:, None], target[:, None], errors[:, None], reduction)
    torch.testing.assert_close(actual, expected)


# Train the actual Lipschitz network with its optional penalty disabled
@pytest.mark.parametrize("penalty", [None, 0.0])
def test_lipschitz_training_without_penalty(penalty):
    torch.manual_seed(5)
    model = LZMLP(1, 1, mlp_dim=(4,))
    optimizer = torch.optim.SGD(model.parameters(), lr=0.01)
    scheduler = torch.optim.lr_scheduler.ExponentialLR(optimizer, gamma=1.0)
    x = np.linspace(-1.0, 1.0, 8, dtype=np.float32)[:, None]
    selected, stats = neural.train_torch(
        model,
        optimizer,
        scheduler,
        neural.mse_error_loss,
        x[:6],
        x[:6] ** 2,
        x[6:],
        x[6:] ** 2,
        num_epochs=2,
        batch_size=3,
        lipschitz=penalty,
    )
    assert len(stats["train_loss"]) == 2
    assert np.isfinite(stats["train_loss"]).all()
    assert torch.isfinite(selected(torch.tensor(x))).all()


# Align column predictions and vector targets before reducing regression errors
@pytest.mark.parametrize("loss", [neural.mse_error_loss, neural.mae_error_loss])
def test_column_predictions_match_vector_targets(loss):
    predictions = torch.tensor([[0.0], [2.0], [4.0]], requires_grad=True)
    targets = torch.tensor([0.0, 2.0, 4.0])
    errors = torch.tensor([0.5, 1.0, 2.0])
    value = loss(predictions, targets, errors, "mean")
    value.backward()
    torch.testing.assert_close(value, torch.zeros_like(value))
    torch.testing.assert_close(predictions.grad, torch.zeros_like(predictions))


# Preserve deterministic inference when validation selects the initial neural state
def test_neural_zero_epoch_refit_is_frozen():
    model = torch.nn.Sequential(torch.nn.Linear(1, 2), torch.nn.Dropout(0.5), torch.nn.Linear(2, 1))
    data = np.ones((4, 1), dtype=np.float32)
    result = neural.train_torch_fixed_epochs(model, neural.mse_error_loss, data, data, num_epochs=0)
    assert not result.training
    assert all(not parameter.requires_grad for parameter in result.parameters())
    torch.testing.assert_close(result(torch.tensor(data)), result(torch.tensor(data)))


# A failed full-data refit must not yield a frozen model with invalid predictions
def test_refit_nonfinite_training_loss():
    model = torch.nn.Linear(1, 1)
    x = np.ones((2, 1), dtype=np.float32)
    y = np.full_like(x, np.nan)
    with pytest.raises(FloatingPointError, match="training loss"):
        neural.train_torch_fixed_epochs(model, neural.mse_error_loss, x, y, num_epochs=1)


# Invalid validation predictions must not silently select the initial model
def test_training_nonfinite_validation():
    model = torch.nn.Linear(1, 1)
    optimizer = torch.optim.SGD(model.parameters(), lr=0.01)
    scheduler = torch.optim.lr_scheduler.ExponentialLR(optimizer, gamma=1.0)
    x = np.ones((2, 1), dtype=np.float32)
    with pytest.raises(FloatingPointError, match="validation loss"):
        neural.train_torch(
            model, optimizer, scheduler, neural.mse_error_loss, x, x, x, np.full_like(x, np.nan), num_epochs=1
        )
