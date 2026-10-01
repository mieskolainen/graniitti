# Deep MLP neural network layers
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


import torch
import torch.nn as nn


class Multiply(nn.Module):
    """
    Multiplication with a non-learnable constant alpha
    """

    # Initialize one fixed output multiplier
    def __init__(self, alpha):
        super().__init__()
        self.alpha = alpha

    # Scale one tensor by the configured constant
    def forward(self, x):
        return torch.mul(x, self.alpha)


# Compute one configured Torch activation module
def get_act(act: str = "relu"):
    activations = {
        "elu": nn.ELU,
        "gelu": nn.GELU,
        "relu": nn.ReLU,
        "silu": nn.SiLU,
        "softplus": nn.Softplus,
        "tanh": nn.Tanh,
    }
    if act not in activations:
        raise Exception(f'Unknown act "{act}" chosen')
    return activations[act]()


# Append activation, normalization and dropout in the requested order
def append_post_linear_layers(
    layers: list,
    *,
    width: int,
    activation: str,
    layer_norm: bool,
    batch_norm: bool,
    dropout: float,
    act_after_norm: bool,
) -> None:
    if not act_after_norm:
        layers.append(get_act(activation))
    if layer_norm:
        layers.append(nn.LayerNorm(width))
    if batch_norm:
        layers.append(nn.BatchNorm1d(width))
    if dropout > 0:
        layers.append(nn.Dropout(dropout, inplace=False))
    if act_after_norm:
        layers.append(get_act(activation))


class LinearLayer(torch.nn.Module):
    # Initialize one dense block while retaining checkpoint field names
    def __init__(
        self,
        dim_in,
        dim_out,
        skip_connections=False,
        activation: str = "silu",
        layer_norm: bool = False,
        batch_norm: bool = False,
        dropout: float = 0.0,
        act_after_norm=True,
    ):
        super().__init__()

        self.layer = torch.nn.Linear(dim_in, dim_out)
        self.skip_connections = skip_connections
        self.act_after_norm = act_after_norm
        self.act = get_act(activation)
        operations = (
            (not act_after_norm, "act", self.act),
            (layer_norm, "ln", nn.LayerNorm(dim_out)),
            (batch_norm, "bn", nn.BatchNorm1d(dim_out)),
            (dropout > 0, "do", nn.Dropout(dropout, inplace=False)),
            (act_after_norm, "act", self.act),
        )
        self.operation_names = []
        for enabled, name, operation in operations:
            if enabled:
                setattr(self, name, operation)
                self.operation_names.append(name)

    # Apply the dense block and its optional residual connection
    def forward(self, x):
        y = self.layer(x)
        for name in self.operation_names:
            y = getattr(self, name)(y)
        if self.skip_connections and (y.shape[-1] == x.shape[-1]):
            return y + x
        return y


# Construct a dense multilayer perceptron
def MLP(
    layers: list[int],
    activation: str = "relu",
    layer_norm: bool = False,
    batch_norm: bool = False,
    dropout: float = 0.0,
    last_act: bool = False,
    skip_connections=False,
    act_after_norm=True,
):
    print(
        __name__
        + f".MLP: {layers} | activation {activation} | layer_norm {layer_norm} | batch_norm {batch_norm} | dropout {dropout} | skip_connections = {skip_connections} | act_after_norm {act_after_norm} | last_act {last_act}"
    )

    block_stop = len(layers) if last_act else len(layers) - 1
    blocks = [
        LinearLayer(
            dim_in=layers[index - 1],
            dim_out=layers[index],
            activation=activation,
            layer_norm=layer_norm,
            batch_norm=batch_norm,
            dropout=dropout,
            skip_connections=skip_connections,
            act_after_norm=act_after_norm,
        )
        for index in range(1, block_stop)
    ]
    if not last_act:
        blocks.append(nn.Linear(layers[-2], layers[-1]))
    return nn.Sequential(*blocks)


class DMLP(nn.Module):
    # Initialize the configurable dense MLP wrapper
    def __init__(
        self,
        in_dim,
        out_dim,
        mlp_dim=(128, 64),
        activation="relu",
        layer_norm=False,
        batch_norm=False,
        dropout=0.0,
        skip_connections=True,
        act_after_norm=True,
        last_tanh=False,
        last_tanh_scale=10.0,
        **kwargs,
    ):
        super().__init__()

        self.mlp = MLP(
            [in_dim, *mlp_dim, out_dim],
            activation=activation,
            skip_connections=skip_connections,
            layer_norm=layer_norm,
            batch_norm=batch_norm,
            dropout=dropout,
            act_after_norm=act_after_norm,
        )

        # Add extra final squeezing activation and post-scale aka "soft clipping"
        if last_tanh:
            self.mlp.add_module("tanh", nn.Tanh())
            self.mlp.add_module("scale", Multiply(last_tanh_scale))

    # Evaluate the wrapped dense network
    def forward(self, x):
        return self.mlp(x)
