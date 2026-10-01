# Lipschitz MLP neural network
#
# https://arxiv.org/abs/2202.08345
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math

import torch
import torch.nn as nn

from core.inference.dmlp import Multiply, append_post_linear_layers


class LipschitzLinear(torch.nn.Module):
    """Lipschitz linear layer"""

    # Initialize one Lipschitz-constrained affine layer
    def __init__(self, in_features, out_features):
        super().__init__()

        self.in_features = in_features
        self.out_features = out_features
        self.weight = torch.nn.Parameter(torch.empty((out_features, in_features), requires_grad=True))
        self.bias = torch.nn.Parameter(torch.empty((out_features), requires_grad=True))
        self.c = torch.nn.Parameter(torch.empty((1), requires_grad=True))
        self.softplus = torch.nn.Softplus()
        self.initialize_parameters()

    # Initialize affine weights and the learnable Lipschitz bound
    def initialize_parameters(self):
        stdv = 1.0 / math.sqrt(self.weight.size(1))
        with torch.no_grad():
            self.weight.uniform_(-stdv, stdv)
            self.bias.uniform_(-stdv, stdv)
            bound = self.weight.abs().sum(1).max()
            self.c.copy_(bound + torch.log(-torch.expm1(-bound)))

    # Compute the positive learned Lipschitz constant
    def get_lipschitz_constant(self):
        return self.softplus(self.c)

    # Apply the row-normalized Lipschitz affine transformation
    def forward(self, input):
        lipc = self.softplus(self.c)
        row_norm = torch.abs(self.weight).sum(1)
        denominator = torch.maximum(row_norm, lipc).clamp_min(torch.finfo(self.weight.dtype).tiny)
        scale = lipc / denominator
        return torch.nn.functional.linear(input, self.weight * scale.unsqueeze(1), self.bias)


class LZMLP(torch.nn.Module):
    # Initialize a Lipschitz MLP with shared dense post-processing rules
    def __init__(
        self,
        in_dim,
        out_dim,
        mlp_dim=(128, 64),
        activation="relu",
        layer_norm=False,
        batch_norm=False,
        dropout=0.0,
        last_tanh=False,
        last_tanh_scale=10.0,
        act_after_norm=True,
        **kwargs,
    ):
        super().__init__()
        dimensions = [in_dim, *mlp_dim, out_dim]
        layers = []
        for index, (dim_in, dim_out) in enumerate(zip(dimensions[:-1], dimensions[1:], strict=True)):
            layers.append(LipschitzLinear(dim_in, dim_out))
            if index < len(dimensions) - 2:
                append_post_linear_layers(
                    layers,
                    width=dim_out,
                    activation=activation,
                    layer_norm=layer_norm,
                    batch_norm=batch_norm,
                    dropout=dropout,
                    act_after_norm=act_after_norm,
                )
        self.mlp = nn.Sequential(*layers)
        if last_tanh:
            self.mlp.add_module("tanh", nn.Tanh())
            self.mlp.add_module("scale", Multiply(last_tanh_scale))

    # Compute the product Lipschitz regularizer across affine layers
    def get_lipschitz_loss(self):
        return math.prod(layer.get_lipschitz_constant() for layer in self.mlp if type(layer) is LipschitzLinear)

    # Evaluate the Lipschitz network
    def forward(self, x):
        return self.mlp(x)
