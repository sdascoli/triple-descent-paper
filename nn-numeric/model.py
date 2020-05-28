import torch
import torch.nn.functional as F
import torch.nn as nn

class MulConstant(nn.Module):
    def __init__(self, constant=10):
        super(MulConstant, self).__init__()
        self.constant = constant

    def forward(self, x):
        return x * self.constant

    def backward(self, g): # is this necessary?
        return g * self.constant, None

class FullyConnected(nn.Module):

    def __init__(self, sigma_0=1.4142135381698608, in_dim=28*28, width=512, n_layers=1, out_dim=10, bias=True, activation='relu', ntk_scaling=False):
        super(FullyConnected, self).__init__()
        self.in_dim = in_dim
        self.width = width
        self.n_layers = n_layers
        self.bias = bias
        self.out_dim = out_dim
        
        if activation=='linear':
            self.activation = nn.Identity()
        if activation=='abs':
            self.activation = nn.PReLU(init=-1)
            self.activation.weight.requires_grad=False
        if activation=='relu':
            self.activation = nn.ReLU()
        elif activation=='tanh':
            self.activation = nn.Tanh()
        
        self.layers = []
        self.layers.append(nn.Linear(self.in_dim, self.width, bias=self.bias))
        if ntk_scaling:
            self.layers.append(MulConstant( 1 / (self.in_dim ** 0.5)))
        self.layers.append(self.activation)
        for i in range(self.n_layers-1):
            self.layers.append(nn.Linear(self.width, self.width, bias=self.bias))
            if ntk_scaling:
                self.layers.append(MulConstant( 1 / (self.width ** 0.5)))
            self.layers.append(self.activation)
        self.layers.append(nn.Linear(self.width, self.out_dim, bias=self.bias),)
        if ntk_scaling:
            self.layers.append(MulConstant( 1 / (self.width ** 0.5)))
        self.net = nn.Sequential(*self.layers)

        # NTK initialization
        if ntk_scaling:
            with torch.no_grad():
                for m in self.net.modules():
                    if isinstance(m, nn.Linear):
                        nn.init.normal_(m.weight, mean=0, std=sigma_0)
                        if m.bias is not None:
                            nn.init.constant_(m.bias, 0)
                        
    def forward(self, x):
        
        x = self.net(x)
        
        return x
