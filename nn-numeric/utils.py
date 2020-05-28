# some useful functions
import os
import shutil
import math
import torch
import torchvision
from torch.utils.data import Dataset, DataLoader
import numpy as np

def hinge_regression(output, target, epsilon=.1, type='quadratic'):
    power = 1 if type=='linear' else 2
    delta = (output-target).abs()-epsilon
    loss = torch.nn.functional.relu(delta)*delta.pow(power)
    return loss.mean()

def hinge_classification(output,target,epsilon=.5, type='quadratic'):
    power = 1 if type=='linear' else 2
    output_size=output.size(1)
    if output_size==1:
        target = 2*target.double()-1
        print(target,output)
        return 0.5*(epsilon-output*target).mean()
    delta = torch.zeros(output.size(0))
    for i,(out,tar) in enumerate(zip(output,target)):
        tar = int(tar)
        delta[i] = epsilon + torch.cat((out[:tar],out[tar+1:])).max() - out[tar]
    loss = 0.5 * torch.nn.functional.relu(delta).pow(power).mean()
    return loss
    
def normalize(x):
    mean = x.mean(dim=0, keepdim=True)
    std = x.std(dim=0, keepdim=True)
    std[std==0]=1
    return (x-mean)/std

def get_data(dataset, task, n_batches, bs, d, noise, var=.5, n_classes=None, teacher=None, train=True):

    n = n_batches*bs
    if dataset=='random':
        data = torch.randn(n, d)
    else:
        assert d**0.5 == int(d**0.5)
        transforms = torchvision.transforms.Compose([
                    torchvision.transforms.Resize(int(d**0.5)), 
                    torchvision.transforms.ToTensor()]) 
        dataset = getattr(torchvision.datasets, dataset.upper())('~/data', train=train, download=True, transform=transforms)
        data = normalize(torch.cat([dataset[mu][0] for mu in range(n)]).view(n,d))
    data = data*d**0.5/data.norm(dim=-1,keepdim=True)
        
    with torch.no_grad():
        dataset = []
        # Gaussian Mixture
        # if task=='classification':
        #     for i in range(n_batches):
        #         vectors = torch.randn(n_classes, d)
        #         if n_classes==2:
        #             vectors[0] = torch.ones(d)
        #             vectors[1] = -torch.ones(d)                    
        #         labels = torch.randint(n_classes,(bs,))
        #         x = torch.ones(bs,d)
        #         y = torch.ones(bs)
        #         for j, label in enumerate(labels):
        #             x[j] = vectors[label] + var * torch.randn(d)
        #             y[j] = label if np.random.random()>noise else np.random.randint(n_classes)           
        #         dataset.append((x,y.long()))
        if task=='classification':
            for i in range(n_batches):
                x = data[i*bs:(i+1)*bs]
                y = teacher(x).max(1)[1].squeeze()
                for j in range(len(y)):
                    if np.random.random()<noise:
                        y[j]= np.random.randint(n_classes)
                dataset.append((x,y.long()))
        elif task=='regression':
            for i in range(n_batches):
                x = data[i*bs:(i+1)*bs]
                y = teacher(x)+noise*torch.randn((bs,1))
                dataset.append((x,y))
    return dataset

def rot(x, th):
    with torch.no_grad(): 
        rotation = torch.eye(len(x))
        rotation[:2,:2] = torch.Tensor([[np.cos(th),np.sin(th)],[-np.sin(th), np.cos(th)]]) 
        return rotation @ x

def who_am_i():
    import subprocess
    whoami = subprocess.run(['whoami'], stdout=subprocess.PIPE)
    whoami = whoami.stdout.decode('utf-8')
    whoami = whoami.strip('\n')
    return whoami

def copy_py(dst_folder):
    # and copy all .py's into dst_folder
    if not os.path.exists(dst_folder):
        print("Folder doesn't exist!")
        return
    for f in os.listdir():
        if f.endswith('.py'):
            shutil.copy2(f, dst_folder)

class FastMNIST(torchvision.datasets.MNIST):
    def __init__(self, *args, **kwargs): 
        super().__init__(*args, **kwargs)

        if torch.cuda.is_available():
            self.device = torch.device('cuda')
        else:
            self.device = torch.device('cpu')

        self.data = self.data.unsqueeze(1).float().div(255)
        self.data = self.data.sub_(self.data.mean()).div_(self.data.std())
        self.data, self.targets = self.data.to(self.device), self.targets.to(self.device)
    def __getitem__(self, index):
        img, target = self.data[index], self.targets[index]
        return img, target, index

class FastFashionMNIST(torchvision.datasets.FashionMNIST):
    def __init__(self, *args, **kwargs): 
        super().__init__(*args, **kwargs)

        if torch.cuda.is_available():
            self.device = torch.device('cuda')
        else:
            self.device = torch.device('cpu')

        self.data = self.data.unsqueeze(1).float().div(255)
        self.data = self.data.sub_(self.data.mean()).div_(self.data.std())
        self.data, self.targets = self.data.to(self.device), self.targets.to(self.device)
    def __getitem__(self, index):
        img, target = self.data[index], self.targets[index]
        return img, target, index
    
class FastCIFAR10(torchvision.datasets.CIFAR10):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if torch.cuda.is_available():
            self.device = torch.device('cuda')
        else:
            self.device = torch.device('cpu')
        self.data = torch.from_numpy(self.data).float().div(255)
        self.data = self.data.sub_(self.data.mean()).div_(self.data.std())
        self.data, self.targets = self.data.to(self.device), torch.LongTensor(self.targets).to(self.device)
    def __getitem__(self, index):
        img, target = self.data[index], self.targets[index]
        return img, target, index

def get_pca(tr_data, te_data, input_size, normalized = True):
    device = tr_data.device
    
    tr_data.data = tr_data.data.view(tr_data.data.size(0),-1)
    te_data.data = te_data.data.view(te_data.data.size(0),-1)
    x = tr_data.data.cpu()
    # DATA IS ALREADY NORMALIZED
    # m = x.mean(0).expand_as(x)
    # u,s,v = torch.svd(torch.t(x-m))
    u,s,v = torch.svd(torch.t(x))
    if normalized:
        tr_data.data = (tr_data.data) @ u[:, :input_size].to(device) / s[:input_size].to(device) ** 0.5
        te_data.data = (te_data.data) @ u[:, :input_size].to(device) / s[:input_size].to(device) ** 0.5
    else:
        tr_data.data = (tr_data.data) @ u[:, :input_size].to(device).to(device) ** 0.5
        te_data.data = (te_data.data) @ u[:, :input_size].to(device).to(device) ** 0.5

    return tr_data, te_data
