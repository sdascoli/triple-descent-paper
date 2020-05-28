import numpy as np
import torch
import torch.nn as nn
import matplotlib.pyplot as plt
import scipy
from collections import defaultdict
from utils import rot
from utils import get_data, hinge_regression, hinge_classification
from model import FullyConnected
import argparse

def train_and_test(model, tr_data, te_data, crit, task, opt, epochs, checkpoints):
    
    tr_losses = []
    te_losses = []
    te_accs = []
    for epoch in range(epochs):
        epoch_loss = 0
        for (x,y) in tr_data:
            x, y = x.to(device), y.to(device)
            opt.zero_grad()
            out = model(x)
            loss = crit(out, y)
            loss.backward()
            epoch_loss += loss.item()/len(tr_data)
            opt.step()
        if epoch in checkpoints:
            tr_losses.append(epoch_loss)
            te_loss, te_acc = test(model, te_data, crit, task)
            te_losses.append(te_loss)
            te_accs.append(te_acc)
    return tr_losses, te_losses, te_accs

def test(model, te_data, crit, task):
    with torch.no_grad():
        for (x,y) in te_data:
            x, y = x.to(device), y.to(device)
            out = model(x)
            test_loss = crit(out, y).item()
            if task=='classification':
                preds = out.max(1)[1]
                test_acc = preds.eq(y).sum().float()/len(y)
                test_acc = test_acc.item()
            else:
                test_acc = 0
            break
    return test_loss, test_acc

def test_ensemble(models, te_data, crit, task):
    with torch.no_grad():
        for (x,y) in te_data:
            x, y = x.to(device), y.to(device)
            outs = torch.stack([model(x) for model in models])
            out = outs.mean(dim=0)
            test_loss = crit(out, y).item()
            if task=='classification':
                preds = out.max(1)[1]
                test_acc = preds.eq(y).sum().float()/len(y)
                test_acc = test_acc.item()
            else:
                test_acc = 0
            break
    return test_loss, test_acc

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--name', default=None, type=str)
    parser.add_argument('--num_seeds', default=1, type=int)
    parser.add_argument('--no_cuda', default=False, type=bool)

    parser.add_argument('--task', default='classification', type=str)
    parser.add_argument('--dataset', default='random', type=str)
    parser.add_argument('--loss_type', default='default', type=str)

    parser.add_argument('--depth', default=1, type=int)
    parser.add_argument('--teacher_width', default=100, type=int)
    parser.add_argument('--teacher_depth', default=2, type=int)
    parser.add_argument('--width', default=20, type=int)
    parser.add_argument('--activation', default='relu', type=str)
    
    parser.add_argument('--epochs', default=1000, type=int)
    parser.add_argument('--d', default=2, type=int)
    parser.add_argument('--n', default=100, type=int)
    parser.add_argument('--n_test', default=1000, type=int)
    parser.add_argument('--noise', default=0.1, type=float)
    parser.add_argument('--test_noise', default=True, type=bool)
    parser.add_argument('--n_classes', default=None, type=int)
    parser.add_argument('--lr', default=0.01, type=float)
    parser.add_argument('--mom', default=0.9, type=float)
    parser.add_argument('--wd', default=0., type=float)
    parser.add_argument('--bs', default=1000000, type=int)
    
    args = parser.parse_args()

    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    if args.no_cuda: device='cpu'

    if not args.n_classes: args.n_classes = 1 if args.task=='regression' else 2
    if not args.loss_type: args.loss_type = 'mse' if args.task=='regression' else 'nll'
    
    if args.task=='classification':
        if args.loss_type == 'linear_hinge':
            crit = lambda x,y : hinge_classification(x,y, type='linear')
        elif args.loss_type == 'quadratic_hinge':
            crit = lambda x,y : hinge_classification(x,y, type='quadratic')
        elif args.loss_type == 'nll':
            crit = nn.CrossEntropyLoss()
        else:
            raise NotImplementedError

    elif args.task=='regression':
        if args.loss_type == 'linear_hinge':
            crit = lambda x,y : hinge_regression(x,y, epsilon=args.epsilon, type='linear')
        elif args.loss_type == 'quadratic_hinge':
            crit = lambda x,y : hinge_regression(x,y, epsilon=args.epsilon, type='quadratic')
        elif args.loss_type == 'mse':
            crit = nn.MSELoss()
        else:
            raise NotImplementedError
    else:
        raise 

    torch.manual_seed(0)
    with torch.no_grad():
        teacher = FullyConnected(width=args.teacher_width, n_layers=args.teacher_depth, in_dim=args.d, out_dim=args.n_classes, activation=args.activation).to(device)

    bs = min(args.bs, args.n)
    n_batches = int(args.n/bs)
    tr_data = get_data(args.dataset, args.task, n_batches, bs, args.d, args.noise, n_classes=args.n_classes, teacher=teacher)
    test_noise = args.noise if args.test_noise else 0
    te_data = get_data(args.dataset, args.task, 1, args.n_test, args.d, test_noise, n_classes=args.n_classes, teacher=teacher)

    tr_losses = []
    te_losses  = []
    te_accs    = []

    students = []
    checkpoints = np.unique(np.logspace(0,np.log10(args.epochs),20).astype(int))
    for seed in range(args.num_seeds):
        torch.manual_seed(seed)
        student = FullyConnected(width=args.width, n_layers=args.depth, in_dim=args.d, out_dim=args.n_classes, activation=args.activation).to(device)
        opt = torch.optim.SGD(student.parameters(), lr=args.lr, momentum=args.mom, weight_decay=args.wd)
        tr_loss_hist, te_loss_hist, te_acc_hist = train_and_test(student, tr_data, te_data, crit, args.task, opt, args.epochs, checkpoints)
        tr_losses.append(tr_loss_hist)
        te_losses.append(te_loss_hist)
        te_accs.append(te_acc_hist)
        students.append(student)

    tr_losses, te_losses, te_accs = np.array(tr_losses), np.array(te_losses), np.array(te_accs)
    tr_loss, te_loss, te_acc = np.mean(tr_losses, axis=0), np.mean(te_losses, axis=0), np.mean(te_accs, axis=0)
    te_loss_ens, te_acc_ens = test_ensemble(students, te_data, crit, args.task)   
    
    dic = {'args':args, 'checkpoints':checkpoints,
           'tr_loss':tr_loss, 'te_loss':te_loss, 'te_acc':te_acc,
           'te_loss_ens':te_loss_ens, 'te_acc_ens':te_acc_ens}
    print(dic)
    torch.save(dic, args.name+'.pyT')
