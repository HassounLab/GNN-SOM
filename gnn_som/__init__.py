import torch.nn as nn
import torch_geometric.nn as gnn
from torch_geometric import __version__ as pygVersion

def createGnnSom(convName, width, depth, featureCount):
    if convName == 'cheb':
        conv = lambda a, b: gnn.ChebConv(a, b, K=5)
    elif convName == 'cheb10k':
        conv = lambda a, b: gnn.ChebConv(a, b, K=10)
    elif convName == 'cheb15k':
        conv = lambda a, b: gnn.ChebConv(a, b, K=15)
    else:
        raise Exception('Unknown convolution operator ' + convName)

    layerSizes = [featureCount] + [width] * depth + [1]

    modules = []
    for i in range(len(layerSizes) - 1):
        modules.append((conv(layerSizes[i], layerSizes[i + 1]), 'x, edge_index -> x'))
        if i != len(layerSizes) - 2:
            modules.append(nn.ReLU())
            modules.append(nn.Dropout(0.5))
    return gnn.Sequential('x, edge_index', modules)

def loadGnnSomState(model, state):
    if pygVersion.startswith('2.'): # convert pre-2.0.0 state format to current
        newState = {}
        for name, value in state.items():
            nns, index, weightOrBias = name.split('.')
            if nns != 'nns':
                raise Exception('Unexpected state parameter')
            index = int(index)
            if weightOrBias == 'weight':
                for k in range(value.shape[0]):
                    newState['module_%d.lins.%d.weight' % (index, k)] = value[k].t()
            elif weightOrBias == 'bias':
                newState['module_%d.bias' % index] = value
            else:
                raise Exception('Unexpected state parameter')
        state = newState
    model.load_state_dict(state)
