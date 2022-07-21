# Using Graph Neural Networks for Site-of-Metabolism Prediction and its Applications to Ranking Promiscuous Enzymatic Products

This repository contains data and code for running the GNN-SOM models for site-of-metabolism prediction.

## Environment Setup

We recommend using [conda](https://docs.conda.io/en/latest/) for managing this environment. Our implementation requires the following software toolkits to be installed:
    * [PyTorch](https://pytorch.org/)
    * [PyTorch Geometric](https://pytorch-geometric.readthedocs.io/)
    * [RDKit](https://www.rdkit.org/)

## Usage

The example code for making SOM predictions on a given enzyme-molecule pair is presented as a Jupyter notebook, which can be found in [GNN-SOM.ipynb](GNN-SOM.ipynb). This notebook makes use of various commonly-used functions provided in the [gnn_som](gnn_som) directory as well as the model state and configuration files in the [data](data) directory.

## License

This project is licensed under the MIT license.
