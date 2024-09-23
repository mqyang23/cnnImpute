# cnnImpute: missing value recovery for single cell RNA sequencing data


Wenjuan Zhang, Brandon Huckaby, John Talburt, Sherman Weissman, Mary Qu Yang.
"cnnImpute: missing value recovery for single cell RNA sequencing data. Scientific Reports, 14(1), 3946"
https://www.nature.com/articles/s41598-024-53998-x

DeepImpute has been implemented in Python3 and R. 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine.

### Installing

You can install cnnImpute's latest release using pip with the following command:


```bash
git clone https://github.com/mqyang23/cnnImpute
```

## Copy conda environment
If you haven't installed conda, please install it first. After that, run this code to copy my conda environment:
```bash
conda env create -f environment_cnnimpute.yml
```

## install R dependencies
First, activate conda environment py123:
```bash
conda activate py123
```


### Usage

cnnImpute can be used as a Python package.

```python
from cnnImpute.cnnimpute import cnnImpute
cnnImpute(inputFile="/home/user/raw.csv",output="/home/user/cnnImputed.csv")
```



```
Positional parameters:
  inputFile             Absolute path to input data.
  output                Absolute path with designed name to output data.

Parameters:
  inputFile                  Path to input data. 
  output                     Path to output data. Default: cnnImpute/imputed.csv
  learning-rate              Learning rate. Default: 0.0001
  batch-size                 Batch size. Default: 32
  max-epochs                 Maz epochs while training the model. Default: 150.
  drop-rate                  Dropout rate for the hidden dropout layer. Default: 0.3
  dropoutRate-threshold      Dropout tate threshold for identify missing values. Default: 0.5
  inputGene-threshold        Threshold of dropout percentage for choosing input genes. Default: 0.5
```

