# cnnImpute: Missing Value Recovery for Single Cell RNA Sequencing Data

## 1. Installing

    git clone https://github.com/mqyang23/cnnimpute

## 2. Usage

Place your input dataset into the directory cnnimpute/input, and ensure that the filename is provided in the inputFile parameter.

    from cnnimpute.cnnImpute import cnnImpute
    cnnImpute(inputFile="raw.csv")

## 3. Optional arguments:

| Argument | Description |
| -------- | ----------- |
| output | Sets the name of the output file. By default, the output file is named imputed.csv. |
| learning_rate | Sets the learning rate in the neural network. |
| batch_size | Sets the batch size during training. |
| max_epochs | Sets the number of epochs used in training. |
| dropout_rate | Sets the dropout rate for the dropout layer. |
| dropoutRate_threshold | Sets the threshold of dropout probability. |
| inputGene_threshold | Sets the missing-value percentage threshold for input genes. |

## 4. Output file:

The output file will be located in the directory cnnimpute/output/.
