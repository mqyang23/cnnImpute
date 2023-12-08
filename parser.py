import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description="scRNA-seq data imputation using cnnImpute."
    )
    parser.add_argument("--inputFile", type=str, help="Path to input data.")
    parser.add_argument(
        "-o", "--output",
        type=str,
        default="imputed.csv",
        help="Path to output data counts. Default: ./filename_imputed.csv",
    )

    parser.add_argument(
        "--learning-rate",
        type=float,
        default=0.0001,
        help="Learning rate. Default: 0.0001"
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=32,
        help="Batch size. Default: 64"
    )
    parser.add_argument(
        "--max-epochs",
        type=int,
        default=150,
        help="Maximum number of epochs. Default: 500"
    )
    parser.add_argument(
        "--dropout-rate",
        type=float,
        default=0.3,
        help="Dropout rate for the hidden dropout layer (0<rate<1). Default: 0.2"
    )
    parser.add_argument(
        "--dropoutRate-threshold",
        type=float,
        default=0.5,
        help="Dropout rate threshold for identify missing values. Default: 0.5"
    )
    parser.add_argument(
        "--inputGene-threshold",
        type=float,
        default=0.5,
        help="Threshold of dropout percentage for choosing input genes. Default: 0.5"
    )
    args = parser.parse_args()

    return args