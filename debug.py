import miner2
import miner2.coexpression
import miner2.preprocess
import miner2.mechanistic_inference
import pandas as pd
import os

EXPRESSION_DATA_FILE = '/home/aliu/omics4tb2/aliu/projects/causalAssociation/results/expected/GSM3587977_AML707B-D97.dem.txt'
expression_data, conversion_table = miner2.miner2.preprocess.main(EXPRESSION_DATA_FILE)
