import sys
import miner2 
import miner2.preprocess
import miner2.mechanistic_inference
import miner2.miner
import pandas as pd
import os
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import dill
import pickle
import seaborn as sb
import time

dill.load_session('./output/' + '2019-07-30-17:52:28.dill')
print(alpha)