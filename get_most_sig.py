#!/bin/env/python3

import numpy as np
import pandas as pd
import sys

source_filename = sys.argv[1]

df = pd.read_csv(source_filename, compression='gzip', sep='\t')
print(df)