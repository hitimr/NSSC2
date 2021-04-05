print("hello")

import pandas as pd
import numpy


dataframe=pd.read_csv('testlog.csv')

print(dataframe.head((0)))


for i in range(len(dataframe.head(0))):
    print(dataframe.head(i))





#python3 plots.py
