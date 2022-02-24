import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#import the excel data
df = pd.read_excel ('C:\RA_COVID_immunity\work\FW__Human_challenge_studies\COVHIC001_qPCR_TS.xlsx')

#print the column headers
print('colummn headers',list(df.columns.values))

#find all the possible subject IDs
Subject_ID = df['Subject ID'].tolist()
Subject_ID_vals = list(set(Subject_ID))
print('Subject_ID_vals',Subject_ID_vals,'length Subject_ID_vals',len(Subject_ID_vals))