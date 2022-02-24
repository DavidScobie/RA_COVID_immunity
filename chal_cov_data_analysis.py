import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn

#import the excel data
df = pd.read_excel ('C:\RA_COVID_immunity\work\FW__Human_challenge_studies\COVHIC001_FFA_MTS.xlsx')

#print the column headers
print('colummn headers',list(df.columns.values))

#only take the first 8 columns of dataframe as the rest is NaNs
df1 = df.iloc[:, 0:8]
#print ('df1',df1)

# drop all rows with any NaN values in the viral load
df2 = df1.dropna()
#print('df2',df2)

#print the data types in the dat frame
#print('df2 dtypes',df2.dtypes)

#just the floats
# df2_floats = df2[df2['Virus Titre (Log10 FFU/mL)'].apply(lambda x: isinstance(x,float))]
# df2_floats_list = df2_floats['Virus Titre (Log10 FFU/mL)'].tolist()
# print('df2_floats_list',df2_floats_list)

#convert the values with floats into strings
df['Virus Titre (Log10 FFU/mL)'] = df['Virus Titre (Log10 FFU/mL)'].astype(str)
df2_floats_to_str_list = df2['Virus Titre (Log10 FFU/mL)'].tolist()
#print('df2_floats_to_str_list',df2_floats_to_str_list)

#just the strings
#df2_str = df2[df2['Virus Titre (Log10 FFU/mL)'].apply(lambda x: isinstance(x,str))]

#sort the string type dataframe by length
s=df2['Virus Titre (Log10 FFU/mL)'].str.len().sort_values().index
df2_str_sorted = df2.reindex(s)
#print('sorted',df2_str_sorted)

#print the list of virus amounts
vir_list = df2_str_sorted['Virus Titre (Log10 FFU/mL)'].tolist()
#print('vir_list',vir_list)

#next get rid of the NaN's, then chop off the DETECTED and NOT DETECTED, combine str and float dataframes, and then plot against days

#get rid of the 'N/A '
df2_str_sorted = df2_str_sorted[df2_str_sorted['Virus Titre (Log10 FFU/mL)'].str.contains("N/A ") == False]
vir_list_Non_NA = df2_str_sorted['Virus Titre (Log10 FFU/mL)'].tolist()
#print('vir_list_Non_NA',vir_list_Non_NA)

#get rid of the DETECTED AND NOT DETECTED
df2_str_sorted['length'] = df2_str_sorted['Virus Titre (Log10 FFU/mL)'].str.len()
df2_str_sorted = df2_str_sorted[df2_str_sorted.length < 7]

#convert the strings to numbers
df2_str_sorted['Virus Titre (Log10 FFU/mL)'] = pd.to_numeric(df2_str_sorted['Virus Titre (Log10 FFU/mL)'], downcast="float")
df2_str_sorted['Study Day'] = pd.to_numeric(df2_str_sorted['Study Day'], downcast="float")

#sort the dataframe by day
#df2_str_sorted = df2_str_sorted.sort_values(['Study Day'], ascending=True)

#convert the virus numbers to a list
vir_list_Non_DET = df2_str_sorted['Virus Titre (Log10 FFU/mL)'].tolist()
print('vir_list_Non_DET',vir_list_Non_DET)

#convert the day numbers to a list
day_list = df2_str_sorted['Study Day'].tolist()
print('length day list',len(day_list),'day_list',day_list)

#plot the virus against day
#df2_str_sorted.plot(x='Study Day', y='Virus Titre (Log10 FFU/mL)', style='o')
plt.figure(0)
plt.plot(day_list,vir_list_Non_DET,'bx')
plt.xlabel('Study Day')
plt.ylabel('Virus Titre (Log10 FFU/mL)')

#plot the means

#find all the possible day values
all_day_vals = list(set(day_list))
print('all_day_vals',all_day_vals)

#find the occurences of each of the days
occ=np.zeros(len(all_day_vals))
for j in day_list:
    for i in all_day_vals:
        if i==j:
            occ[int(i-min(all_day_vals))]+=1
print('occ',occ)

#divide virus amount by number of counts on that day
div_vir_list=[]
k=0
for j in day_list:
    for i in all_day_vals:
        if i==j:
            div_vir_list.append(vir_list_Non_DET[int(k)]/occ[int(i-min(all_day_vals))])
            k+=1

#sum the virus amounts on their specific day
div_vir_list_sum = np.zeros(len(all_day_vals))
k=0
for j in day_list:
    for i in all_day_vals:
        if i==j:
            div_vir_list_sum[int(i-min(all_day_vals))]+=div_vir_list[int(k)]
            k+=1
print('div_vir_list_sum',div_vir_list_sum)

plt.plot(all_day_vals,div_vir_list_sum,'-rx')

#how many patients do we have? do the patients get sick or stay healthy or both? (out of the 36)

#find all the possible subject IDs
Subject_ID = df2_str_sorted['Subject ID'].tolist()
Subject_ID_vals = list(set(Subject_ID))
print('Subject_ID',Subject_ID)

#plot the subjects in different colours
df2_str_sorted['Subject ID'] = df2_str_sorted['Subject ID'].astype(str)
seaborn.relplot(data=df2_str_sorted, x='Study Day', y='Virus Titre (Log10 FFU/mL)', hue='Subject ID')
plt.figure(2)
seaborn.pointplot(data=df2_str_sorted, x='Study Day', y='Virus Titre (Log10 FFU/mL)', hue='Subject ID', ci=None)

plt.show()














"""

#print the column for the day
day = df2[['Study Day']]
print('day',day)

#print the column for the viral load
virus = df2['Virus Titre (Log10 FFU/mL)']
print('virus',virus)

#only take the virus values which are numbers
#virus_num = virus.select_dtypes(include='float')
#virus_num = virus.select_dtypes(include=np.number)
#virus_num = virus._get_numeric_data()
#print('virus_num',virus_num)

df3 = virus.astype(str, errors = 'raise')
print('df3',df3)

print('virus dtypes',df3.dtypes)

print(df2.sort_values(by='Virus Titre (Log10 FFU/mL)') )

#plot the virus titre per day
# df.plot(x='Study Day', y='Virus Titre (Log10 FFU/mL)')
# plt.show()
"""