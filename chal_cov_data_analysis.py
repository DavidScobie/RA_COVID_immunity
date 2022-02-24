import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn

#import the excel data
df = pd.read_excel ('C:\Research_Assistant\work\FW__Human_challenge_studies\COVHIC001_FFA_MTS.xlsx')

#print the column headers
print('colummn headers',list(df.columns.values))

#only take the first 8 columns of dataframe as the rest are NaNs
df2 = df.iloc[:, 0:8]

#get rid of nan
df2 = df2.dropna()

#convert the values with floats into strings
df2['Virus Titre (Log10 FFU/mL)'] = df2['Virus Titre (Log10 FFU/mL)'].astype(str)

# drop all rows with any 'N/A' and 'N/A ' values in the viral load
df2 = df2[df2['Virus Titre (Log10 FFU/mL)'].str.contains("N/A") == False]
df2 = df2[df2['Virus Titre (Log10 FFU/mL)'].str.contains("N/A ") == False]
#vir_list_Non_NA = df2['Virus Titre (Log10 FFU/mL)'].tolist()

#sort the string type dataframe by length
s=df2['Virus Titre (Log10 FFU/mL)'].str.len().sort_values().index
df2_str_sorted = df2.reindex(s)

#get rid of the DETECTED AND NOT DETECTED
df2_str_sorted['length'] = df2_str_sorted['Virus Titre (Log10 FFU/mL)'].str.len()
df2_str_sorted = df2_str_sorted[df2_str_sorted.length < 7]

#convert the strings to numbers
df2_str_sorted['Virus Titre (Log10 FFU/mL)'] = pd.to_numeric(df2_str_sorted['Virus Titre (Log10 FFU/mL)'], downcast="float")
df2_str_sorted['Study Day'] = pd.to_numeric(df2_str_sorted['Study Day'], downcast="float")

#print all the virus titre to 2dp
vir_list = df2_str_sorted['Virus Titre (Log10 FFU/mL)'].tolist()
round_to_tenths = [round(num, 2) for num in vir_list]
round_to_tenths.sort()
print('length vir_list',len(vir_list),'vir_list',round_to_tenths)

##create the 'effective study day' which takes AM and PM into account

#convert study day and Timepoint into lists
day_list = df2_str_sorted['Study Day'].tolist()
print('length day list',len(day_list),'day_list',day_list)
Timepoint_list = df2_str_sorted['Timepoint'].tolist()
print('Timepoint_list length',len(Timepoint_list),'Timepoint_list',Timepoint_list)

#create 'effective day' list and add 0.5 to values of study day if it is PM, (keep same as study day if AM)
effective_day = np.zeros(len(day_list))
for i in range (len(day_list)):
    if Timepoint_list[i] == "AM" or Timepoint_list[i] == "AM ":
        effective_day[i] = day_list[i]
    elif Timepoint_list[i] == "PM" or Timepoint_list[i] == "PM ":
        effective_day[i] = day_list[i] + 0.5
    else:
        print('i',i) #for checking if there is an error

print('effective_day',effective_day)


#convert the virus numbers to a list
vir_list_Non_DET = df2_str_sorted['Virus Titre (Log10 FFU/mL)'].tolist()
print('vir_list_Non_DET',len(vir_list_Non_DET),'vir_list_Non_DET',vir_list_Non_DET)

#convert the day numbers to a list
# day_list = df2_str_sorted['Study Day'].tolist()
# print('length day list',len(day_list),'day_list',day_list)

#plot the virus against day
#df2_str_sorted.plot(x='Study Day', y='Virus Titre (Log10 FFU/mL)', style='o')
plt.figure(0)
plt.plot(effective_day,vir_list_Non_DET,'bx')
plt.xlabel('Study Day')
plt.ylabel('Virus Titre (Log10 FFU/mL)')

"""

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

#plot individual patients on different days

"""

plt.show()

