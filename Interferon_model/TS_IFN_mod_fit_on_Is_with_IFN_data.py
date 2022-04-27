import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn
from lmfit import minimize, Parameters, Parameter, report_fit
from scipy.integrate import odeint

#import the excel data

df = pd.read_excel(r'C:\Research_Assistant\work\data\interferon\Nasal_mediators_to_sd4_for_MN.xlsx')  #r makes it a raw string

#print the column headers
print('colummn headers',list(df.columns.values))

#convert times to a list
time = df['Time'].tolist()
print('time',time)

#convert Donor (SUBJECT ID) to list
Subject_ID = df['Donor'].tolist()
print('Subject_ID',Subject_ID)

seaborn.relplot(data=df, x='Time', y='IFN-α', hue='Donor')

#find all the possible subject IDs
Subject_ID_vals = list(set(Subject_ID))

#plot individual patients on different days
k=0
Subject_ID_vals_short = Subject_ID_vals[0:3]   #just plotting the first patient as a check up
for j in Subject_ID_vals:
    k+=1
    #plt.figure()
    df_Subj_ID_sorted = df.loc[df['Donor'] == j]  #make a subset of the dataframe based on the subject ID
    df_Subj_ID_sub_eff_sort = df_Subj_ID_sorted.sort_values(["Time"], ascending=True) #sort the values of the dataframe based on the time
    # df_Subj_ID_sub_eff_sort.plot(x='Time', y='IFN-α',kind='line') #plot the subject points as a line plot

    # plt.title('Subject ID=%i' %j)
    # plt.xlabel('Study Day')
    # plt.ylabel('IFN-α')

#########plot the means

#convert the IFN numbers to a list
IFN_alpha_list = df['IFN-α'].tolist()

#find all the possible effective_day values
eff_day_vals = list(set(time))
eff_day_vals.sort()  #THIS ONLY WORKS IF THERE IS AT LEAST 1 COUNT AT EACH TIME POINT
#print('eff_day_vals',eff_day_vals)

#find the occurences of each of the days
occ=np.zeros(len(eff_day_vals))
for j in time:
    for i in eff_day_vals:
        if i==j:
            occ[int(2*(i-min(eff_day_vals)))]+=1   #0.5 gap between vals, and begin at min val
#print('occ',occ)


#divide virus amount by number of counts on that day
div_IFN_alpha_list=[]
k=0
for j in time:
    for i in eff_day_vals:
        if i==j:
            div_IFN_alpha_list.append(IFN_alpha_list[int(k)]/occ[int(2*(i-min(eff_day_vals)))])
            k+=1
#print('div_vir_list',div_vir_list)


#sum the virus amounts on their specific day
div_IFN_alpha_list_sum = np.zeros(len(eff_day_vals))
k=0
for j in time:
    for i in eff_day_vals:
        if i==j:
            div_IFN_alpha_list_sum[int(2*(i-min(eff_day_vals)))]+=div_IFN_alpha_list[int(k)]
            k+=1
#print('div_vir_list_sum',div_vir_list_sum)

plt.figure()
plt.plot(eff_day_vals,div_IFN_alpha_list_sum,'-rx')
plt.xlabel('Time/days')
plt.ylabel('IFN-α')

##########################
#Do the same process but for the IFN beta

seaborn.relplot(data=df, x='Time', y='IFN-β', hue='Donor')

#convert the IFN numbers to a list
IFN_beta_list = df['IFN-β'].tolist()

#find all the possible effective_day values
eff_day_vals = list(set(time))
eff_day_vals.sort()  #THIS ONLY WORKS IF THERE IS AT LEAST 1 COUNT AT EACH TIME POINT
#print('eff_day_vals',eff_day_vals)

#find the occurences of each of the days
occ=np.zeros(len(eff_day_vals))
for j in time:
    for i in eff_day_vals:
        if i==j:
            occ[int(2*(i-min(eff_day_vals)))]+=1   #0.5 gap between vals, and begin at min val
#print('occ',occ)


#divide virus amount by number of counts on that day
div_IFN_beta_list=[]
k=0
for j in time:
    for i in eff_day_vals:
        if i==j:
            div_IFN_beta_list.append(IFN_beta_list[int(k)]/occ[int(2*(i-min(eff_day_vals)))])
            k+=1
#print('div_vir_list',div_vir_list)


#sum the virus amounts on their specific day
div_IFN_beta_list_sum = np.zeros(len(eff_day_vals))
k=0
for j in time:
    for i in eff_day_vals:
        if i==j:
            div_IFN_beta_list_sum[int(2*(i-min(eff_day_vals)))]+=div_IFN_beta_list[int(k)]
            k+=1
#print('div_vir_list_sum',div_vir_list_sum)

plt.figure()
plt.plot(eff_day_vals,div_IFN_beta_list_sum,'-rx')
plt.xlabel('Time/days')
plt.ylabel('IFN-α')

#plot the sum of (IFN alpha + IFN beta)
plt.figure()
plt.plot(eff_day_vals,div_IFN_beta_list_sum + div_IFN_alpha_list_sum,'-rx')
plt.xlabel('Time/days')
plt.ylabel('IFN-α + IFN-β')

#plot the log of the sum of (IFN alpha + IFN beta)
plt.figure()
plt.plot(eff_day_vals,np.log10(div_IFN_beta_list_sum + div_IFN_alpha_list_sum),'-rx')
plt.xlabel('Time/days')
plt.ylabel('log(IFN-α + IFN-β)')

#plot the log of IFN alpha
plt.figure()
plt.plot(eff_day_vals,np.log10(div_IFN_alpha_list_sum),'-rx')
plt.xlabel('Time/days')
plt.ylabel('log(IFN-α)')

#plot the log of IFN beta
plt.figure()
plt.plot(eff_day_vals,np.log10(div_IFN_beta_list_sum),'-rx')
plt.xlabel('Time/days')
plt.ylabel('log(IFN-β)')

plt.show()

