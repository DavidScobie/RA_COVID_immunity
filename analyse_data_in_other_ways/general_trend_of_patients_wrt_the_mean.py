import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn
from lmfit import minimize, Parameters, Parameter, report_fit
from scipy.integrate import odeint

#import the excel data

df = pd.read_excel ('C:\Research_Assistant\work\data\FW__Human_challenge_studies\COVHIC001_qPCR_TS.xlsx')

#print the column headers
print('colummn headers',list(df.columns.values))

#only take the first 8 columns of dataframe as the rest are NaNs
df2 = df.iloc[:, 0:8]

#get rid of nan
df2 = df2.dropna()

#convert the values with floats into strings
df2['Virus Titre (Log10 copies/mL)'] = df2['Virus Titre (Log10 copies/mL)'].astype(str)

# drop all rows with any 'N/A' and 'N/A ' values in the viral load
df2 = df2[df2['Virus Titre (Log10 copies/mL)'].str.contains("N/A") == False]
df2 = df2[df2['Virus Titre (Log10 copies/mL)'].str.contains("N/A ") == False]
#vir_list_Non_NA = df2['Virus Titre (Log10 copies/mL)'].tolist()

#sort the string type dataframe by length
s=df2['Virus Titre (Log10 copies/mL)'].str.len().sort_values().index
df2_str_sorted = df2.reindex(s)

############ only keep patients who have more than 15 virus values which arent NON DETECTED

#find all the possible subject IDs
Subject_ID = df2_str_sorted['Subject ID'].tolist()
Subject_ID_vals = list(set(Subject_ID))
k=0 #counter for the subject ID's
m=0 #counter for the subjects who have less than 15 NON DETs

Subject_ID_vals_short = Subject_ID_vals[0:3]   #just plotting the first patients as a check up
for j in Subject_ID_vals:

    #print('j',str(j))
    k+=1
    df2_Subj_ID_sorted = df2_str_sorted.loc[df2_str_sorted['Subject ID'] == j] #make a subset of the dataframe based on the subject ID
    #print('df2_Subj_ID_sorted',df2_Subj_ID_sorted)

    #get rid of the DETECTED
    df2_Subj_ID_sorted = df2_Subj_ID_sorted.assign(length = df2_Subj_ID_sorted['Virus Titre (Log10 copies/mL)'].str.len())
    df2_Subj_ID_sorted_just_nums = df2_Subj_ID_sorted #make a copy of dataframe to be filled with just numbers
    df2_Subj_ID_sorted_just_nums = df2_Subj_ID_sorted[df2_Subj_ID_sorted.length < 7] #DETECTED is 8 characters so we remove DETECTED here

    #make NON DETECTED = 0
    df2_Subj_ID_sorted_NON_DET = df2_Subj_ID_sorted  #initialise the dataframe yet to be filled with NON DET values
    df2_Subj_ID_sorted_NON_DET = df2_Subj_ID_sorted[df2_Subj_ID_sorted['Virus Titre (Log10 copies/mL)'].str.len() > 10] #the only lengths greater than 10 are the ones which say 'NON DETECTED'
    df2_Subj_ID_sorted_NON_DET = df2_Subj_ID_sorted_NON_DET.assign(temp = np.zeros(len(df2_Subj_ID_sorted_NON_DET['Virus Titre (Log10 copies/mL)'].tolist())))  #make zeros for NON DETECTED
    df2_Subj_ID_sorted_NON_DET = df2_Subj_ID_sorted_NON_DET.drop('Virus Titre (Log10 copies/mL)', 1)#remove the column of virus
    df2_Subj_ID_sorted_NON_DET = df2_Subj_ID_sorted_NON_DET.rename(columns={"temp": "Virus Titre (Log10 copies/mL)"}) #rename back to original name
    number_NON_DETs = len(df2_Subj_ID_sorted_NON_DET['Virus Titre (Log10 copies/mL)'].values.tolist())

    df2_Subj_ID_sorted_comb = pd.concat([df2_Subj_ID_sorted_just_nums, df2_Subj_ID_sorted_NON_DET], axis=0) #combine the virus numbers and the NON DET zeros into datframe

    print('Subject ID',j,'number of NON DETs',number_NON_DETs,'num time points',df2_Subj_ID_sorted_comb.shape[0])

    #if number of NON DETs less than 15, then cut this out from the bunch. We only want PCR +ve people in the analysis
    if number_NON_DETs < 15: #the case where I want to keep the data
        print('LESS THAN 15 NON DETs')
        m+=1
        #print('m',m)
        if m == 1: #the first time through this loop we just create the dataframe
            print('THE BEGINNING')
            df2_cut_out_many = df2_Subj_ID_sorted_comb
        else: #for all of the subsequent patients with few NON DETs, we append their data to the dataframe
            print('APPENDING')
            df2_cut_out_many = df2_cut_out_many.append(df2_Subj_ID_sorted_comb)
    else:
        print('MORE THAN 15 NON DETs')

print('df2_cut_out_many end',df2_cut_out_many)

#convert the strings to numbers
df2_str_sorted = df2_cut_out_many
df2_str_sorted['Virus Titre (Log10 copies/mL)'] = pd.to_numeric(df2_str_sorted['Virus Titre (Log10 copies/mL)'], downcast="float")
df2_str_sorted['Study Day'] = pd.to_numeric(df2_str_sorted['Study Day'], downcast="float")

###################

##create the 'effective study day' which takes AM and PM into account

#convert study day and Timepoint into lists
day_list = df2_str_sorted['Study Day'].tolist()
#print('length day list',len(day_list),'day_list',day_list)
Timepoint_list = df2_str_sorted['Timepoint'].tolist()
#print('Timepoint_list length',len(Timepoint_list),'Timepoint_list',Timepoint_list)

#create 'effective_day' list and add 0.5 to values of study day if it is PM, (keep same as study day if AM)
effective_day = np.zeros(len(day_list))
for i in range (len(day_list)):
    if Timepoint_list[i] == "AM" or Timepoint_list[i] == "AM ":
        effective_day[i] = day_list[i]
    elif Timepoint_list[i] == "PM" or Timepoint_list[i] == "PM ":
        effective_day[i] = day_list[i] + 0.5
    else:
        print('i',i) #for checking if there is an error

#print('effective_day',effective_day)

#convert the virus numbers to a list
vir_list_Non_DET = df2_str_sorted['Virus Titre (Log10 copies/mL)'].tolist()
#print('vir_list_Non_DET',len(vir_list_Non_DET),'vir_list_Non_DET',vir_list_Non_DET)

#plot the virus against day
"""
plt.figure()
plt.plot(effective_day,vir_list_Non_DET,'bx')
plt.xlabel('Study Day')
plt.ylabel('Virus Titre (Log10 copies/mL)')
"""
#plot the means

#find all the possible effective_day values
eff_day_vals = list(set(effective_day))
eff_day_vals.sort()  #THIS ONLY WORKS IF THERE IS AT LEAST 1 COUNT AT EACH TIME POINT
#print('eff_day_vals',eff_day_vals)

#find the occurences of each of the days
occ=np.zeros(len(eff_day_vals))
for j in effective_day:
    for i in eff_day_vals:
        if i==j:
            occ[int(2*(i-min(eff_day_vals)))]+=1   #0.5 gap between vals, and begin at min val
#print('occ',occ)


#divide virus amount by number of counts on that day
div_vir_list=[]
k=0
for j in effective_day:
    for i in eff_day_vals:
        if i==j:
            div_vir_list.append(vir_list_Non_DET[int(k)]/occ[int(2*(i-min(eff_day_vals)))])
            k+=1
#print('div_vir_list',div_vir_list)


#sum the virus amounts on their specific day
div_vir_list_sum = np.zeros(len(eff_day_vals))
k=0
for j in effective_day:
    for i in eff_day_vals:
        if i==j:
            div_vir_list_sum[int(2*(i-min(eff_day_vals)))]+=div_vir_list[int(k)]
            k+=1
#print('div_vir_list_sum',div_vir_list_sum)

plt.plot(eff_day_vals,div_vir_list_sum,'-rx')
plt.title('Average of all patients')
plt.xlabel('Study Day')
plt.ylabel('Virus Titre (Log10 copies/mL)')
plt.xlim([1,18.5])
plt.ylim([-0.1,10.4])

everyone_average_virus_titre = np.sum(div_vir_list_sum) / len(div_vir_list_sum)

#how many patients do we have? do the patients get sick or stay healthy or both? (out of the 36)

#find all the possible subject IDs
Subject_ID = df2_str_sorted['Subject ID'].tolist()
Subject_ID_vals = list(set(Subject_ID))
#print('Subject_ID',Subject_ID)

#append effective_day to the dataframe
df2_str_sorted['effective_day'] = effective_day

#plot the subjects in different colours
df2_str_sorted['Subject ID'] = df2_str_sorted['Subject ID'].astype(str)

#seaborn.relplot(data=df2_str_sorted, x='effective_day', y='Virus Titre (Log10 copies/mL)', hue='Subject ID')

plt.figure()
seaborn.pointplot(data=df2_str_sorted, x='effective_day', y='Virus Titre (Log10 copies/mL)', hue='Subject ID', ci=None)

pat_average_virus_titre = []

#plot individual patients on different days
Subject_ID_vals_short = Subject_ID_vals[0:3]   #just plotting the first patient as a check up
for j in Subject_ID_vals:
    k+=1
    #plt.figure()
    df2_Subj_ID_sorted = df2_str_sorted[df2_str_sorted['Subject ID'].str.contains(str(j)) == True]  #make a subset of the dataframe based on the subject ID
    df2_Subj_ID_sub_eff_sort = df2_Subj_ID_sorted.sort_values(["effective_day"], ascending=True) #sort the values of the dataframe based on the effective_day
    # df2_Subj_ID_sub_eff_sort.plot(x='effective_day', y='Virus Titre (Log10 copies/mL)',kind='line',xlim=[1,18.5],ylim=[-0.1,10.4]) #plot the subject points as a line plot

    # plt.title('Subject ID=%i' %j)
    # plt.xlabel('Study Day')
    # plt.ylabel('Virus Titre (Log10 copies/mL)')

    # print('subject_ID',j)
    # print('df2_Subj_ID_sub_eff_sort',df2_Subj_ID_sub_eff_sort)

    #collecting data for the histogram
    number_of_virus_vals = df2_Subj_ID_sub_eff_sort.shape[0]
    sum_of_virus_vals = np.sum(df2_Subj_ID_sub_eff_sort['Virus Titre (Log10 copies/mL)'].tolist())
    pat_average_virus_titre.append(sum_of_virus_vals / number_of_virus_vals)

print('pat_average_virus_titre',pat_average_virus_titre)
#histogram of the average of the virus titres
plt.figure()
plt.hist(pat_average_virus_titre, density=False, bins=10)
plt.axvline(x=everyone_average_virus_titre,color='red')
plt.xlabel('Average virus titre over full time course')
plt.ylabel('Number of patients')

#############make plot of the diff = (number of points above average curve) - (number of points below average curve) for each patient

plt.show()

