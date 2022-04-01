import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn
from lmfit import minimize, Parameters, Parameter, report_fit
from scipy.integrate import odeint

############ Read data in and do the basic preprocessing

#import the excel data

df = pd.read_excel('C:\Research_Assistant\work\FW__Human_challenge_studies\COVHIC001_FFA_TS.xlsx')

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
#vir_list_Non_NA = df2['Virus Titre (Log10 copies/mL)'].tolist()

#sort the string type dataframe by length
s=df2['Virus Titre (Log10 FFU/mL)'].str.len().sort_values().index
df_str_sorted = df2.reindex(s)

#get rid of the DETECTED and turn the NON DETECTED into zeros
df_str_sorted['length'] = df_str_sorted['Virus Titre (Log10 FFU/mL)'].str.len()
df_str_sorted_just_nums = df_str_sorted
df_str_sorted_just_nums = df_str_sorted[df_str_sorted.length < 7]

df_str_sorted_NON_DET = df_str_sorted
df_str_sorted_NON_DET = df_str_sorted[df_str_sorted.length > 10]
df_str_sorted_NON_DET['Virus Titre (Log10 FFU/mL)'] = np.zeros(len(df_str_sorted_NON_DET['Virus Titre (Log10 FFU/mL)'].tolist()))

df2_str_sorted = pd.concat([df_str_sorted_just_nums, df_str_sorted_NON_DET], axis=0)

print('df2_str_sorted virus',df2_str_sorted['Virus Titre (Log10 FFU/mL)'].tolist(),'length virus',len(df2_str_sorted['Virus Titre (Log10 FFU/mL)'].tolist()))

#convert the strings to numbers
df2_str_sorted['Virus Titre (Log10 FFU/mL)'] = pd.to_numeric(df2_str_sorted['Virus Titre (Log10 FFU/mL)'], downcast="float")
df2_str_sorted['Study Day'] = pd.to_numeric(df2_str_sorted['Study Day'], downcast="float")

#################################################################   #Only considering patients with more than 4 datapoints

#find all the possible subject IDs
Subject_ID = df2_str_sorted['Subject ID'].tolist()
Subject_ID_vals = list(set(Subject_ID))
#print('Subject_ID',Subject_ID)

#only consider patients who have more than 4 datapoints
vir_list_Non_DET = []
day_list = []
Timepoint_list = []
Subj_ID_list = []
for j in Subject_ID_vals:
    df2_Subj_ID_sorted = df2_str_sorted[df2_str_sorted['Subject ID'].astype(str).str.contains(str(j)) == True]  #make a subset of the dataframe based on the subject ID
    if len(df2_Subj_ID_sorted['Virus Titre (Log10 FFU/mL)'].tolist()) > 4:
        print('j',j) #the subjects who have at least 5 datapoints
        #convert virus numbers to list
        vir_list_Non_DET.append(df2_Subj_ID_sorted['Virus Titre (Log10 FFU/mL)'].tolist())
        day_list.append(df2_Subj_ID_sorted['Study Day'].tolist())
        Timepoint_list.append(df2_Subj_ID_sorted['Timepoint'].tolist())
        Subj_ID_list.append(df2_Subj_ID_sorted['Subject ID'].tolist())

#make all the lists flat
vir_list_Non_DET = [item for sublist in vir_list_Non_DET for item in sublist]
day_list = [item for sublist in day_list for item in sublist]
Timepoint_list = [item for sublist in Timepoint_list for item in sublist]
Subj_ID_list = [item for sublist in Subj_ID_list for item in sublist]
print('vir_list_Non_DET',vir_list_Non_DET)
print('day_list',day_list)
print('Timepoint_list',Timepoint_list)
print('Subj_ID_list',Subj_ID_list)

#create the 'effective study day' which takes AM and PM into account
#create 'effective_day' list and add 0.5 to values of study day if it is PM, (keep same as study day if AM)
effective_day = np.zeros(len(day_list))
for i in range (len(day_list)):
    if Timepoint_list[i] == "AM" or Timepoint_list[i] == "AM ":
        effective_day[i] = day_list[i]
    elif Timepoint_list[i] == "PM" or Timepoint_list[i] == "PM ":
        effective_day[i] = day_list[i] + 0.5
    else:
        print('i',i) #for checking if there is an error

#make dataframe of the data from people who have at least 5 datapoints
df_over_4_len_ppl = pd.DataFrame(list(zip(vir_list_Non_DET, Subj_ID_list,day_list,Timepoint_list,effective_day)), columns =['Virus Titre (Log10 FFU/mL)','Subject ID','Study Day','Timepoint','eff_day'])

##we need at least one point at each half day, cut the data off above this threshold where we dont have the data

#sort effective day values into ascending order
sort_eff_day = np.sort(effective_day)

#remove duplicates in the list of effective days
uni_sort_eff_day = np.unique(sort_eff_day)

#find if difference between effective day values is greater than 0.5, if so then use this indicie as threshold for later steps
diff_uni_sort_eff_day = np.diff(uni_sort_eff_day)



not_zerop5 = np.where(diff_uni_sort_eff_day != 0.5)[0]
print('not_zerop5',not_zerop5)

#if difference is 0.5 between all the effective days then thresh is the last possible day. If not then choose day accordingly
if len(not_zerop5) == 0: #case where we have data at every half day
    thresh = uni_sort_eff_day[-1]
else: #case where we dont have data at every half day
    thresh = uni_sort_eff_day[int(not_zerop5[0])] + 0.5

#get rid of the rows of dataframe that have effective day above the effective day threshold (as above this day data is discontinuous so harder to plot etc..)
df_over_4_len_ppl_less_thresh = df_over_4_len_ppl[df_over_4_len_ppl.eff_day <= thresh]
print('df_over_4_len_ppl_less_thresh',df_over_4_len_ppl_less_thresh)

##plotting the patients in different colours

plt.figure()
seaborn.pointplot(data=df_over_4_len_ppl_less_thresh, x='eff_day', y='Virus Titre (Log10 FFU/mL)', hue='Subject ID', ci=None)
plt.title('Plotting individual patients with more than 4 datapoints in different colours')

#convert study day and Timepoint into lists
day_list = df_over_4_len_ppl_less_thresh['Study Day'].tolist()
print('length day list',len(day_list),'day_list',day_list)
Timepoint_list = df_over_4_len_ppl_less_thresh['Timepoint'].tolist()
print('Timepoint_list length',len(Timepoint_list),'Timepoint_list',Timepoint_list)

#convert the virus numbers to a list
vir_list_Non_DET = df_over_4_len_ppl_less_thresh['Virus Titre (Log10 FFU/mL)'].tolist()
print('vir_list_Non_DET',len(vir_list_Non_DET),'vir_list_Non_DET',vir_list_Non_DET)

effective_day = df_over_4_len_ppl_less_thresh['eff_day'].tolist()
print('effective_day',len(effective_day),'effective_day',effective_day)
"""
#plot the virus against day
# plt.figure()
# plt.plot(effective_day,vir_list_Non_DET,'bx')
# plt.xlabel('Study Day')
# plt.ylabel('Virus Titre (Log10 FFU/mL)')
"""
##plot the means

#find all the possible effective_day values
eff_day_vals = list(set(effective_day))
eff_day_vals.sort()  #THIS ONLY WORKS IF THERE IS AT LEAST 1 COUNT AT EACH TIME POINT
print('eff_day_vals',eff_day_vals)
print('eff_day_vals',eff_day_vals,'length eff_day_vals',len(eff_day_vals))

#find the occurences of each of the days
occ=np.zeros(len(eff_day_vals))
for j in effective_day:
    for i in eff_day_vals:
        if i==j:
            occ[int(2*(i-min(eff_day_vals)))]+=1   #0.5 gap between vals, and begin at min val
print('occ',occ,'length occ',len(occ))


#divide virus amount by number of counts on that day
div_vir_list=[]
k=0
for j in effective_day:
    for i in eff_day_vals:
        if i==j:
            div_vir_list.append(vir_list_Non_DET[int(k)]/occ[int(2*(i-min(eff_day_vals)))])
            k+=1
print('div_vir_list',div_vir_list,'len div_vir_list',len(div_vir_list))     #CORRECT


#sum the virus amounts on their specific day
div_vir_list_sum = np.zeros(len(eff_day_vals))
k=0
for j in effective_day:
    for i in eff_day_vals:
        if i==j:
            div_vir_list_sum[int(2*(i-min(eff_day_vals)))]+=div_vir_list[int(k)]
            k+=1
print('div_vir_list_sum',div_vir_list_sum,'length div_vir_list_sum',len(div_vir_list_sum))

##if the mean value is zero, then discard it from the array
print('div_vir_list_sum before',div_vir_list_sum)
div_vir_list_sum_end_chopped = np.trim_zeros(div_vir_list_sum, trim='b')
print('div_vir_list_sum after',div_vir_list_sum_end_chopped)

#chop the zeros off the end of the days and virus arrays (as at the end we just have zeros from NON DETECTED)
len_with_zeros = len(div_vir_list_sum)
len_without_zeros = len(div_vir_list_sum_end_chopped)
eff_day_vals_end_chopped = eff_day_vals[:-(len_with_zeros - len_without_zeros)]
print('eff_day_vals_end_chopped',eff_day_vals_end_chopped,'length',len(eff_day_vals_end_chopped),'eff_day_vals',eff_day_vals,'length',len(eff_day_vals),'len with zeros',len_with_zeros,'len without zeros',len_without_zeros)

#chop the zeros off the start of days and virus days (as at the start we just have zeros from NON DETECTED)
div_vir_list_sum_front_and_end_chopped = np.trim_zeros(div_vir_list_sum_end_chopped, trim='f')
difference =  len(div_vir_list_sum_end_chopped) - len(div_vir_list_sum_front_and_end_chopped)
print('difference',difference)
eff_day_vals_front_and_end_chopped = eff_day_vals_end_chopped[difference:]
print('eff_day_vals_front_and_end_chopped',eff_day_vals_front_and_end_chopped)

plt.figure()
plt.plot(eff_day_vals_front_and_end_chopped,div_vir_list_sum_front_and_end_chopped,'-rx')  #something is going on with these indicies here, find out what
plt.xlabel('Days Post Infection')
plt.ylabel('Virus Titre (Log10 FFU/mL)')
plt.title('patients with at least 5 datapoints log scale')

#plot actual virus amount (instead of log10 of virus amount)
act_div_vir_list_sum = np.zeros(len(div_vir_list_sum))
for i in range (len(div_vir_list_sum)):
    act_div_vir_list_sum[i] = 10**(div_vir_list_sum[i])

plt.figure()
plt.plot(eff_day_vals,act_div_vir_list_sum,'-rx')
plt.xlabel('Days Post Infection')
plt.ylabel('Virus Titre (copies/mL)')
plt.title('patients with at least 5 datapoints linear scale')

#np.save('FFA_TS_V_measured_NON_DET_eq_zero', 10**div_vir_list_sum_front_and_end_chopped)
#np.save('FFA_TS_t_measured_NON_DET_eq_zero', eff_day_vals_front_and_end_chopped)

plt.show()