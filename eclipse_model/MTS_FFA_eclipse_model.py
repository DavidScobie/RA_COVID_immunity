import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn
from lmfit import minimize, Parameters, Parameter, report_fit
from scipy.integrate import odeint

############ Read data in and do the basic preprocessing

#import the excel data

df = pd.read_excel('C:\Research_Assistant\work\FW__Human_challenge_studies\COVHIC001_FFA_MTS.xlsx')

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
print('Subj_ID_list',Subj_ID_list)  #fine up to here

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
print('sort_eff_day',sort_eff_day) #this is fine

#remove duplicates in the list of effective days
uni_sort_eff_day = np.unique(sort_eff_day)
print('uni_sort_eff_day',uni_sort_eff_day)

#find if difference between effective day values is greater than 0.5, if so then use this indicie as threshold for later steps
diff_uni_sort_eff_day = np.diff(uni_sort_eff_day)
print('diff_uni_sort_eff_day',diff_uni_sort_eff_day,'len(diff_uni_sort_eff_day)',len(diff_uni_sort_eff_day))

################find length of diff_uni_sort_eff_day and use it to work out whether gap greater than 0.5 happens at start or end of eff_day
start_or_end_thresh = int(np.floor(len(diff_uni_sort_eff_day)/3)) #here we use 3 because peak is roughly a 3rd of the way into data through time
print('start_or_end_thresh',start_or_end_thresh)

not_zerop5 = np.where(diff_uni_sort_eff_day != 0.5) #this is an array of the places where the differences between days is not 0.5
print('not_zerop5',not_zerop5)

count=0 #this is used to count if there is a gap greater than 0.5 in first 3rd so this is factored in if there is also a gap greater than 0.5 in the last 2 thirds in making the dataframe

if len(not_zerop5[0]) == 0: #case where we have data at every half day
    thresh = uni_sort_eff_day[-1]
    #get rid of the rows of dataframe that have effective day above the effective day threshold (as above this day data is discontinuous so harder to plot etc..)
    df_over_4_len_ppl_less_thresh = df_over_4_len_ppl[df_over_4_len_ppl.eff_day <= thresh]
    print('df_over_4_len_ppl_less_thresh',df_over_4_len_ppl_less_thresh)

elif not_zerop5[-1] < start_or_end_thresh: #if there is a gap greater than 0.5 in the first 3rd of the array
    count+=1
    thresh = uni_sort_eff_day[int(not_zerop5[-1])] + 0.5 #need to index last value of array here
    #get rid of the rows of dataframe that have effective day below the effective day threshold (as above this day data is discontinuous so harder to plot etc..)
    df_over_4_len_ppl_less_thresh = df_over_4_len_ppl[df_over_4_len_ppl.eff_day >= thresh]
    print('df_over_4_len_ppl_less_thresh',df_over_4_len_ppl_less_thresh)


if len(not_zerop5[0]) == 0: #case where we have data at every half day
    thresh = uni_sort_eff_day[-1]
    #get rid of the rows of dataframe that have effective day above the effective day threshold (as above this day data is discontinuous so harder to plot etc..)
    df_over_4_len_ppl_less_thresh = df_over_4_len_ppl[df_over_4_len_ppl.eff_day <= thresh]
    print('df_over_4_len_ppl_less_thresh',df_over_4_len_ppl_less_thresh)

elif not_zerop5[0] > start_or_end_thresh: #if there is a gap of greater than 0.5 in the last two thirds of the array then we need to cut the data off after this point
        thresh = uni_sort_eff_day[int(not_zerop5[0])] + 0.5
        if count == 0: #if there is only a gap greater than 0.5 in last two thirds of array
            #get rid of the rows of dataframe that have effective day above the effective day threshold (as above this day data is discontinuous so harder to plot etc..)
            df_over_4_len_ppl_less_thresh = df_over_4_len_ppl[df_over_4_len_ppl.eff_day <= thresh]
            print('df_over_4_len_ppl_less_thresh',df_over_4_len_ppl_less_thresh)
        else: #if there is also a gap greater than 0.5 in first third of array
            #get rid of the rows of dataframe that have effective day above the effective day threshold (as above this day data is discontinuous so harder to plot etc..)
            df_over_4_len_ppl_less_thresh = df_over_4_len_ppl_less_thresh[df_over_4_len_ppl_less_thresh.eff_day <= thresh]
            print('df_over_4_len_ppl_less_thresh',df_over_4_len_ppl_less_thresh)

print('thresh',thresh)
#############################

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

act_div_vir_list_sum_f_e_chop = np.zeros(len(div_vir_list_sum_front_and_end_chopped))
for i in range (len(div_vir_list_sum_front_and_end_chopped)):
    act_div_vir_list_sum_f_e_chop[i] = 10**(div_vir_list_sum_front_and_end_chopped[i])

#######################################################

def f(y, t, paras):
    """
    Your system of differential equations
    """
    #the U, I1 and I2
    U = y[0]
    I2 = y[1]
    I1 = y[2]

    #the parameters alpha, beta, gamma, delta
    try:
        alpha = paras['alpha'].value
        beta = paras['beta'].value
        gamma = paras['gamma'].value
        delta = paras['delta'].value

    except KeyError:
        alpha, beta, gamma, delta = paras

    # the model equations
    f0 = - alpha * U * I2   #dU/dt
    f1 = (gamma * I1) - (delta * I2)    #d(I2)/dt
    f2 = (alpha * U * I2) - (beta * I1) - (gamma * I1)    #d(I1)/dt
    return [f0, f1, f2]

def g(t, x0, paras):
    """
    Solution to the ODE x'(t) = f(t,x,k) with initial condition x(0) = x0
    """
    x = odeint(f, x0, t, args=(paras,))
    return x


def residual(paras, t, data):

    """
    compute the residual between actual data and fitted data
    """

    x0 = paras['U0'].value, paras['I20'].value, paras['I10'].value
    model = g(t, x0, paras)

    # we now have data for the sum of I2 and I1
    I2_I1_model = model[:, 1] + model[:, 2]

    #want to find the residual between the log of the virus measured and fitted data
    log_I2_I1_model = np.log10(I2_I1_model)
    log_data = np.log10(data)

    #so we are now minimising the residual between the data and the sum of I1 and I2
    return (log_I2_I1_model - log_data).ravel()


# initial conditions

#extrapolation to find initial conditions...
def best_fit(X, Y):

    xbar = sum(X)/len(X)
    ybar = sum(Y)/len(Y)
    n = len(X) # or len(Y)

    numer = sum([xi*yi for xi,yi in zip(X, Y)]) - n * xbar * ybar
    denum = sum([xi**2 for xi in X]) - n * xbar**2

    b = numer / denum
    a = ybar - b * xbar

    print('best fit line:\ny = {:.2f} + {:.2f}x'.format(a, b))

    return a, b

# Get the indices of maximum element in eff_day_vals
max_indic_arr = np.where(act_div_vir_list_sum_f_e_chop == np.amax(act_div_vir_list_sum_f_e_chop[:8])) #here we want to extrapolate from the first peak (which is thought to be in the first 8 values)
max_indic = int(max_indic_arr[0])

a, b = best_fit(eff_day_vals_front_and_end_chopped[:max_indic+1],np.log10(act_div_vir_list_sum_f_e_chop)[:max_indic+1])

# plt.figure()
# plt.scatter(eff_day_vals[:max_indic+1], np.log10(act_div_vir_list_sum)[:max_indic+1], marker='o', color='red', label='measured V data', s=75)
# yfit = [a + b * xi for xi in eff_day_vals[:max_indic+1]]
# print('yfit',yfit)
# plt.plot(eff_day_vals[:max_indic+1], yfit)
# plt.xlabel('Days Post Infection')
# plt.ylabel('Virus Titre (Log10 copies/mL)')
# plt.xlim(left=0)
# plt.ylim(bottom=0)

#add the point at time=0, virus=933 to the eff_day_vals and act_div_vir_list_sum arrays
v1 = 0
v2 = 10**a #THIS IS 10**2.97. 2.97 is the y intercept of the line of best fit
eff_day_vals_front_and_end_chopped = np.insert(eff_day_vals_front_and_end_chopped, 0, v1, axis=0)
act_div_vir_list_sum_f_e_chop = np.insert(act_div_vir_list_sum_f_e_chop, 0, v2, axis=0)
print('eff_day_vals_front_and_end_chopped',eff_day_vals_front_and_end_chopped,'act_div_vir_list_sum_f_e_chop',act_div_vir_list_sum_f_e_chop)
"""
#paper initial conditions
U0 = 4*(10**(8))  #the number of cells in an adult is 4x10^8
V0 = 0.31   #cannot be measured as it is below detectable levels. Previous work has shown use <50 copies/ml
I0 = 0   #Should be zero
"""
#my optimised initial conditions
U0 = 6.8*(10**(5))  #the number of cells in an adult is 4x10^8
I20 = act_div_vir_list_sum_f_e_chop[0] / 2  #just taking the first measured value
#I20 = 1
I10 = act_div_vir_list_sum_f_e_chop[0] / 2
#I10 = 1
y0 = [U0, I20, I10]

##cut off the datapoints after day 15 because these are just noise
#only do this if the effective day goes up to 15
exists = 15 in eff_day_vals_front_and_end_chopped
if exists == True:
    fifteenth_indic_arr = np.where(eff_day_vals_front_and_end_chopped == 15)
    fifteenth_indic = int(fifteenth_indic_arr[0])
    print('fifteenth_indic',fifteenth_indic)

    # measured data
    t_measured = eff_day_vals_front_and_end_chopped[:fifteenth_indic]
    V_measured = act_div_vir_list_sum_f_e_chop[:fifteenth_indic]

else:
    # measured data
    t_measured = eff_day_vals_front_and_end_chopped
    V_measured = act_div_vir_list_sum_f_e_chop

#plt.figure()
fig, (ax1, ax2, ax3) = plt.subplots(1,3)
ax1.scatter(t_measured[1:], 10**(-6)*V_measured[1:], marker='o', color='red', label='measured (I2+I1) data', s=75) #the first point is found by extrapolation. Therefore it is not physical so dont plot it.

# set parameters including bounds; you can also fix parameters (use vary=False)
params = Parameters()
params.add('U0', value=U0, vary=False)
params.add('I20', value=I20, vary=False)
params.add('I10', value=I10, vary=False)
"""
#parameters optimised on first 6 days of data
params.add('alpha', value=4.24*(10**(-7)), min=4.23*(10**(-7)), max=4.25*(10**(-7)))   #rate that viral particles infect susceptible cells
params.add('beta', value=61.2, min=61.1, max=61.3)    #Clearance rate of infected cells
params.add('gamma', value=1.83, min=1.82, max=1.84)        #Infected cells release virus at rate gamma
params.add('delta', value=1.45, min=1.44, max=1.46)     #clearance rate of virus particles
"""
#my optimised parameters
params.add('alpha', value=20*(10**(-5)), min=16.9*(10**(-5)), max=45.1*(10**(-5)))   #rate that viral particles infect susceptible cells
params.add('beta', value=50, min=49.9, max=50.1)    #Clearance rate of infected cells
params.add('gamma', value=2.8, min=2.5, max=5.3)        #Infected cells release virus at rate gamma
params.add('delta', value=5, min=4.9, max=5.1)     #clearance rate of virus particles

# fit model
result = minimize(residual, params, args=(t_measured, V_measured), method='leastsq')  # leastsq nelder
# check results of the fit
data_fitted = g(t_measured, y0, result.params)

# plot fitted data
#plt.figure()
ax1.plot(t_measured, 10**(-6)*data_fitted[:, 1], '-', linewidth=2, color='green', label='fitted I2 data')
ax1.plot(t_measured, 10**(-6)*data_fitted[:, 2], '-', linewidth=2, color='blue', label='fitted I1 data')
ax1.plot(t_measured, (10**(-6)*data_fitted[:, 2]) + (10**(-6)*data_fitted[:, 1]), '-', linewidth=2, color='red', label='fitted (I1 + I2) data')
ax1.legend()
ax1.set_xlim([0, max(t_measured)])
#ax1.set_ylim([0, 1.1 * 10**(-6)*max(V_measured)])
ax1.set_xlabel('Days Post Infection')
ax1.set_ylabel('Virus Titre Concentration (million copies/mL)')
ax1.set_title('a)')
# display fitted statistics
report_fit(result)

#plot the fitted data and the model for log(virus) against day
log_V_measured = np.log10(V_measured)
log_Id_fitted = np.log10(data_fitted[:, 1])
log_Is_fitted = np.log10(data_fitted[:, 2])
log_Id_Is_fitted = np.log10(data_fitted[:, 1] + data_fitted[:, 2])
#plt.figure()
ax2.scatter(t_measured[1:], log_V_measured[1:], marker='o', color='red', label='measured (I2+I1) data', s=75) #the first point is found by extrapolation. Therefore it is not physical so dont plot it.
ax2.plot(t_measured, log_Id_fitted, '-', linewidth=2, color='green', label='fitted I2 data')
ax2.plot(t_measured, log_Is_fitted, '-', linewidth=2, color='blue', label='fitted I1 data')
ax2.plot(t_measured, log_Id_Is_fitted, '-', linewidth=2, color='red', label='fitted (I2 + I1) data')

#plotting curve from day -3
first_vir_vals = np.array([(-3*b)+a,a])

first_eff_day_vals = np.array([-3,0])

print('first_vir_vals',first_vir_vals,'first_eff_day_vals',first_eff_day_vals)

ax2.plot(first_eff_day_vals, first_vir_vals, '-', linewidth=2, color='red')
ax2.set_xlim(left=-3)
ax2.legend()
ax2.set_xlabel('Days Post Infection')
ax2.set_ylabel('Virus Titre Concentration (Log10 copies/mL)')
ax2.set_title('b)')

print('virus val day minus 3: ',(-3*b)+a,'virus val day minus 2: ',(-2*b)+a,'virus val day minus 1: ',(-1*b)+a)

#plot the measured data, along with the fitted model for V, I and U
#plt.figure()
I2_fitted = data_fitted[:, 1]
ax3.scatter(t_measured[1:], 10**(-6)*V_measured[1:], marker='o', color='red', label='measured (I2 + I1) data', s=75) #the first point is found by extrapolation. Therefore it is not physical so dont plot it.
ax3.plot(t_measured, 10**(-6)*I2_fitted, '-', linewidth=2, color='green', label='fitted I2 data')
U_fitted = data_fitted[:, 0]
I1_fitted = data_fitted[:, 2]
ax3.plot(t_measured, 10**(-6)*U_fitted, '-', linewidth=2, color='black', label='fitted U data')
ax3.plot(t_measured, 10**(-6)*I1_fitted, '-', linewidth=2, color='blue', label='fitted I1 data')
#plt.ylim(bottom=0.9 * min(log_V_measured))
ax3.set_xlim(left=0)
#ax3.set_ylim([0, 1.1 * 10**(-6)*max(Id_measured)])
ax3.legend()
ax3.set_xlabel('Days Post Infection')
ax3.set_ylabel('Concentration (million copies/mL)')
ax3.set_title('c)')

############## find area under the Is and Id curves
# plt.figure()
# plt.plot(t_measured, data_fitted[:, 1], '-', linewidth=2, color='green', label='fitted Id data')
# print('Id points',data_fitted[:, 1])
# plt.plot(t_measured, data_fitted[:, 2], '-', linewidth=2, color='blue', label='fitted Is data')
# print('Is points',data_fitted[:, 2])
# plt.legend()
# plt.xlabel('Days Post Infection')
# plt.ylabel('Virus Titre (copies/mL)')

Is_area = np.trapz(data_fitted[:, 2], dx=0.5)
Id_area = np.trapz(data_fitted[:, 1], dx=0.5)
print('Is_area',Is_area,'Id_area',Id_area)

#np.save('qPCR_MTS_V_measured_NON_DET_eq_zero_fit_Id+Is', V_measured)
#np.save('qPCR_MTS_t_measured_NON_DET_eq_zero_fit_Id+Is', t_measured)

plt.show()

