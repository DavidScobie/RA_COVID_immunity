import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn
from lmfit import minimize, Parameters, Parameter, report_fit
from scipy.integrate import odeint

#import the excel data

df = pd.read_excel ('C:\Research_Assistant\work\FW__Human_challenge_studies\COVHIC001_qPCR_MTS.xlsx')

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
df_str_sorted = df2.reindex(s)

#get rid of the DETECTED AND NOT DETECTED
df_str_sorted['length'] = df_str_sorted['Virus Titre (Log10 copies/mL)'].str.len()
df_str_sorted_just_nums = df_str_sorted
df_str_sorted_just_nums = df_str_sorted[df_str_sorted.length < 7]

df_str_sorted_NON_DET = df_str_sorted
df_str_sorted_NON_DET = df_str_sorted[df_str_sorted.length > 10]
df_str_sorted_NON_DET['Virus Titre (Log10 copies/mL)'] = np.zeros(len(df_str_sorted_NON_DET['Virus Titre (Log10 copies/mL)'].tolist()))

df2_str_sorted = pd.concat([df_str_sorted_just_nums, df_str_sorted_NON_DET], axis=0)

print('df2_str_sorted virus',df2_str_sorted['Virus Titre (Log10 copies/mL)'].tolist(),'length virus',len(df2_str_sorted['Virus Titre (Log10 copies/mL)'].tolist()))

#convert the strings to numbers
df2_str_sorted['Virus Titre (Log10 copies/mL)'] = pd.to_numeric(df2_str_sorted['Virus Titre (Log10 copies/mL)'], downcast="float")
df2_str_sorted['Study Day'] = pd.to_numeric(df2_str_sorted['Study Day'], downcast="float")

#print all the virus titre to 2dp
vir_list = df2_str_sorted['Virus Titre (Log10 copies/mL)'].tolist()
round_to_tenths = [round(num, 2) for num in vir_list]
round_to_tenths.sort()
#print('length vir_list',len(vir_list),'vir_list',round_to_tenths)

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
"""
plt.plot(eff_day_vals,div_vir_list_sum,'-rx')
"""

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
"""
plt.figure()
seaborn.pointplot(data=df2_str_sorted, x='effective_day', y='Virus Titre (Log10 copies/mL)', hue='Subject ID', ci=None)

#plot individual patients on different days
Subject_ID_vals_short = Subject_ID_vals[0:3]   #just plotting the first patient as a check up
for j in Subject_ID_vals:
    k+=1
    #plt.figure()
    df2_Subj_ID_sorted = df2_str_sorted[df2_str_sorted['Subject ID'].str.contains(str(j)) == True]  #make a subset of the dataframe based on the subject ID
    df2_Subj_ID_sub_eff_sort = df2_Subj_ID_sorted.sort_values(["effective_day"], ascending=True) #sort the values of the dataframe based on the effective_day
    df2_Subj_ID_sub_eff_sort.plot(x='effective_day', y='Virus Titre (Log10 copies/mL)',kind='line',xlim=[1,18.5],ylim=[2.8,10.4]) #plot the subject points as a line plot

    plt.title('Subject ID=%i' %j)
    plt.xlabel('Study Day')
    plt.ylabel('Virus Titre (Log10 copies/mL)')
"""


#plot actual virus amount (instead of log10 of virus amount)
act_div_vir_list_sum = np.zeros(len(div_vir_list_sum))
for i in range (len(div_vir_list_sum)):
    act_div_vir_list_sum[i] = 10**(div_vir_list_sum[i])

plt.figure()
plt.plot(eff_day_vals,div_vir_list_sum,'-rx')
plt.xlabel('Days Post Infection')
plt.ylabel('Virus Titre (copies/mL)')

#######################################################

def f(y, t, paras):
    """
    Your system of differential equations
    """
    #the U, I1, I2, V and F
    U = y[0]
    I1 = y[1]
    I2 = y[2]
    V = y[3]
    F = y[4]

    #the parameters alpha, beta, gamma, delta, mu, s, tau and epsilon
    try:
        alpha = paras['alpha'].value
        beta = paras['beta'].value
        gamma = paras['gamma'].value
        delta = paras['delta'].value
        mu = paras['mu'].value
        s = paras['s'].value
        tau = paras['tau'].value
        epsilon = paras['epsilon'].value

    except KeyError:
        alpha, beta, gamma, delta, mu, s, tau, epsilon = paras

    # the model equations
    f0 = - alpha * V * U   #dU/dt
    f1 = (alpha * U * V) - (beta * I1)    #d(I1)/dt
    f2 = (beta * I1) - (mu * I2)    #d(I2)/dt
    f3 = (gamma * I2) - (delta * V) #d(V)/dt
    f4 = (s * I2 * (t - tau)) - (epsilon * F) #d(F)/dt
    return [f0, f1, f2, f3, f4]

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

    x0 = paras['U0'].value, paras['I10'].value, paras['I20'].value, paras['V0'].value, paras['F0'].value
    model = g(t, x0, paras)

    # we have data for one of our variables
    V_model = model[:, 3]

    #want to find the residual between the log of the virus measured and fitted data
    log_V_model = np.log10(V_model)
    log_data = np.log10(data)

    #so we are now minimising the residual between the data and V
    return (log_V_model - log_data).ravel()


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
max_indic_arr = np.where(act_div_vir_list_sum == np.amax(act_div_vir_list_sum[:8])) #here we want to extrapolate from the first peak (which is thought to be in the first 8 values)
max_indic = int(max_indic_arr[0])

a, b = best_fit(eff_day_vals[:max_indic+1],np.log10(act_div_vir_list_sum)[:max_indic+1])

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
eff_day_vals = np.insert(eff_day_vals, 0, v1, axis=0)
act_div_vir_list_sum = np.insert(act_div_vir_list_sum, 0, v2, axis=0)
print('eff_day_vals',eff_day_vals,'act_div_vir_list_sum',act_div_vir_list_sum)
"""
#paper initial conditions
U0 = 4*(10**(8))  #the number of cells in an adult is 4x10^8
V0 = 0.31   #cannot be measured as it is below detectable levels. Previous work has shown use <50 copies/ml
I0 = 0   #Should be zero
"""
#my optimised initial conditions
U0 = 4*(10**(6))  #the number of cells in an adult is 4x10^8
#I10 = act_div_vir_list_sum_f_e_chop[0] / 2
I10 = 0
#I20 = act_div_vir_list_sum_f_e_chop[0] / 2  #just taking the first measured value
I20 = 0
V0 = act_div_vir_list_sum[0]
F0 = 0
y0 = [U0, I10, I20, V0, F0]

##cut off the datapoints after day 15 because these are just noise
#only do this if the effective day goes up to 15
exists = 15 in eff_day_vals
if exists == True:
    fifteenth_indic_arr = np.where(np.array(eff_day_vals) == 15)

    print('fifteenth_indic_arr[0]',fifteenth_indic_arr[0],'fifteenth_indic_arr',fifteenth_indic_arr,'eff_day_vals',eff_day_vals)

    fifteenth_indic = int(fifteenth_indic_arr[0])
    print('fifteenth_indic',fifteenth_indic)

    # measured data
    t_measured = eff_day_vals[:fifteenth_indic]
    V_measured = act_div_vir_list_sum[:fifteenth_indic]

else:
    # measured data
    t_measured = eff_day_vals
    V_measured = act_div_vir_list_sum

#plt.figure()
fig, (ax1, ax2, ax3) = plt.subplots(1,3)
ax1.scatter(t_measured[1:], 10**(-6)*V_measured[1:], marker='o', color='red', label='measured V data', s=75) #the first point is found by extrapolation. Therefore it is not physical so dont plot it.

# set parameters including bounds; you can also fix parameters (use vary=False)
params = Parameters()
params.add('U0', value=U0, vary=False)
params.add('I10', value=I10, vary=False)
params.add('I20', value=I20, vary=False)
params.add('V0', value=V0, vary=False)
params.add('F0', value=F0, vary=False)
"""
#parameters optimised on first 6 days of data
params.add('alpha', value=4.24*(10**(-7)), min=4.23*(10**(-7)), max=4.25*(10**(-7)))   #rate that viral particles infect susceptible cells
params.add('beta', value=61.2, min=61.1, max=61.3)    #Clearance rate of infected cells
params.add('gamma', value=1.83, min=1.82, max=1.84)        #Infected cells release virus at rate gamma
params.add('delta', value=1.45, min=1.44, max=1.46)     #clearance rate of virus particles
"""
#my optimised parameters
params.add('alpha', value=4*(10**(-4)), min=5*(10**(-5)), max=9*(10**(-4)))   #rate that viral particles infect susceptible cells
params.add('beta', value=0.55, min=0, max=100)    #Clearance rate of infected cells
params.add('gamma', value=0.24, min=0, max=10)        #Infected cells release virus at rate gamma
params.add('delta', value=5, min=0, max=10)     #clearance rate of virus particles
params.add('mu', value=6, min=0, max=10)
params.add('s', value=1, min=0.99, max=1.01)
params.add('epsilon', value=1, min=0.9, max=1.1)
params.add('tau', value=0.5, min=0.49, max=0.51)

# fit model
result = minimize(residual, params, args=(t_measured, V_measured), method='leastsq')  # leastsq nelder
# check results of the fit
data_fitted = g(t_measured, y0, result.params)

# plot fitted data
#plt.figure()
# ax1.plot(t_measured, 10**(-6)*data_fitted[:, 1], '-', linewidth=2, color='green', label='fitted I1 data')
# ax1.plot(t_measured, 10**(-6)*data_fitted[:, 2], '-', linewidth=2, color='blue', label='fitted I2 data')
ax1.plot(t_measured, 10**(-6)*data_fitted[:, 3], '-', linewidth=2, color='red', label='fitted V data')
#ax1.plot(t_measured, 10**(-6)*data_fitted[:, 4], '-', linewidth=2, color='brown', label='fitted F data')
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
log_I1_fitted = np.log10(data_fitted[:, 1])
log_I2_fitted = np.log10(data_fitted[:, 2])
log_V_fitted = np.log10(data_fitted[:, 3])
log_F_fitted = np.log10(data_fitted[:, 4])
#plt.figure()
ax2.scatter(t_measured[1:], log_V_measured[1:], marker='o', color='red', label='measured V data', s=75) #the first point is found by extrapolation. Therefore it is not physical so dont plot it.
ax2.plot(t_measured, log_I1_fitted, '-', linewidth=2, color='green', label='fitted I1 data')
ax2.plot(t_measured, log_I2_fitted, '-', linewidth=2, color='blue', label='fitted I2 data')
ax2.plot(t_measured, log_V_fitted, '-', linewidth=2, color='red', label='fitted V data')
ax2.plot(t_measured, log_F_fitted, '-', linewidth=2, color='brown', label='fitted F data')

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
I1_fitted = data_fitted[:, 1]
ax3.scatter(t_measured[1:], 10**(-6)*V_measured[1:], marker='o', color='red', label='measured V data', s=75) #the first point is found by extrapolation. Therefore it is not physical so dont plot it.
ax3.plot(t_measured, 10**(-6)*I1_fitted, '-', linewidth=2, color='green', label='fitted I1 data')
U_fitted = data_fitted[:, 0]
I2_fitted = data_fitted[:, 2]
V_fitted = data_fitted[:, 3]
F_fitted = data_fitted[:, 4]
ax3.plot(t_measured, 10**(-6)*U_fitted, '-', linewidth=2, color='black', label='fitted U data')
ax3.plot(t_measured, 10**(-6)*I2_fitted, '-', linewidth=2, color='blue', label='fitted I2 data')
ax3.plot(t_measured, 10**(-6)*V_fitted, '-', linewidth=2, color='red', label='fitted V data')
ax3.plot(t_measured, 10**(-6)*F_fitted, '-', linewidth=2, color='brown', label='fitted F data')
#plt.ylim(bottom=0.9 * min(log_V_measured))
ax3.set_xlim(left=0)
#ax3.set_ylim([0, 1.1 * 10**(-6)*max(Id_measured)])
ax3.legend()
ax3.set_xlabel('Days Post Infection')
ax3.set_ylabel('Concentration (million copies/mL)')
ax3.set_title('c)')

#try to recreate figure from the paper
# plt.figure()
# plt.scatter(t_measured, log_V_measured, marker='o', color='black', label='measured V data', s=75)
# plt.plot(t_measured, log_V_fitted, '-', linewidth=2, color='black', label='fitted V data')
# plt.plot(t_measured, log_F_fitted, '-', linewidth=2, color='red', label='fitted F data')
# plt.legend()

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
