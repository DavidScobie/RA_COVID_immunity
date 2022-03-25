import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn
from lmfit import minimize, Parameters, Parameter, report_fit
from scipy.integrate import odeint

##################################### process the TS data

#import the excel data
df = pd.read_excel ('C:\Research_Assistant\work\FW__Human_challenge_studies\COVHIC001_qPCR_TS.xlsx')

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

#get rid of the DETECTED AND NOT DETECTED
df2_str_sorted['length'] = df2_str_sorted['Virus Titre (Log10 copies/mL)'].str.len()
df2_str_sorted = df2_str_sorted[df2_str_sorted.length < 7]

#convert the strings to numbers
df2_str_sorted['Virus Titre (Log10 copies/mL)'] = pd.to_numeric(df2_str_sorted['Virus Titre (Log10 copies/mL)'], downcast="float")
df2_str_sorted['Study Day'] = pd.to_numeric(df2_str_sorted['Study Day'], downcast="float")

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

#convert the virus numbers to a list
vir_list_Non_DET = df2_str_sorted['Virus Titre (Log10 copies/mL)'].tolist()
#print('vir_list_Non_DET',len(vir_list_Non_DET),'vir_list_Non_DET',vir_list_Non_DET)

##plot the means

#find all the possible effective_day values
eff_day_vals_TS = list(set(effective_day))
eff_day_vals_TS.sort()  #THIS ONLY WORKS IF THERE IS AT LEAST 1 COUNT AT EACH TIME POINT
#print('eff_day_vals_TS',eff_day_vals_TS)

#find the occurences of each of the days
occ=np.zeros(len(eff_day_vals_TS))
for j in effective_day:
    for i in eff_day_vals_TS:
        if i==j:
            occ[int(2*(i-min(eff_day_vals_TS)))]+=1   #0.5 gap between vals, and begin at min val
#print('occ',occ)


#divide virus amount by number of counts on that day
div_vir_list=[]
k=0
for j in effective_day:
    for i in eff_day_vals_TS:
        if i==j:
            div_vir_list.append(vir_list_Non_DET[int(k)]/occ[int(2*(i-min(eff_day_vals_TS)))])
            k+=1
#print('div_vir_list',div_vir_list)


#sum the virus amounts on their specific day
div_vir_list_sum = np.zeros(len(eff_day_vals_TS))
k=0
for j in effective_day:
    for i in eff_day_vals_TS:
        if i==j:
            div_vir_list_sum[int(2*(i-min(eff_day_vals_TS)))]+=div_vir_list[int(k)]
            k+=1
#print('div_vir_list_sum',div_vir_list_sum)

##how many patients do we have? do the patients get sick or stay healthy or both? (out of the 36)

#find all the possible subject IDs
Subject_ID = df2_str_sorted['Subject ID'].tolist()
Subject_ID_vals = list(set(Subject_ID))
#print('Subject_ID',Subject_ID)

#append effective_day to the dataframe
df2_str_sorted['effective_day'] = effective_day

#plot the subjects in different colours
df2_str_sorted['Subject ID'] = df2_str_sorted['Subject ID'].astype(str)

#plot actual virus amount (instead of log10 of virus amount)
act_div_vir_list_sum_TS = np.zeros(len(div_vir_list_sum))
for i in range (len(div_vir_list_sum)):
    act_div_vir_list_sum_TS[i] = 10**(div_vir_list_sum[i])

plt.figure()
plt.plot(eff_day_vals_TS,act_div_vir_list_sum_TS,'-rx')
plt.xlabel('Days Post Infection')
plt.ylabel('Virus Titre (copies/mL)')

######################################## process the MTS data

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
df2_str_sorted = df2.reindex(s)

#get rid of the DETECTED AND NOT DETECTED
df2_str_sorted['length'] = df2_str_sorted['Virus Titre (Log10 copies/mL)'].str.len()
df2_str_sorted = df2_str_sorted[df2_str_sorted.length < 7]

#convert the strings to numbers
df2_str_sorted['Virus Titre (Log10 copies/mL)'] = pd.to_numeric(df2_str_sorted['Virus Titre (Log10 copies/mL)'], downcast="float")
df2_str_sorted['Study Day'] = pd.to_numeric(df2_str_sorted['Study Day'], downcast="float")

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

#convert the virus numbers to a list
vir_list_Non_DET = df2_str_sorted['Virus Titre (Log10 copies/mL)'].tolist()
#print('vir_list_Non_DET',len(vir_list_Non_DET),'vir_list_Non_DET',vir_list_Non_DET)

##plot the means

#find all the possible effective_day values
eff_day_vals_MTS = list(set(effective_day))
eff_day_vals_MTS.sort()  #THIS ONLY WORKS IF THERE IS AT LEAST 1 COUNT AT EACH TIME POINT
#print('eff_day_vals_MTS',eff_day_vals_MTS)

#find the occurences of each of the days
occ=np.zeros(len(eff_day_vals_MTS))
for j in effective_day:
    for i in eff_day_vals_MTS:
        if i==j:
            occ[int(2*(i-min(eff_day_vals_MTS)))]+=1   #0.5 gap between vals, and begin at min val
#print('occ',occ)


#divide virus amount by number of counts on that day
div_vir_list=[]
k=0
for j in effective_day:
    for i in eff_day_vals_MTS:
        if i==j:
            div_vir_list.append(vir_list_Non_DET[int(k)]/occ[int(2*(i-min(eff_day_vals_MTS)))])
            k+=1
#print('div_vir_list',div_vir_list)


#sum the virus amounts on their specific day
div_vir_list_sum = np.zeros(len(eff_day_vals_MTS))
k=0
for j in effective_day:
    for i in eff_day_vals_MTS:
        if i==j:
            div_vir_list_sum[int(2*(i-min(eff_day_vals_MTS)))]+=div_vir_list[int(k)]
            k+=1
#print('div_vir_list_sum',div_vir_list_sum)

##how many patients do we have? do the patients get sick or stay healthy or both? (out of the 36)

#find all the possible subject IDs
Subject_ID = df2_str_sorted['Subject ID'].tolist()
Subject_ID_vals = list(set(Subject_ID))
#print('Subject_ID',Subject_ID)

#append effective_day to the dataframe
df2_str_sorted['effective_day'] = effective_day

#plot the subjects in different colours
df2_str_sorted['Subject ID'] = df2_str_sorted['Subject ID'].astype(str)

#plot actual virus amount (instead of log10 of virus amount)
act_div_vir_list_sum_MTS = np.zeros(len(div_vir_list_sum))
for i in range (len(div_vir_list_sum)):
    act_div_vir_list_sum_MTS[i] = 10**(div_vir_list_sum[i])

plt.figure()
plt.plot(eff_day_vals_MTS,act_div_vir_list_sum_MTS,'-rx')
plt.xlabel('Days Post Infection')
plt.ylabel('Virus Titre (copies/mL)')


#######################################################
#model the data via differential equations

def f(y, t, paras):

    #Your system of differential equations

    #the U, V and I
    U = y[0]
    I_m = y[1]
    I_n = y[2]
    V_m = y[3]
    V_n = y[4]

    #the parameters alpha, beta, gamma, delta
    try:
        alpha_m = paras['alpha_m'].value
        alpha_n = paras['alpha_n'].value
        beta = paras['beta'].value
        gamma = paras['gamma'].value
        delta = paras['delta'].value

    except KeyError:
        alpha, beta, gamma, delta = paras

    # the model equations
    dUdt = - ( alpha_m * V_m * U ) - ( alpha_n * V_n * U)
    dImdt = (alpha_m * V_m * U) - (beta * I_m)
    dIndt = (alpha_n * V_n * U) - (beta * I_n)
    dVmdt = ( gamma * I_m ) - ( delta * V_m )
    dVndt = ( gamma * I_n ) - (delta * V_n )
    return [dUdt, dImdt, dIndt, dVmdt, dVndt]

def g(t, x0, paras):

    #Solution to the ODE x'(t) = f(t,x,k) with initial condition x(0) = x0

    x = odeint(f, x0, t, args=(paras,))
    return x


def residual(paras, t, data):


    #compute the residual between actual data and fitted data


    x0 = paras['U0'].value, paras['I_m0'].value, paras['I_n0'].value, paras['V_m0'].value, paras['V_n0'].value
    model = g(t, x0, paras)

    # you have data for two of your variables
    V_m_model = model[:, 3]
    V_n_model = model[:, 4]

    #want to find the residual between the log of the virus measured and fitted data for mouth and nose
    log_V_m_model = np.log10(V_m_model)
    log_data = np.log10(data)

    return ( (log_V_m_model - log_data) + () ).ravel()


# initial conditions

#my optimised initial conditions
U0 = 4*(10**(8))  #the number of cells in an adult is 4x10^8
I_m0 = 0   #Should be zero
I_n0 = 0   #Should be zero
V_m0 = act_div_vir_list_sum_TS[0]   #just taking the first measured value
V_n0 = act_div_vir_list_sum_MTS[0]   #just taking the first measured value
y0 = [U0, I_m0, I_n0, V_m0, V_n0]

# measured data
t_measured_TS = eff_day_vals_TS
t_measured_MTS = eff_day_vals_MTS
V_measured_TS = act_div_vir_list_sum_TS
V_measured_MTS = act_div_vir_list_sum_MTS

#plt.figure()
fig, (ax1, ax2, ax3) = plt.subplots(1,3)
ax1.scatter(t_measured_TS, 10**(-6)*V_measured_TS, marker='o', color='red', label='measured V_TS data', s=75)

# set parameters including bounds; you can also fix parameters (use vary=False)
params = Parameters()
params.add('U0', value=U0, vary=False)
params.add('I_m0', value=I_m0, vary=False)
params.add('I_n0', value=I_n0, vary=False)
params.add('V_m0', value=V_m0, vary=False)
params.add('V_n0', value=V_n0, vary=False)

# #parameters optimised on first 6 days of data
# params.add('alpha', value=4.24*(10**(-7)), min=4.23*(10**(-7)), max=4.25*(10**(-7)))   #rate that viral particles infect susceptible cells
# params.add('beta', value=61.2, min=61.1, max=61.3)    #Clearance rate of infected cells
# params.add('gamma', value=1.83, min=1.82, max=1.84)        #Infected cells release virus at rate gamma
# params.add('delta', value=1.45, min=1.44, max=1.46)     #clearance rate of virus particles

#my optimised parameters
params.add('alpha_m', value=6.63*(10**(-7)), min=1*(10**(-8)), max=9*(10**(-6)))   #rate that viral particles infect susceptible cells
params.add('alpha_n', value=6.63*(10**(-7)), min=1*(10**(-8)), max=9*(10**(-6)))   #rate that viral particles infect susceptible cells
params.add('beta', value=56, min=0, max=75)    #Clearance rate of infected cells
params.add('gamma', value=0.66, min=0, max=6)        #Infected cells release virus at rate gamma
params.add('delta', value=0.51, min=0, max=100)     #clearance rate of virus particles

# fit model
result = minimize(residual, params, args=(t_measured_TS, t_measured_MTS, V_measured_TS, V_measured_MTS), method='leastsq')  # leastsq nelder
# check results of the fit
data_fitted = g(t_measured, y0, result.params)

# plot fitted data
#plt.figure()
ax1.plot(t_measured, 10**(-6)*data_fitted[:, 1], '-', linewidth=2, color='red', label='fitted V data')
ax1.legend()
ax1.set_xlim([0, max(t_measured)])
ax1.set_ylim([0, 1.1 * 10**(-6)*max(V_measured)])
ax1.set_xlabel('Days Post Infection')
ax1.set_ylabel('Virus Titre Concentration (million copies/mL)')
ax1.set_title('a)')
# display fitted statistics
report_fit(result)

#collect the parameters from the overall model
overall_alpha=[]
overall_beta=[]
overall_gamma=[]
overall_delta=[]
for name, param in result.params.items():
    print(f'{name:7s} {param.value:11.5f} {param.stderr:11.5f}')
    if name == 'alpha':
        overall_alpha.append(param.value)
    if name == 'beta':
        overall_beta.append(param.value)
    if name == 'gamma':
        overall_gamma.append(param.value)
    if name == 'delta':
        overall_delta.append(param.value)

#plot the fitted data and the model for log(virus) against day
log_V_measured = np.log10(V_measured)
log_V_fitted = np.log10(data_fitted[:, 1])
V_fitted = data_fitted[:, 1]
#plt.figure()
ax2.scatter(t_measured, log_V_measured, marker='o', color='red', label='measured qPCR V data', s=75)
ax2.plot(t_measured, log_V_fitted, '-', linewidth=2, color='red', label='fitted qPCR V data')
ax2.set_xlim(left=0)
ax2.set_xlabel('Days Post Infection')
ax2.set_ylabel('Virus Titre Concentration (Log10 copies/mL)')
ax2.set_title('b)')
print('log_V_measured',log_V_measured, 'LENGTH log_V_measured',len(log_V_measured),'IT IS SAVED HERE')
#np.save('TS_log_V_measured', log_V_measured)
#np.save('TS_t_measured', t_measured)
print('t_measured',t_measured)

#####plot the FFA data on top
# FFA_virus=np.array([78.22279965, 586.58819126, 641.94820653, 2290.86793783, 239.88332219, 194.98447853, 429.86622684, 245.47088348, 218.77614918, 169.04409205, 50.11872428])
# log_FFA_virus = np.log10(FFA_virus)
# FFA_effective_day = np.array([3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0])
# ax2.plot(FFA_effective_day, log_FFA_virus, marker='o', color='black', label='measured FFA V data')
# ax2.legend()
##

#plot the measured data, along with the fitted model for V, I and U
#plt.figure()
ax3.scatter(t_measured, 10**(-6)*V_measured, marker='o', color='red', label='measured V data', s=75)
ax3.plot(t_measured, 10**(-6)*V_fitted, '-', linewidth=2, color='red', label='fitted V data')
U_fitted = data_fitted[:, 0]
I_fitted = data_fitted[:, 2]
ax3.plot(t_measured, 10**(-6)*U_fitted, '-', linewidth=2, color='green', label='fitted U data')
ax3.plot(t_measured, 10**(-6)*I_fitted, '-', linewidth=2, color='blue', label='fitted I data')
#plt.ylim(bottom=0.9 * min(log_V_measured))
ax3.set_xlim(left=0)
ax3.set_ylim([0, 1.1 * 10**(-6)*max(V_measured)])
ax3.legend()
ax3.set_xlabel('Days Post Infection')
ax3.set_ylabel('Concentration (million copies/mL)')
ax3.set_title('c)')

# ax3.scatter(t_measured, log_V_measured, marker='o', color='red', label='measured V data', s=75)
# ax3.plot(t_measured, log_V_fitted, '-', linewidth=2, color='red', label='fitted V data')
# log_U_fitted = np.log10(data_fitted[:, 0])
# log_I_fitted = np.log10(data_fitted[:, 2])
# ax3.plot(t_measured, log_U_fitted, '-', linewidth=2, color='green', label='fitted U data')
# ax3.plot(t_measured, log_I_fitted, '-', linewidth=2, color='blue', label='fitted I data')
# #plt.ylim(bottom=0.9 * min(log_V_measured))
# ax3.set_xlim(left=0)
# ax3.legend()
# ax3.set_xlabel('Days Post Infection')
# ax3.set_ylabel('Concentration (Log10 copies/mL)')
# ax3.set_title('c)')


plt.show()

