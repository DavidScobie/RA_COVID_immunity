import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn
from lmfit import minimize, Parameters, Parameter, report_fit
from scipy.integrate import odeint
from scipy.stats import norm, lognorm
import matplotlib.mlab as mlab
import scipy.stats as stats
import statistics
from matplotlib import cm

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

    k+=1
    df2_Subj_ID_sorted = df2_str_sorted.loc[df2_str_sorted['Subject ID'] == j] #make a subset of the dataframe based on the subject ID

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

    #if number of NON DETs less than 15, then cut this out from the bunch
    if number_NON_DETs < 15: #the case where I want to keep the data
        m+=1
        if m == 1: #the first time through this loop we just create the dataframe
            print('THE BEGINNING')
            df2_cut_out_many = df2_Subj_ID_sorted_comb
        else: #for all of the subsequent patients with few NON DETs, we append their data to the dataframe
            print('APPENDING')
            df2_cut_out_many = df2_cut_out_many.append(df2_Subj_ID_sorted_comb)
    else:
        print('MORE THAN 15 NON DETs')

#convert the strings to numbers
df2_str_sorted = df2_cut_out_many
df2_str_sorted['Virus Titre (Log10 copies/mL)'] = pd.to_numeric(df2_str_sorted['Virus Titre (Log10 copies/mL)'], downcast="float")
df2_str_sorted['Study Day'] = pd.to_numeric(df2_str_sorted['Study Day'], downcast="float")

################### create the 'effective study day' which takes AM and PM into account

day_list = df2_str_sorted['Study Day'].tolist() #convert study day and Timepoint into lists
Timepoint_list = df2_str_sorted['Timepoint'].tolist()

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

###########plot the means

#find all the possible effective_day values
eff_day_vals = list(set(effective_day))
eff_day_vals.sort()  #THIS ONLY WORKS IF THERE IS AT LEAST 1 COUNT AT EACH TIME POINT

#find the occurences of each of the days
occ=np.zeros(len(eff_day_vals))
for j in effective_day:
    for i in eff_day_vals:
        if i==j:
            occ[int(2*(i-min(eff_day_vals)))]+=1   #0.5 gap between vals, and begin at min val

#divide virus amount by number of counts on that day
div_vir_list=[]
k=0
for j in effective_day:
    for i in eff_day_vals:
        if i==j:
            div_vir_list.append(vir_list_Non_DET[int(k)]/occ[int(2*(i-min(eff_day_vals)))])
            k+=1

#sum the virus amounts on their specific day
div_vir_list_sum = np.zeros(len(eff_day_vals))
k=0
for j in effective_day:
    for i in eff_day_vals:
        if i==j:
            div_vir_list_sum[int(2*(i-min(eff_day_vals)))]+=div_vir_list[int(k)]
            k+=1

#############plot the individual subjects

#find all the possible subject IDs
Subject_ID = df2_str_sorted['Subject ID'].tolist()
Subject_ID_vals = list(set(Subject_ID))

#append effective_day to the dataframe
df2_str_sorted['effective_day'] = effective_day

#plot the subjects in different colours
df2_str_sorted['Subject ID'] = df2_str_sorted['Subject ID'].astype(str)

seaborn.relplot(data=df2_str_sorted, x='effective_day', y='Virus Titre (Log10 copies/mL)', hue='Subject ID')

# plt.figure()
# seaborn.pointplot(data=df2_str_sorted, x='effective_day', y='Virus Titre (Log10 copies/mL)', hue='Subject ID', ci=None)
"""
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

###########plot actual virus amount (instead of log10 of virus amount)
act_div_vir_list_sum = np.zeros(len(div_vir_list_sum))
for i in range (len(div_vir_list_sum)):
    act_div_vir_list_sum[i] = 10**(div_vir_list_sum[i])

# plt.figure()
# plt.plot(eff_day_vals,act_div_vir_list_sum,'-rx')
# plt.xlabel('Days Post Infection')
# plt.ylabel('Virus Titre (copies/mL)')

#######################################################

def f(y, t, paras):
    """
    Your system of differential equations
    """
    #the U, Id and Is
    U = y[0]
    Id = y[1]
    Is = y[2]

    #the parameters alpha, beta, gamma, delta
    try:
        alpha = paras['alpha'].value
        beta = paras['beta'].value
        gamma = paras['gamma'].value
        delta = paras['delta'].value

    except KeyError:
        alpha, beta, gamma, delta = paras

    # the model equations
    f0 = - alpha * U * Is   #dU/dt
    f1 = (gamma * Is) - (delta * Id)    #d(Id)/dt
    f2 = (alpha * U * Is) - (beta * Is) - (gamma * Is)    #d(Is)/dt
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

    x0 = paras['U0'].value, paras['Id0'].value, paras['Is0'].value
    model = g(t, x0, paras)

    # we now have data for the sum of Id and Is
    Id_Is_model = model[:, 1] + model[:, 2]

    #want to find the residual between the log of the virus measured and fitted data
    log_Id_IS_model = np.log10(Id_Is_model)
    log_data = np.log10(data)

    #so we are now minimising the residual between the data and the sum of Is and Id
    return (log_Id_IS_model - log_data).ravel()


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
max_indic_arr = np.where(act_div_vir_list_sum == np.amax(act_div_vir_list_sum))
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

#add the point at time=0, virus=initial to the eff_day_vals and act_div_vir_list_sum arrays
v1 = 0
v2 = 10**a #THIS IS 10**2.97. 2.97 is the y intercept of the line of best fit
eff_day_vals = np.insert(eff_day_vals, 0, v1, axis=0)
act_div_vir_list_sum = np.insert(act_div_vir_list_sum, 0, v2, axis=0)
print('eff_day_vals',eff_day_vals,'act_div_vir_list_sum',act_div_vir_list_sum)

#my optimised initial conditions
U0 = 4*(10**(8))  #the number of cells in an adult is 4x10^8
Is0 = act_div_vir_list_sum[0]
Is0_init = act_div_vir_list_sum[0] #for later on when I need to refer back to initial
Id0 = 0
y0 = [U0, Id0, Is0]
y0_init = [U0, Id0, Is0] #for later on when I need to refer back to initial

# measured data
t_measured = eff_day_vals
t_measured_init = t_measured #this is the timepoints for the average of the data across all patients
V_measured = act_div_vir_list_sum
V_measured_init = V_measured #this is the virus for the average of the data across all patients

#########################################################################   ABOVE THIS IS FINE

#plt.figure()
fig, (ax1, ax2, ax3) = plt.subplots(1,3)
ax1.scatter(t_measured, 10**(-6)*V_measured, marker='o', color='red', label='measured V data', s=75)

# set parameters including bounds; you can also fix parameters (use vary=False)
params = Parameters()
params.add('U0', value=U0, vary=False)
params.add('Id0', value=Id0, vary=False)
params.add('Is0', value=Is0, vary=False)

#my optimised parameters
params.add('alpha', value=5.5*(10**(-7)), min=7.99*(10**(-8)), max=8.01*(10**(-7)))   #rate that viral particles infect susceptible cells
params.add('beta', value=1*(10**(-11)), min=0, max=1.1*(10**(-11)))    #Clearance rate of infected cells
params.add('gamma', value=214, min=0, max=400)        #Infected cells release virus at rate gamma
params.add('delta', value=0.86, min=0, max=100)     #clearance rate of virus particles

# fit model
result = minimize(residual, params, args=(t_measured, V_measured), method='leastsq')  # leastsq nelder
# check results of the fit
data_fitted = g(t_measured, y0, result.params)

##############################

# plot fitted data
ax1.plot(t_measured, 10**(-6)*data_fitted[:, 1], '-', linewidth=2, color='green', label='fitted Id data')
ax1.plot(t_measured, 10**(-6)*data_fitted[:, 2], '-', linewidth=2, color='blue', label='fitted Is data')
ax1.plot(t_measured, (10**(-6)*data_fitted[:, 2]) + (10**(-6)*data_fitted[:, 1]), '-', linewidth=2, color='red', label='fitted (Is + Id) data')
ax1.legend()
ax1.set_xlim([0, max(t_measured)])
#ax1.set_ylim([0, 1.1 * 10**(-6)*max(V_measured)])
ax1.set_xlabel('Days Post Infection')
ax1.set_ylabel('Virus Titre Concentration (million copies/mL)')
ax1.set_title('a)')
# display fitted statistics
report_fit(result)

#compute the variance
print('chisqr',result.chisqr)
print('ndata',result.ndata)
overall_variance = (result.chisqr) / (result.ndata) #(chi_squ / N)
print('overall_variance',overall_variance)

#plot the fitted data and the model for log(virus) against day
log_V_measured = np.log10(V_measured)
log_Id_fitted = np.log10(data_fitted[:, 1])
log_Is_fitted = np.log10(data_fitted[:, 2])
log_Id_Is_fitted = np.log10(data_fitted[:, 1] + data_fitted[:, 2])
#plt.figure()
ax2.scatter(t_measured[1:], log_V_measured[1:], marker='o', color='red', label='measured (Id+Is) data', s=75) #the first point is found by extrapolation. Therefore it is not physical so dont plot it.
ax2.plot(t_measured, log_Id_fitted, '-', linewidth=2, color='green', label='fitted Id data')
ax2.plot(t_measured, log_Is_fitted, '-', linewidth=2, color='blue', label='fitted Is data')
ax2.plot(t_measured, log_Id_Is_fitted, '-', linewidth=2, color='red', label='fitted (Id + Is) data')

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
Id_fitted = data_fitted[:, 1]
ax3.scatter(t_measured[1:], 10**(-6)*V_measured[1:], marker='o', color='red', label='measured (Id + Is) data', s=75) #the first point is found by extrapolation. Therefore it is not physical so dont plot it.
ax3.plot(t_measured, 10**(-6)*Id_fitted, '-', linewidth=2, color='green', label='fitted Id data')
U_fitted = data_fitted[:, 0]
Is_fitted = data_fitted[:, 2]
ax3.plot(t_measured, 10**(-6)*U_fitted, '-', linewidth=2, color='black', label='fitted U data')
ax3.plot(t_measured, 10**(-6)*Is_fitted, '-', linewidth=2, color='blue', label='fitted Is data')
#plt.ylim(bottom=0.9 * min(log_V_measured))
ax3.set_xlim(left=0)
#ax3.set_ylim([0, 1.1 * 10**(-6)*max(Id_measured)])
ax3.legend()
ax3.set_xlabel('Days Post Infection')
ax3.set_ylabel('Concentration (million copies/mL)')
ax3.set_title('c)')

############################

I_area = np.trapz(data_fitted[:, 1], dx=0.5)
print('I_area',I_area)

##########################################
#fit models to different patients

Subject_ID_vals_short = Subject_ID_vals[0:-3]  #only want to model the first 15 patients as these are the training dat set
print('Subject_ID_vals_short',Subject_ID_vals_short)

#initialise arrays of patient parameters
alphas=[]
betas=[]
gammas = []
deltas = []
red_chi_squs = []
residuals = []
sum_residuals_squs = []
chi_squs = []
ndatas = []
variances = []
subj_IDs_over_5=[]

for j in Subject_ID_vals_short:

    df2_Subj_ID_sorted = df2_str_sorted[df2_str_sorted['Subject ID'].str.contains(str(j)) == True]  #make a subset of the dataframe based on the subject ID
    df2_Subj_ID_sub_eff_sort = df2_Subj_ID_sorted.sort_values(["effective_day"], ascending=True) #sort the values of the dataframe based on the effective_day

    if len(df2_Subj_ID_sub_eff_sort['Virus Titre (Log10 copies/mL)'].tolist()) > 5: #only use the subjects with more than 5 data points
        k+=1
        #convert the virus and the effective day values to a list
        div_vir_list_sum = df2_Subj_ID_sub_eff_sort['Virus Titre (Log10 copies/mL)'].tolist()
        eff_day_list = df2_Subj_ID_sub_eff_sort['effective_day'].tolist()

        print('SUBJECT ID ',j)

        #compute the actual virus amount (not the log)
        act_div_vir_list_sum = np.zeros(len(div_vir_list_sum))
        for i in range (len(div_vir_list_sum)):
            act_div_vir_list_sum[i] = 10**(div_vir_list_sum[i])

        #we want the first data point to be the extrapolation on the average of everyone, not the extrapolation on the patient in question

        eff_day_list = np.insert(eff_day_list, 0, v1, axis=0)
        act_div_vir_list_sum = np.insert(act_div_vir_list_sum, 0, v2, axis=0)

        #my optimised initial conditions
        U0 = 4*(10**(8))  #the number of target cells in an adult is 4x10^8
        Is0 = act_div_vir_list_sum[0]
        Id0 = 0
        y0 = [U0, Id0, Is0]

        # measured data
        t_measured = eff_day_list
        V_measured = act_div_vir_list_sum

        # plt.figure()
        # plt.scatter(t_measured, V_measured, marker='o', color='red', label='measured V data', s=75)

        # set parameters including bounds; you can also fix parameters (use vary=False)
        params = Parameters()
        params.add('U0', value=U0, vary=False)
        params.add('Id0', value=Id0, vary=False)
        params.add('Is0', value=Is0, vary=False)

        #my optimised parameters
        params.add('alpha', value=5.5*(10**(-7)), min=2*(10**(-8)), max=8.01*(10**(-6)))   #rate that viral particles infect susceptible cells
        params.add('beta', value=1*(10**(-11)), min=0, max=1.1*(10**(-11)))    #Clearance rate of infected cells
        params.add('gamma', value=214, min=0, max=600)        #Infected cells release virus at rate gamma
        params.add('delta', value=0.86, min=0, max=100)     #clearance rate of virus particles


        # fit model
        result = minimize(residual, params, args=(t_measured, V_measured), method='leastsq')  # leastsq nelder
        # check results of the fit
        data_fitted = g(t_measured, y0, result.params)

        # display fitted statistics and append parameters to lists
        subj_IDs_over_5.append(j)
        report_fit(result)
        for name, param in result.params.items():
            if name == 'alpha':
                alphas.append(param.value)
            if name == 'beta':
                betas.append(param.value)
            if name == 'gamma':
                gammas.append(param.value)
            if name == 'delta':
                deltas.append(param.value)

        red_chi_squs.append(result.redchi)
        residuals.append(result.residual)
        sum_residuals_squs.append(sum(((result.residual)**2)))
        chi_squs.append(result.chisqr)
        ndatas.append(result.ndata)

        #plot the fitted data and the model for log(virus) against day
        log_V_measured = np.log10(V_measured)
        log_Id_fitted = np.log10(data_fitted[:, 1])
        log_Is_fitted = np.log10(data_fitted[:, 2])
        log_Id_Is_fitted = np.log10(data_fitted[:, 1] + data_fitted[:, 2])

        #########plot the log of virus amount against time
        plt.figure()
        plt.scatter(t_measured[1:], log_V_measured[1:], marker='o', color='red', label='measured V data', s=75) #the first point is found by extrapolation. Therefore it is not physical so dont plot it.
        plt.plot(t_measured, log_Id_Is_fitted, '-', linewidth=2, color='red', label='fitted I data')
        plt.ylim(bottom=0.9 * min(log_V_measured))
        plt.xlim(left=0)
        plt.legend()
        plt.xlabel('Days Post Infection')
        plt.ylabel('Concentration (Log10 copies/mL)')
        plt.title('Subject ID=%i' %j)


print('alphas',alphas)
print('betas',betas)
print('gammas',gammas)
print('deltas',deltas)
#print('residuals',residuals)
print('sum_residuals_squs',sum_residuals_squs)
print('ndatas',ndatas)
variances = np.array(sum_residuals_squs) / np.array(ndatas)
print('variances',variances)
print('average variance',sum(variances)/len(variances))


#only include patients who have variance less than a value
refined_alphas = []
refined_betas = []
refined_gammas = []
refined_deltas = []

for i in range (len(variances)):
    if variances[i]<=100: #the value here is the cut off for the variance
        refined_alphas.append(alphas[i])
        refined_betas.append(betas[i])
        refined_gammas.append(gammas[i])
        refined_deltas.append(deltas[i])

##############################################################################

####Find the median of the alpha values
alpha_med = np.median(alphas)

########use a normal distribution to compute the random effect and find the new adjusted alpha values (hopefully in a lognormal dist)
adj_alphas = []
for i in range (len(alphas)): #length 18 for all the patients
    exponent = np.random.normal(loc=0.0, scale=0.4) #this is the random effect. Randomly sample from a normal distribution with mean zero and w=0.4
    adj_alph = alpha_med*np.exp(exponent) #add the fixed term onto the random term
    adj_alphas.append(adj_alph) #append it to an array

print('adj_alphas',adj_alphas)
plt.figure()
plt.hist(adj_alphas, density=False, bins=8)
plt.xlabel('Alpha value')
plt.ylabel('Density of alpha values')

####Find the median of the beta values
beta_med = np.median(betas)

########use a normal distribution to compute the random effect and find the new adjusted beta values (hopefully in a lognormal dist)
adj_betas = []
for i in range (len(betas)): #length 18 for all the patients
    exponent = np.random.normal(loc=0.0, scale=0.4) #this is the random effect. Randomly sample from a normal distribution with mean zero and w=0.4
    adj_bet = beta_med*np.exp(exponent) #add the fixed term onto the random term
    adj_betas.append(adj_bet) #append it to an array

print('adj_betas',adj_betas)
plt.figure()
plt.hist(adj_betas, density=False, bins=8)
plt.xlabel('Beta value')
plt.ylabel('Density of beta values')

####Find the median of the gamma values
gamma_med = np.median(gammas)

########use a normal distribution to compute the random effect and find the new adjusted beta values (hopefully in a lognormal dist)
adj_gammas = []
for i in range (len(gammas)): #length 18 for all the patients
    exponent = np.random.normal(loc=0.0, scale=0.4) #this is the random effect. Randomly sample from a normal distribution with mean zero and w=0.4
    adj_gam = gamma_med*np.exp(exponent) #add the fixed term onto the random term
    adj_gammas.append(adj_gam) #append it to an array

print('adj_gammas',adj_gammas)
plt.figure()
plt.hist(adj_gammas, density=False, bins=8)
plt.xlabel('gamma value')
plt.ylabel('Density of gamma values')

####Find the median of the delta values
delta_med = np.median(deltas)

########use a normal distribution to compute the random effect and find the new adjusted beta values (hopefully in a lognormal dist)
adj_deltas = []
for i in range (len(deltas)): #length 18 for all the patients
    exponent = np.random.normal(loc=0.0, scale=0.4) #this is the random effect. Randomly sample from a normal distribution with mean zero and w=0.4
    adj_del = delta_med*np.exp(exponent) #add the fixed term onto the random term
    adj_deltas.append(adj_del) #append it to an array

print('adj_deltas',adj_deltas)
plt.figure()
plt.hist(adj_deltas, density=False, bins=8)
plt.xlabel('delta value')
plt.ylabel('Density of delta values')

###################### Need to rethink this. It will be different for 3 parameters. The space may need to become 3d for alpha, beta and gamma, and have colour coding for BIC?
########alternatively dont make a plot, just store the data in a bigger matrix. Will need an extra nested loop for the gamma. Will need some thinking for this
###############also have a thik about whether you really need that -1/2 in the likelihood equation. Need to scale this properly with the kln(n) part with BIC
"""
sn=1 #the error on each point
k_param=4 # the number of parameters in the model
n_points = len(t_measured_init) #the number of data points for the average of all patients

how_many_points_gam = 2
how_many_points_alph_bet = 5
how_many_points_del = 2

range_alph = np.max(adj_alphas) - np.min(adj_alphas)
range_bet = np.max(adj_betas) - np.min(adj_betas)
range_gam = np.max(adj_gammas) - np.min(adj_gammas)
range_del = np.max(adj_deltas) - np.min(adj_deltas)

proportion = 1 #the proportion of parameter space that we want to explore (use this for unstable models)
shift_alpha = 0 #if we want to shift the parameters away from the centre of the alpha space then change this
shift_beta = 0 #if we want to shift the parameters away from the centre of the beta space then change this

alphas_to_surf = np.linspace(np.min(adj_alphas) + (range_alph*((1-proportion)/(2))) + (range_alph*shift_alpha), np.max(adj_alphas) - (range_alph*((1-proportion)/(2))) + (range_alph*shift_alpha), num=how_many_points_alph_bet)
betas_to_surf = np.linspace(np.min(adj_betas) + (range_bet*((1-proportion)/(2))) + (range_bet*shift_beta), np.max(adj_betas) - (range_bet*((1-proportion)/(2))) + (range_bet*shift_beta), num=how_many_points_alph_bet)
gammas_to_surf = np.linspace(np.min(adj_gammas) + (range_gam*((1-proportion)/(2))), np.max(adj_gammas) - (range_gam*((1-proportion)/(2))), num=how_many_points_gam)
deltas_to_surf = np.linspace(np.min(adj_deltas) + (range_del*((1-proportion)/(2))), np.max(adj_deltas) - (range_del*((1-proportion)/(2))), num=how_many_points_del)

print('alphas_to_surf',alphas_to_surf)
print('betas_to_surf',betas_to_surf)

######old code
# alphas_to_surf = np.linspace(np.min(adj_alphas), np.max(adj_alphas), num=how_many_points)
# betas_to_surf = np.linspace(np.min(adj_betas), np.max(adj_betas), num=how_many_points)
# gammas_to_surf = np.linspace(np.min(adj_gammas), np.max(adj_gammas), num=how_many_points)

BICs_all = [] #array for storing all of the BICs. In order to find the lowest one

print('gammas_to_surf',gammas_to_surf)
print('deltas_to_surf',deltas_to_surf)

for n in (deltas_to_surf):
    for m in (gammas_to_surf):
    #for m in (np.array([gamma_med])):
        #print('gamma = ',m)
        X, Y = np.meshgrid(alphas_to_surf, betas_to_surf)
        # print('X',X)
        BIC = [] #initialise the array of BIC
        sum_diff_squ = []
        for i in (alphas_to_surf):
            #print('alpha = ',i)
            for j in (betas_to_surf):
                #print('beta = ',j)

                params = Parameters()
                params.add('U0', value=U0, vary=False)
                params.add('V0', value=V0_init, vary=False)
                params.add('I0', value=I0, vary=False)

                #my optimised parameters
                params.add('alpha', value=i, min=i - 10**(-11), max=i + 10**(-11))   #rate that viral particles infect susceptible cells
                params.add('beta', value=j, min=j - 10**(-11), max=j + 10**(-11))
                params.add('gamma', value=m, min=m - 10**(-11), max=m + 10**(-11))
                params.add('delta', value=n, min=n - 10**(-11), max=n + 10**(-11))     #clearance rate of virus particles


                result = minimize(residual, params, args=(t_measured_init, V_measured_init), method='leastsq')  # leastsq nelder
                # check results of the fit
                data_fitted = g(t_measured_init, y0_init, result.params)

                #print('result.chisqr',result.chisqr)

                #plot the fitted data and the model for log(virus) against day
                log_V_measured = np.log10(V_measured_init)
                gn_Ftrue_log_I_fitted = np.log10(data_fitted[:, 1])

                #########plot the log of virus amount against time
                # plt.figure()
                # plt.scatter(t_measured_init[1:], log_V_measured[1:], marker='o', color='red', label='measured V data', s=75) #the first point is found by extrapolation. Therefore it is not physical so dont plot it.
                # plt.plot(t_measured_init, gn_Ftrue_log_I_fitted, '-', linewidth=2, color='red', label='fitted I data')
                # plt.ylim(bottom=0.9 * min(log_V_measured), top=9)
                # plt.xlim(left=0)
                # plt.legend()
                # plt.xlabel('Days Post Infection')
                # plt.ylabel('Concentration (Log10 copies/mL)')
                # plt.title("alpha={i}, beta={j}, gamma={m}, delta={n}".format(i=i, j=j, m=m, n=n))

                #print('t_measured_init',len(t_measured_init),'log_V_measured',len(log_V_measured),'gn_Ftrue_log_I_fitted',len(gn_Ftrue_log_I_fitted))

                #find differences between the gnFtrue virus amount and the virus amount for all of the models
                diff = [] #the difference between gn_Ftrue and Dn
                for k in range (len(t_measured_init)):
                    diff.append(log_V_measured[k] - gn_Ftrue_log_I_fitted[k])
                #print('np.square(diff)',np.square(diff))

                log_lik_term = -0.5*np.sum((np.square(diff)/(sn**2)) + np.log(2*np.pi*(sn**2)))
                BIC.append((k_param*np.log(n_points))-(2*log_lik_term))

        print('BIC',BIC)
        BIC_mat = np.reshape(BIC, (len(alphas_to_surf), len(betas_to_surf)))
        print('BIC_mat',BIC_mat)

        BICs_all.append(BIC)

        # Plot the surface.
        #plt.figure()
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        surf = ax.plot_surface(X, Y, BIC_mat, cmap=cm.coolwarm,
                            linewidth=0, antialiased=False)
        ax.set_xlabel('alpha')
        ax.set_ylabel('beta')
        ax.set_zlabel('BIC')
        #ax.set_title('gamma=%.6e' %m) # f represents a float
        ax.set_title("gamma={m}, delta={n}".format(m=m, n=n)) # f represents a float

flat_BICs_all = [food for sublist in BICs_all for food in sublist] #flatten the BIC list

min_indic = np.argmin(flat_BICs_all) #finding the lowest BIC
print('flat_BICs_all',flat_BICs_all,'min_indic',min_indic,'np.min(flat_BICs_all)',np.min(flat_BICs_all))
"""
plt.show()