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
    #the U, I
    U = y[0]
    I = y[1]

    #the parameters alpha, beta
    try:
        alpha = paras['alpha'].value
        beta = paras['beta'].value

    except KeyError:
        alpha, beta = paras

    # the model equations
    f0 = - alpha * U * I   #dU/dt
    f1 = (alpha * U * I) - (beta * I)    #dI/dt
    return [f0, f1]

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

    x0 = paras['U0'].value, paras['I0'].value
    model = g(t, x0, paras)

    # We have data for I (this is the virus data)
    I_model = model[:, 1]

    #want to find the residual between the log of the virus measured and fitted data
    log_I_model = np.log10(I_model)
    log_data = np.log10(data)

    #so we are now minimising the residual between the data and I
    return (log_I_model - log_data).ravel()


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
I0 = act_div_vir_list_sum[0]
I0_init = act_div_vir_list_sum[0] #for later on when I need to refer back to initial
y0 = [U0, I0]
y0_init = [U0, I0] #for later on when I need to refer back to initial

# measured data
t_measured = eff_day_vals
t_measured_init = t_measured #this is the timepoints for the average of the data across all patients
V_measured = act_div_vir_list_sum
V_measured_init = V_measured #this is the virus for the average of the data across all patients

#plt.figure()
fig, (ax1, ax2, ax3) = plt.subplots(1,3)
ax1.scatter(t_measured[1:], 10**(-6)*V_measured[1:], marker='o', color='red', label='measured V data', s=75) #the first point is found by extrapolation. Therefore it is not physical so dont plot it.

# set parameters including bounds; you can also fix parameters (use vary=False)
params = Parameters()
params.add('U0', value=U0, vary=False)
params.add('I0', value=I0, vary=False)

#my optimised parameters
params.add('alpha', value=1.9*(10**(-8)), min=1*(10**(-9)), max=6.3*(10**(-7)))   #rate that viral particles infect susceptible cells
params.add('beta', value=1.2*(10**(0)), min=0, max=1.1*(10**(2)))

# fit model
result = minimize(residual, params, args=(t_measured, V_measured), method='leastsq')  # leastsq nelder
# check results of the fit
data_fitted = g(t_measured, y0, result.params)

# plot fitted data
ax1.plot(t_measured, 10**(-6)*data_fitted[:, 1], '-', linewidth=2, color='red', label='fitted I data')
ax1.legend()
ax1.set_xlim([0, max(t_measured)])
ax1.set_xlabel('Days Post Infection')
ax1.set_ylabel('Virus Titre Concentration (million copies/mL)')
ax1.set_title('a)')
# display fitted statistics
report_fit(result)

#collect the parameters from the overall model
overall_alpha=[]
overall_beta=[]
for name, param in result.params.items():
    if name == 'alpha':
        overall_alpha.append(param.value)
    if name == 'beta':
        overall_beta.append(param.value)

#compute the variance
overall_variance = (result.chisqr) / (result.ndata) #(chi_squ / N)
print('overall_variance',overall_variance)

#plot the fitted data and the model for log(virus) against day
log_V_measured = np.log10(V_measured)
log_V_measured_init = np.log10(V_measured)
log_I_fitted = np.log10(data_fitted[:, 1])

#plt.figure()
ax2.scatter(t_measured[1:], log_V_measured[1:], marker='o', color='red', label='measured V data', s=75) #the first point is found by extrapolation. Therefore it is not physical so dont plot it.
ax2.plot(t_measured, log_I_fitted, '-', linewidth=2, color='red', label='fitted I data')

#################plotting curve from day -3
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

##################plot the measured data, along with the fitted model for V, I and U
I_fitted = data_fitted[:, 1]
ax3.scatter(t_measured[1:], 10**(-6)*V_measured[1:], marker='o', color='red', label='measured V data', s=75) #the first point is found by extrapolation. Therefore it is not physical so dont plot it.
ax3.plot(t_measured, 10**(-6)*I_fitted, '-', linewidth=2, color='red', label='fitted I data')
U_fitted = data_fitted[:, 0]
ax3.plot(t_measured, 10**(-6)*U_fitted, '-', linewidth=2, color='black', label='fitted U data')
ax3.set_xlim(left=0)
ax3.legend()
ax3.set_xlabel('Days Post Infection')
ax3.set_ylabel('Concentration (million copies/mL)')
ax3.set_title('c)')

############## find area under the Is and Id curves
# plt.figure()
# plt.plot(t_measured, data_fitted[:, 1], '-', linewidth=2, color='green', label='fitted Id data')
# # print('Id points',data_fitted[:, 1])
# plt.plot(t_measured, data_fitted[:, 2], '-', linewidth=2, color='blue', label='fitted Is data')
# # print('Is points',data_fitted[:, 2])
# plt.legend()
# plt.xlabel('Days Post Infection')
# plt.ylabel('Virus Titre (copies/mL)')

I_area = np.trapz(data_fitted[:, 1], dx=0.5)
print('I_area',I_area,'I_area',I_area)

#save arrays
#np.save('qPCR_MTS_V_measured_NON_DET_eq_zero_fit_Id+Is', V_measured)
#np.save('qPCR_MTS_t_measured_NON_DET_eq_zero_fit_Id+Is', t_measured)

##########################################
#fit models to different patients

Subject_ID_vals_short = Subject_ID_vals[0:-3]  #only want to model the first 15 patients as these are the training dat set
print('Subject_ID_vals_short',Subject_ID_vals_short)

#initialise arrays of patient parameters
alphas=[]
betas=[]
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

        ###extrapolation to find first data point
        max_indic_arr = np.where(act_div_vir_list_sum == np.amax(act_div_vir_list_sum)) # Get the indices of maximum element in eff_day_vals
        max_indic = int(max_indic_arr[0])

        a, b = best_fit(eff_day_list[:max_indic+1],np.log10(act_div_vir_list_sum)[:max_indic+1])

        #add the point at time=0, virus=10**a to the eff_day_vals and act_div_vir_list_sum arrays
        v1 = 0
        v2 = 10**a #THIS IS 10**a. where a is the y intercept of the line of best fit
        eff_day_list = np.insert(eff_day_list, 0, v1, axis=0)
        act_div_vir_list_sum = np.insert(act_div_vir_list_sum, 0, v2, axis=0)

        #my optimised initial conditions
        U0 = 4*(10**(8))  #the number of target cells in an adult is 4x10^8
        I0 = act_div_vir_list_sum[0]
        y0 = [U0, I0]

        # measured data
        t_measured = eff_day_list
        V_measured = act_div_vir_list_sum

        # plt.figure()
        # plt.scatter(t_measured, V_measured, marker='o', color='red', label='measured V data', s=75)

        # set parameters including bounds; you can also fix parameters (use vary=False)
        params = Parameters()
        params.add('U0', value=U0, vary=False)
        params.add('I0', value=I0, vary=False)

        #my optimised parameters
        params.add('alpha', value=1.9*(10**(-8)), min=1*(10**(-9)), max=6.3*(10**(-7)))   #rate that viral particles infect susceptible cells
        params.add('beta', value=1.2*(10**(0)), min=0, max=1.1*(10**(2)))

        # fit model
        result = minimize(residual, params, args=(t_measured, V_measured), method='leastsq')  # leastsq nelder
        # check results of the fit
        data_fitted = g(t_measured, y0, result.params)

        #####plot absolute values
        # plt.plot(t_measured, data_fitted[:, 1], '-', linewidth=2, color='red', label='fitted V data')
        # plt.legend()
        # plt.xlim([0, max(t_measured)])
        # plt.ylim([0, 1.1 * max(V_measured)])
        # plt.xlabel('Days Post Infection')
        # plt.ylabel('Virus Titre (Log10 copies/mL)')
        # plt.title('Subject ID=%i' %j)

        # display fitted statistics and append parameters to lists
        subj_IDs_over_5.append(j)
        report_fit(result)
        for name, param in result.params.items():
            if name == 'alpha':
                alphas.append(param.value)
            if name == 'beta':
                betas.append(param.value)

        red_chi_squs.append(result.redchi)
        residuals.append(result.residual)
        sum_residuals_squs.append(sum(((result.residual)**2)))
        chi_squs.append(result.chisqr)
        ndatas.append(result.ndata)

        #plot the fitted data and the model for log(virus) against day
        log_V_measured = np.log10(V_measured)
        log_I_fitted = np.log10(data_fitted[:, 1])

        ##########plot the log of virus amount against time
        # plt.figure()
        # plt.scatter(t_measured[1:], log_V_measured[1:], marker='o', color='red', label='measured V data', s=75) #the first point is found by extrapolation. Therefore it is not physical so dont plot it.
        # plt.plot(t_measured, log_I_fitted, '-', linewidth=2, color='red', label='fitted I data')
        #plt.ylim(bottom=0.9 * min(log_V_measured))
        #plt.xlim(left=0)
        # plt.legend()
        # plt.xlabel('Days Post Infection')
        # plt.ylabel('Concentration (Log10 copies/mL)')
        # plt.title('Subject ID=%i' %j)


print('alphas',alphas)
print('betas',betas)
#print('residuals',residuals)
print('sum_residuals_squs',sum_residuals_squs)
print('ndatas',ndatas)
variances = np.array(sum_residuals_squs) / np.array(ndatas)
print('variances',variances)
print('average variance',sum(variances)/len(variances))

#only include patients who have variance less than a value
refined_alphas = []
refined_betas = []

for i in range (len(variances)):
    if variances[i]<=100: #the value here is the cut off for the variance
        refined_alphas.append(alphas[i])
        refined_betas.append(betas[i])

########################### plot histograms of alpha, beta and kappa

n_bins = 50   #number of bins for the histogram

##############alphas

#plot the histogram of data
plt.figure()
plt.hist(refined_alphas, density=False, bins=n_bins,color = "skyblue",label="individual patient alpha values")
plt.ylabel('Density of patients')
plt.xlabel('Alpha')
plt.title('Histogram of alpha values across individual patients')

#plot the overall alpha (across all the patients) over the top
y, x, _ = plt.hist(refined_alphas, density=True, bins=n_bins,color = "skyblue")
X = [overall_alpha, overall_alpha]
Y = [0, y.max()]
plt.plot(X,Y,color='red',label="alpha value from model fit on average of patients")

# fit a lognormal distribution to the data
s, loc, scale = lognorm.fit(refined_alphas)
mu_alpha = loc
sigma_alpha = lognorm.std(s, loc=loc, scale=scale)
print('alpha s',s,'loc',loc,'scale',scale,'sigma_alpha',sigma_alpha)

#plot this lognormal distribution over the range of values required
x=np.linspace(0,np.amax(refined_alphas),500)
solu=lognorm.pdf(x, s, loc, scale)
plt.plot(x, solu,color="black",label="log normal of individual patient alpha values")

###plot the median of the alpha values over the top
median_alpha = statistics.median(refined_alphas)
print('median_alpha',median_alpha)
X = [median_alpha, median_alpha]
Y = [0, y.max()]
plt.plot(X,Y,color='green',label="median alpha value from histogram")
plt.legend()

#plot with labels on the points
fig, ax = plt.subplots()
ax.scatter(x,solu)
for i, txt in enumerate(x):
    ax.annotate(txt, (x[i], solu[i]))
#print('x',x,'len(x)',len(x),'solu',solu,'len(solu)',len(solu))

###########betas
plt.figure()
plt.hist(refined_betas, density=False, bins=n_bins,color = "skyblue",label="individual patient beta values")
plt.ylabel('Number of patients')
plt.xlabel('Beta')
plt.title('Histogram of beta values across individual patients')

#plot the overall beta (across all the patients) over the top
y, x, _ = plt.hist(refined_betas, density=True, bins=n_bins,color = "skyblue")
X = [overall_beta, overall_beta]
Y = [0, y.max()]
plt.plot(X,Y,color='red',label="beta value from model fit on average of patients")

# fit a histogram to the beta data

# best fit of data
s, loc, scale = lognorm.fit(refined_betas)
mu_beta = loc
sigma_beta = lognorm.std(s, loc=loc, scale=scale)
print('beta s',s,'loc',loc,'scale',scale,'sigma_beta',sigma_beta)

x=np.linspace(0,np.amax(refined_betas),500)
solu=lognorm.pdf(x, s, loc, scale)
plt.plot(x, 4*solu,color="black",label="log normal of individual patient beta values")

#plot the median of the beta values over the top
median_beta = statistics.median(refined_betas)
print('median_beta',median_beta)
X = [median_beta, median_beta]
Y = [0, y.max()]
plt.plot(X,Y,color='green',label="median beta value from histogram")
plt.legend()

##############################
#Compare gn(Ftrue) to the plot of the average over all patients

#firstly have at the optimal model
params = Parameters()
params.add('U0', value=U0, vary=False)
params.add('I0', value=I0, vary=False)

#use the mean values
params.add('alpha', value=mu_alpha, min=mu_alpha - 10**(-11), max=mu_alpha + 10**(-11))   #rate that viral particles infect susceptible cells
params.add('beta', value=mu_beta, min=mu_beta - 10**(-11), max=mu_beta + 10**(-11))
"""
# #use the median values
# params.add('alpha', value=median_alpha, min=median_alpha - 10**(-11), max=median_alpha + 10**(-11))   #rate that viral particles infect susceptible cells
# params.add('beta', value=median_beta, min=median_beta - 10**(-11), max=median_beta + 10**(-11))
# params.add('kappa', value=median_kappa, min=median_kappa - 10**(-11), max=median_kappa + 10**(-11))
"""
# fit model
result = minimize(residual, params, args=(t_measured_init, V_measured_init), method='leastsq')  # we put in the mean over all patients time and V_measured vals
report_fit(result)
data_fitted = g(t_measured_init, y0_init, result.params)
gn_Ftrue_log_I_fitted = np.log10(data_fitted[:, 1])

#plotting the model fit (found by taking the peak of the gaussian distribution of params for all patients)
plt.figure()
plt.plot(t_measured_init, gn_Ftrue_log_I_fitted, '-', linewidth=2, color='blue', label='Gn(Ftrue)')

#plot this with the scatterplot of the mean of the data for all patients
plt.scatter(t_measured_init[1:], log_V_measured_init[1:], marker='o', color='red', label='mean V data for all patients', s=75) #the first point is found by extrapolation. Therefore it is not physical so dont plot it.
plt.xlabel('Days Post Infection')
plt.ylabel('Virus Titre Concentration (Log10 copies/mL)')
plt.title('Scatterplot of average patient data with model fit from means of the lognormal distributions')
plt.legend()

###################################
#compute the posterior for validation patients

sn = 1 #for now. Compute better estimate later

######
#fit models to validation set of patients

##only want to model the last 3 patients as these are the validation data set
Subject_ID_vals_short = Subject_ID_vals[-3:]
print('Subject_ID_vals_short',Subject_ID_vals_short)

posteriors =[]

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

        # measured data
        t_measured = eff_day_list
        V_measured = act_div_vir_list_sum
        log_V_measured = np.log10(V_measured)

        #######compute first term of loss function

        print('subject ID',j,'t_measured',t_measured,'len(t_measured)',len(t_measured),'log_V_measured',log_V_measured,'len(log_V_measured)',len(log_V_measured))
        print('time points average over everyone',t_measured_init,'length',len(t_measured_init),'len(gn_Ftrue_log_I_fitted)',len(gn_Ftrue_log_I_fitted))

        indices = np.where(np.in1d(t_measured_init, t_measured))[0]
        print('indices',indices)

        #find differences between the gnFtrue virus amount and the virus amount for this patient
        diff = []
        for k in range (len(t_measured_init)):
            for m in range (len(t_measured)):
                if t_measured[m] == t_measured_init[k]:
                    diff.append(log_V_measured[m] - gn_Ftrue_log_I_fitted[k])
        print('diff',diff)

        #plot line plot of gnFtrue against t_measured_init in blue, and scatterplot of log_V_measured against t_measured in red, title of patient ID
        plt.figure()
        plt.title('Subject ID=%i' %j)
        plt.scatter(t_measured, log_V_measured, marker='o', color='red', label='patient V data', s=75)
        plt.plot(t_measured_init, gn_Ftrue_log_I_fitted, '-', linewidth=2, color='blue', label='Gn_Ftrue')
        plt.legend()
        plt.xlabel('Days Post Infection')
        plt.ylabel('Concentration (Log10 copies/mL)')

        ####################compute the likelihood term
        lik_term = -0.5*np.sum((np.square(diff)/(sn**2)) + np.log(2*np.pi*(sn**2)))
        print('lik_term',lik_term)

        ##################compute the prior term
        alpha_first_bit = 1/(((2*np.pi)**0.5)*(sigma_alpha**2))
        alpha_exponent = -((mu_alpha - overall_alpha)**2)/(2*(sigma_alpha**2))
        alpha_second_bit = np.exp(alpha_exponent)
        alpha_term = alpha_first_bit * alpha_second_bit

        beta_first_bit = 1/(((2*np.pi)**0.5)*(sigma_beta**2))
        beta_exponent = -((mu_beta - overall_beta)**2)/(2*(sigma_beta**2))
        beta_second_bit = np.exp(beta_exponent)
        beta_term = beta_first_bit * beta_second_bit

        prior = np.log(alpha_term * beta_term)
        print('prior',prior)

        ####################compute the posterior
        posterior = lik_term + prior
        posteriors.append(posterior)
        print('posterior',posterior)

print('posteriors',posteriors)
avg_posterior = np.sum(posteriors)/len(posteriors)
print('avg_posterior',avg_posterior)

plt.show()

