import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn
from lmfit import minimize, Parameters, Parameter, report_fit
from scipy.integrate import odeint
from numpy import arange
from pandas import read_csv
from scipy.optimize import curve_fit

#plot the FFA data
FFA_virus = np.load('FFA_MTS_V_measured_NON_DET_eq_zero.npy')
print('FFA_virus',FFA_virus)
log_FFA_virus = np.log10(FFA_virus)
FFA_effective_day = np.load('FFA_MTS_t_measured_NON_DET_eq_zero.npy')
plt.figure()
plt.plot(FFA_effective_day, log_FFA_virus, marker='o', color='black', label='measured FFA V data')

#plot the qPCR data
log_qPCR_data = np.log10(np.load('qPCR_MTS_V_measured_NON_DET_eq_zero.npy'))[:-2]
qPCR_effective_day = np.load('qPCR_MTS_t_measured_NON_DET_eq_zero.npy')[:-2]
print('len(log_qPCR_data)',len(log_qPCR_data),'len(qPCR_effective_day)',len(qPCR_effective_day))
plt.plot(qPCR_effective_day, log_qPCR_data, marker='o', color='red', label='measured qPCR V data')
plt.ylim([0, 1.1 *max(log_qPCR_data)])
plt.xlim([0, 20])
plt.legend()

print('FFA_effective_day',FFA_effective_day,'length',len(FFA_effective_day))
print('log_FFA_virus',log_FFA_virus,'length',len(log_FFA_virus))
print('qPCR_effective_day',qPCR_effective_day,'length',len(qPCR_effective_day))
print('log_qPCR_data',log_qPCR_data,'length',len(log_qPCR_data))

#qPCR data is longer than FFA data, so we need to cut qPCR data off at the right points
qPCR_lower_bound = np.where(qPCR_effective_day == FFA_effective_day[0])
qPCR_upper_bound = np.where(qPCR_effective_day == FFA_effective_day[-1])
print('qPCR_lower_bound',int(qPCR_lower_bound[0]),'qPCR_upper_bound',int(qPCR_upper_bound[0]))
print('qPCR_effective_day',qPCR_effective_day)

#keep only the data for each type on the same days
log_qPCR_data_short = log_qPCR_data[int(qPCR_lower_bound[0]):int(qPCR_upper_bound[0])+1] #have to add 1 because of how python does indexing
qPCR_effective_day_short = qPCR_effective_day[int(qPCR_lower_bound[0]):int(qPCR_upper_bound[0])+1]
print('qPCR_effective_day_short',qPCR_effective_day_short)

print('log_FFA_virus',log_FFA_virus,'length',len(log_FFA_virus))
print('log_qPCR_data_short',log_qPCR_data_short,'length',len(log_qPCR_data_short))

fig, ax = plt.subplots()
ax.scatter(log_FFA_virus, log_qPCR_data_short)  #this is the scatterplot of the points through which the line of best fit is drawn

for i, txt in enumerate(FFA_effective_day):  #this is just labelling all the scatterpoints
    ax.annotate(txt, (log_FFA_virus[i], log_qPCR_data_short[i]))

ax.set_xlim(left=0)
ax.set_ylim(bottom=0)

#####################

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

a1, b1 = best_fit(log_FFA_virus, log_qPCR_data_short)

yfit = [a1 + b1 * xi for xi in log_FFA_virus]
print('yfit',yfit)

#add the point at time=0, virus=933 to the eff_day_vals and act_div_vir_list_sum arrays
v1 = 0
v2 = a1 #THIS IS 10**a. where a is the y intercept of the line of best fit
log_FFA_virus_inc0 = np.insert(log_FFA_virus, 0, v1, axis=0)
yfit = np.insert(yfit, 0, v2, axis=0)
#print('eff_day_vals',eff_day_vals,'act_div_vir_list_sum',act_div_vir_list_sum)

plt.plot(log_FFA_virus_inc0, yfit)
plt.xlabel('log_FFA_virus')
plt.ylabel('log_qPCR_data_short')

##########################
new_FFA_log_virus = a1 + (b1*log_FFA_virus)

print('FFA_effective_day',FFA_effective_day)
print('new_FFA_log_virus',new_FFA_log_virus)

plt.figure()
plt.plot(FFA_effective_day, new_FFA_log_virus, marker='o', color='blue', label='transformed FFA V data')
plt.plot(qPCR_effective_day, log_qPCR_data, marker='o', color='red', label='measured qPCR V data')
plt.plot(FFA_effective_day, log_FFA_virus, marker='o', color='black', label='measured FFA V data')
plt.ylim([0, 1.1 *max(log_qPCR_data_short)])
plt.xlim([0, 20])
plt.xlabel('Effective Day')
plt.ylabel('Virus Titre Concentration (Log10 copies/ml)')
plt.legend()

#######scatterplot of the absolute values of the qPCR and the FFA data

#then make a scatterplot of the FFA and the qPCR data
qPCR_data = 10**log_qPCR_data

#firstly, check that the plots of FFA and qPCR virus against day look as they should
plt.figure()
plt.plot(FFA_effective_day, FFA_virus, marker='o', color='black', label='measured FFA V data')
plt.xlabel('effective day')
plt.ylabel('Virus Titre Concentration (million copies/ml)')
plt.figure()
plt.plot(qPCR_effective_day, qPCR_data, marker='o', color='red', label='measured qPCR V data')
plt.xlabel('effective day')
plt.ylabel('Virus Titre Concentration (million copies/ml)')

qPCR_data_short = qPCR_data[int(qPCR_lower_bound[0]):int(qPCR_upper_bound[0])+1] #have to add 1 because of how python does indexing

fig, ax = plt.subplots()
ax.scatter(FFA_virus, qPCR_data_short)  #this is the scatterplot of the points through which the line of best fit is drawn

for i, txt in enumerate(FFA_effective_day):  #this is just labelling all the scatterpoints
    ax.annotate(txt, (FFA_virus[i], qPCR_data_short[i]))

ax.set_xlim(left=0)
ax.set_ylim(bottom=0)

plt.xlabel('FFA_virus')
plt.ylabel('qPCR_data_short')

###############
#regression to find relationship between the FFA and qPCR data
# fit a fifth degree polynomial to the data

# define the true objective function
def objective(x, a, b, c, d, e, f):
	return (a * x) + (b * x**2) + (c * x**3) + (d * x**4) + (e * x**5) + f

# choose the input and output variables
x, y = FFA_virus, qPCR_data_short
# curve fit
popt, _ = curve_fit(objective, x, y)
# summarize the parameter values
a, b, c, d, e, f = popt
print('y = %.5f * x + %.5f * x^2 + %.5f * x^3 + %.5f * x^4 + %.5f * x^5 + %.5f' % (a, b, c, d, e, f))
# plot input vs output
plt.scatter(x, y)
# define a sequence of inputs between the smallest and largest known inputs
x_line = arange(min(x), max(x), 1)
# calculate the output for the range
y_line = objective(x_line, a, b, c, d, e, f)
# create a line plot for the mapping function
plt.plot(x_line, y_line, '--', color='red')

#### plot the log proportionality that we though the relationship was previously
y_log_fit = 10**(a1)*((x_line**b1))
plt.plot(x_line, y_log_fit, '--', color='blue')

plt.show()