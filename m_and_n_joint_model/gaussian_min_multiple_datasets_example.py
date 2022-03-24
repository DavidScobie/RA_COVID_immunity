import matplotlib.pyplot as plt
import numpy as np

from lmfit import Parameters, minimize, report_fit


def gauss(x, amp, cen, sigma): #amplitude, centre and sigma define the shape of the gaussian
    """Gaussian lineshape."""   #this function makes the gaussian
    return amp * np.exp(-(x-cen)**2 / (2.*sigma**2))


def gauss_dataset(params, i, x):
    """Calculate Gaussian lineshape from parameters for data set."""
    #this function takes the 3 parameters (with the correct names of 'amp_', 'cen_' and 'sig_) and assigns to amp, cen and sig
    #then it returns the gaussian at certain x coordiantes (from -1 to 2 in this example)
    amp = params[f'amp_{i+1}']
    cen = params[f'cen_{i+1}']
    sig = params[f'sig_{i+1}']
    return gauss(x, amp, cen, sig)


def objective(params, x, data):
    """Calculate total residual for fits of Gaussians to several data sets."""
    #objective function here is analagous to residual function in my code (in that it is first argument of minimize)
    ndata, _ = data.shape #this is just finding the height, ndata (5) and width, _ (151) of the shape of the data
    resid = 0.0*data[:]  #this is just making an empty array (of same shape as data) ready for loop below

    # make residual per data set
    for i in range(ndata): #range(5)
        #this is just subtracting the model fit away from the measured data
        resid[i, :] = data[i, :] - gauss_dataset(params, i, x)

    # now flatten this to a 1D array, as minimize() needs
    return resid.flatten()

np.random.seed(2021)  #this somehow keeps the measured data the exact same every time the code is run
x = np.linspace(-1, 2, 151) #x goes between -1 and 2, and has 151 points along the way
data = []
for _ in np.arange(5): #[0, 1, 2, 3, 4]
    #here we are
    params = Parameters()
    amp = 0.60 + 9.50*np.random.rand()
    print('amp',amp)
    cen = -0.20 + 1.20*np.random.rand()
    sig = 0.25 + 0.03*np.random.rand()
    dat = gauss(x, amp, cen, sig) + np.random.normal(size=x.size, scale=0.1)
    data.append(dat)
data = np.array(data)

fit_params = Parameters()
for iy, y in enumerate(data):
    fit_params.add(f'amp_{iy+1}', value=0.5, min=0.0, max=200)
    fit_params.add(f'cen_{iy+1}', value=0.4, min=-2.0, max=2.0)
    fit_params.add(f'sig_{iy+1}', value=0.3, min=0.01, max=3.0)

for iy in (2, 3, 4, 5):
    fit_params[f'sig_{iy}'].expr = 'sig_1'

out = minimize(objective, fit_params, args=(x, data))
report_fit(out.params)

plt.figure()
for i in range(5):
    y_fit = gauss_dataset(out.params, i, x)
    plt.plot(x, data[i, :], 'o', x, y_fit, '-')

plt.show()