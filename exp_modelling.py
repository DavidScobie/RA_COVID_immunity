import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# function that returns dy/dt
def model(y,t):
    k = 0.3
    dydt = -k * y
    return dydt

# initial condition
y0 = 5

# time points
t = np.linspace(0,20)

# solve ODE
y = odeint(model,y0,t)

# plot results
plt.figure(0)
plt.plot(t,y)
plt.xlabel('time')
plt.ylabel('y(t)')

###############################

# function that returns values
def new_model(z,t):
    I = z[0]
    V = z[1]
    U = z[2]
    alpha = 4*(10**(-8))
    beta = 1.07
    gamma = 3.07
    delta = 2.4
    dIdt = (alpha*V*U) - (beta*I)
    dVdt = (gamma*I) - (delta*V)
    dUdt = -alpha*V*U
    return [dIdt,dVdt,dUdt]

# initial conditions
z0 = [0,0.31,4*(10**(8))]

#time points
t = np.linspace(0,16)

#solve ODE
z = odeint(new_model,z0,t)

I = z[:,0]
V = z[:,1]
U = z[:,2]

#find logs of the results
logI = np.log10(I)
logV = np.log10(V)
logU = np.log10(U)

#plot results
plt.figure(1)
plt.plot(t,logI,'b--')
plt.plot(t,logV,'r-')
plt.plot(t,logU,'g--')
plt.xlabel('time')
plt.ylabel('log10(Virus)')

plt.show()