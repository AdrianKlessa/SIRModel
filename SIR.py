import numpy as np
import matplotlib.pyplot as plt

beta =0.5 #Transmission rate (function of the rate of contact and probability of transmission given contact)
gamma = 0.1 #recovery rate

time_period = 200 #days
dt = 0.1

periods = (int)(time_period/dt)

S = np.zeros(periods+1)
I = np.zeros(periods+1)
R = np.zeros(periods+1)

S[0]=1000
I[0]=1
#R[0] is already 0
N=S[0]+I[0]+R[0]

X = np.linspace(0,time_period,periods+1)

def dSdt(beta, infected, susceptible, overall):
    return (-1*((beta*infected*susceptible)/overall))

def dIdt(beta, infected, susceptible, overall, gamma):
    return (((beta*infected*susceptible)/overall)-(gamma*infected))

def dRdt(gamma, infected):
    return (gamma*infected)

for x in range(1, periods+1):
    S[x]=S[x-1]+dt*dSdt(beta,I[x-1],S[x-1],N)
    I[x]=I[x-1]+dt*dIdt(beta,I[x-1],S[x-1],N,gamma)
    R[x]=R[x-1]+dt*dRdt(gamma,I[x-1])

plt.plot(X,S, label="susceptible")
plt.plot(X,I, label="infected")
plt.plot(X,R, label="recovered")
#N is the overall population, so it has to be the sum of SIR since we don't model deaths
print("Most infected at day: ", np.argmax(I)*dt, "(",np.amax(I),")")
plt.scatter(np.argmax(I)*dt,np.amax(I), label="Max infected")
plt.legend()



plt.show()