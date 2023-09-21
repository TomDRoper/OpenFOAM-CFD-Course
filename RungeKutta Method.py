# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 08:52:15 2021

@author: baileycd2
"""
#this code will show how the Runge Kutta method is used in python
import numpy as np
import matplotlib.pyplot as plt

def RungeKutta(f,x0,t):
    x = np.zeros(len(t)) #len(t) gives the length of the variable t
    x[0] = x0
    for n in range (0,len(t)-1):
        e1 = h*f(x[n],t[n]) #defines all of the parts of the runge kutta
        e2 = h*f(x[n]+ 0.5*h,t[n] + 0.5*e1)
        e3 = h*f(x[n]+ 0.5*h,t[n] + 0.5*e2)
        e4 = h*f(x[n] + h,t[n] + e3)
        x[n+1] = x[n] + ((e1 +(2*e2)+(2*e3)+e4)/6) #puts them all together
    return x

t0 = 0 #intial time point
tf = 5 #final time point
n = 0.001 #used to change the time step
h = 10 ** (-n) # this creates the time step
print ('Time Step = ', h)


t = np.arange(t0,tf+h,h) #creates an array of the time step
x0 = 1 #intial value
#lambda is an anyomous function, this means it behaves like def but does not have a name in the variable explorer window
f = lambda x, t: -x #defines f for the function, this is used to put back into our RungeKutta function above
x = RungeKutta(f,x0,t) #runs the rungekutta function
x_true = np.exp(-t) #has the true value of the function that the runge kutta will be compared to

E = abs(x - x_true) #error between them
print(E) #This would print the difference in the values of the functions

Ex = np.log(E)

tx = np.log(t)

#this will plot the true value versus the runge kutta and you can see the difference
plt.plot(t,x,'b.-',t,x_true,'r-')
# plt.plot(tx,Ex)
# plt.plot(t,E)
plt.legend(['RungeKutta','True'])
plt.grid(True)
plt.title("Solution of $x'=-x, x(0)=1$, Time step = "+str(h))
plt.xlabel ('Time')
plt.ylabel ('Value')
plt.show()
