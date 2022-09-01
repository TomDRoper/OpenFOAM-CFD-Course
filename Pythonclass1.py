# -*- coding: utf-8 -*-
"""
Spyder Editor

"""
#Imports workbooks into python to be used later
import numpy as np
import scipy as sp
from scipy import optimize
import time
from matplotlib import pyplot as plt

x = 'cat'
x = 1
x = 2.0
a = x + x
d = x - x
h = x * x 
w = x/x
n = x**x

#list
y = ['cat', 1 , 2.0]
y[0] = 'dog'

#tuple
g = ('cat', 1, 2.0)
#g[0] = 'dog'

#Array
z = np.array([1,2,3])
#Adding columns to array
q = np.array([[1,2,3],[4,5,6]])

o = np.zeros(25)
f = np.ones(25)
t = np.ones([25,25],dtype=int)

#add np.linspace example

#For loop
start = time.process_time()
p = np.zeros(1)
for k in range(1,4):
   #p = np.append(p,k)
   #p = np.append(p,k**2)
   print(p)

end = time.process_time()
print(end-start)

#while Loop
r = 10
s = 1
while s <= r:
    print(s*r)
    s = s + 1
    #if loop
    # if s > r:
    #     print('while loop completed')
    # else:
    #     print('while loop still going')

#Create your own function in python
o = 2
def func(o):
   return o**2 + o + 5 
    
t = func(o)
print(t)

#Scipy has a minimize function that may be useful in the future
def fun(paramt):
    # print(paramt)  # As you can see, params is an array in NumPy.
    x, y, z = paramt # You might want to give the component variables names to make them easier to read.
    return x**2 + y**3 + z**3

first_guess = [0.5, 0.5, 0.5]
res = optimize.minimize(fun, first_guess)

res.x

#How to graph
# x axis values 
x = [1,2,3] 
# corresponding y axis values 
y = [2,4,1] 
    
# plotting the points  
plt.plot(x, y) 
    
# naming the x axis 
plt.xlabel('x - axis') 
# naming the y axis 
plt.ylabel('y - axis') 
    
# giving a title to my graph 
plt.title('The graph works!') 
    
# function to show the plot 
plt.show() 
        