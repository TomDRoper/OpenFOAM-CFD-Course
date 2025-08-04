#Header
#---------------------------------------------------
# Python Introduction 
# Code Created By Colin Bailey
#---------------------------------------------------

#Section 0 - Workbooks
#---------------------------------------------------

#Imports workbooks into python to be used later
import numpy as np
import scipy as sp
from scipy import optimize
import time
from matplotlib import pyplot as plt



#Section 1 - Varibles, Operations and Groups
#---------------------------------------------------
# ================================
# Variable Overwriting in Python
# ================================

# Variable initially set as a string
x = 'cat'

# Reassigned to an integer
x = 1

# Reassigned again to a float
x = 2.0
# At this point, x is a float, so all operations below are based on float math

# ================================
# Simple Math Operators in Python
# ================================

# All arithmetic operations are done using x = 2.0
a = x + x     # Addition → 2.0 + 2.0 = 4.0
d = x - x     # Subtraction → 2.0 - 2.0 = 0.0
h = x * x     # Multiplication → 2.0 * 2.0 = 4.0
w = x / x     # Division → 2.0 / 2.0 = 1.0
n = x ** x    # Exponentiation → 2.0^2.0 = 4.0

# Print results for demonstration
print("Addition:", a)
print("Subtraction:", d)
print("Multiplication:", h)
print("Division:", w)
print("Exponentiation:", n)

# ================================
# Lists in Python
# ================================

# Lists can contain multiple data types
y = ['cat', 1, 2.0]

# Lists are mutable (you can change elements)
y[0] = 'dog'  # Replaces 'cat' with 'dog'
print("Modified list:", y)


# ================================
# Tuples in Python
# ================================

# Tuples can also store multiple data types
g = ('cat', 1, 2.0)

g[0] = 'dog'  




#Section 2 - Arrays and Loops
#----------------------------------------------------
# ===================================
# Array Creation Examples
# ===================================

# 1D array
z = np.array([1, 2, 3])

# 2D array (2 rows, 3 columns)
q = np.array([[1, 2, 3],
              [4, 5, 6]])

# 1D array of 25 zeros
o = np.zeros(25)

# 1D array of 25 ones
f = np.ones(25)

# 2D array (25x25) of ones with integer data type
t = np.ones((25, 25), dtype=int)

# Generate 5 evenly spaced values between 2.0 and 3.0
r = np.linspace(2.0, 3.0, num=5)
print("Linspace result:", r)


# ===================================
# For Loop Example with Timer
# ===================================

# Start timing
start = time.process_time()

# Initialize array
p = np.zeros(1)

# Loop from 1 to 3 (inclusive)
for k in range(1, 4):
    print("Current value of p:", p)
    # Uncomment below to actually append values to the array
    p = np.append(p, k)
    # p = np.append(p, k**2)

# End timing
end = time.process_time()
print("For loop elapsed time:", end - start)


# ===================================
# While Loop Example
# ===================================

r = 10  # Constant multiplier
s = 1   # Loop counter

# Loop while s is less than or equal to r
while s <= r:
    print(f"{s} * {r} = {s*r}")
    s += 1

# Optional if/else structure to notify when loop ends (uncomment to use)
# if s > r:
#     print('While loop completed')
# else:
#     print('While loop still going')


#Section 3 - Functions
#-------------------------------------------------
# ================================
# Example 1: Create your own function
# ================================

# Define an input value
o = 2

# Define a custom function that returns o^2 + o + 5
def func(o):
    return o**2 + o + 5

# Call the function with o = 2
t = func(o)

# Print the result
print("Result of func(o):", t)  # Output: 11


# ================================
# Example 2: Using scipy.optimize.minimize
# ================================

# Define a function of multiple variables to minimize
def fun(params):
    # Unpack the parameter list into named variables for clarity
    x, y, z = params
    
    # Return the value of the objective function: x² + y³ + z³
    return x**2 + y**3 + z**3

# Provide an initial guess for x, y, and z
first_guess = [0.5, 0.5, 0.5]

# Use scipy's minimize function to find the values of x, y, z that minimize `fun`
res = optimize.minimize(fun, first_guess)

# Output the optimized parameters
print("Optimized parameters:", res.x)


#Section 4 - Graphing
#-------------------------------------------------
# ==== User Adjustable Parameters ====
m = 3              # Slope
b = -1              # Y-intercept
x_min = -10        # Minimum x-value
x_max = 10         # Maximum x-value
num_points = 100   # Number of points in the graph

line_thickness = 5       # Thickness of the plotted line
line_style = '-'         # Style of the line ('-', '--', '-.', ':')
line_color = 'blue'      # Color of the line

font_size = 14           # Font size for labels and title
title_text = "Graph of y = mx + b"
show_grid = True         # Toggle to show grid (True/False)

# ==== Generate Data ====
x = np.linspace(x_min, x_max, num_points)
y = m * x + b

# ==== Plot ====
plt.figure(figsize=(8, 6))
plt.plot(x, y, line_style, linewidth=line_thickness, color=line_color)

# Add axes labels and title
plt.xlabel('x-axis', fontsize=font_size)
plt.ylabel('y-axis', fontsize=font_size)
plt.title(title_text, fontsize=font_size + 2)

# Show or hide grid
if show_grid:
    plt.grid(True)

# Optional: Add annotation for the equation
# eq_text = f"y = {m}x + {b}"
# plt.text(0.05, 0.9, eq_text, transform=plt.gca().transAxes,
#          fontsize=font_size, bbox=dict(facecolor='white', edgecolor='gray'))

# Show the plot
plt.show()

