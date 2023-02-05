import matplotlib.pyplot as plt
import numpy as np
import math
import csv
from matplotlib import cm

# Subroutine to read Pareto Front data in file.
# Read y1,y2,x1,x2 from .csv file format.
def read_Pareto_Front_In_File(filename):  
    y1 = []
    y2 = []
    x1 = []
    x2 = []
    
    with open(filename) as csvDataFile:
        csvReader = csv.reader(csvDataFile,delimiter=';')
        for row in csvReader:
            y1.append(row[0])
            y2.append(row[1])
            x1.append(row[2])
            x2.append(row[3])
    return y1,y2,x1,x2

# Main : Plot Pareto Front

# Set path
filepath="results/"

# Set the figure size
plt.rcParams["figure.figsize"] = [12.0, 5.0]
plt.rcParams["figure.autolayout"] = True

# Set the Pareto Front for OF_test01
n=100
sx1 = np.linspace(0.1,1.0,n)
sx2 = 0.2*np.ones(len(sx1))
sy1 = sx1
sy2 = np.ones(len(sx1))
for i in range (0,n):
    g = 2.0 -math.exp(-((sx2[i]-0.2)/0.004)**2)-0.8*math.exp(-((sx2[i]-0.6)/0.4)**2)
    sy2[i] = g/sx1[i]

# Plot all runs
for i in range (1,11):

    # Set data for run "i"
    y1,y2,x1,x2=read_Pareto_Front_In_File(filepath+'MOPSOL_test1_%.2d.csv' %i)
    y1 = [float(i) for i in y1]
    y2 = [float(i) for i in y2]
    x1 = [float(i) for i in x1]
    x2 = [float(i) for i in x2]

    # Set point size and color
    datalen = len(y1)
    pointsize = 6
    ps = pointsize*np.ones(datalen) 

    # Set figure and Plot y1 x y2 graph for run "i".
    fig1 = plt.figure("Pareto Front (10 runs)")
    plt.subplot(2,5,i)
    plt.plot(sy1, sy2, c="green",linewidth=0.5,label="Pareto Front")
    plt.scatter(y1, y2,c="black",s=ps,marker='.',label="RUN_%.2d" %i);
    plt.legend(loc='upper right')
    plt.xlabel('OF_01')
    plt.ylabel('OF_02')

    # Set point size and color
    datalen = len(x1)
    pointsize = 6
    ps = pointsize*np.ones(datalen) 

    # Set figure and Plot x1 x x2 graph for run "i".
    fig2 = plt.figure("Domain of Pareto Front (10 runs)")
    plt.subplot(2,5,i)
    plt.scatter(x1, x2,c="blue",s=ps,marker='.',label="RUN_%.2d" %i);
    plt.legend(loc='upper right')
    plt.xlabel('x1')
    plt.ylabel('x2')

# Plot merged Pareto Front
# Set data
y1,y2,x1,x2=read_Pareto_Front_In_File(filepath+'MOPSOLight_final.csv')
y1 = [float(i) for i in y1]
y2 = [float(i) for i in y2]
x1 = [float(i) for i in x1]
x2 = [float(i) for i in x2]

# Set point size and color
datalen = len(y1)
pointsize = 6
ps = pointsize*np.ones(datalen) 
    
# Set the figure size
plt.rcParams["figure.figsize"] = [4.5, 4.0]
plt.rcParams["figure.autolayout"] = True

# Set figure and Plot y1 X y2 graph for merged Pareto Front.
fig3 = plt.figure("Merged Pareto Front (10 runs)")
plt.plot(sy1, sy2, c="green",linewidth=0.5,label="Pareto Front")
plt.scatter(y1, y2,c="black",s=ps,marker='.',label="Merged Pareto Front");
plt.legend(loc='upper right')
plt.xlabel('OF 1')
plt.ylabel('OF 2')

# Set figure and Plot x1 X x2 graph for domain of merged Pareto Front.
fig4 = plt.figure("Domain of merged Pareto Front (10 runs)")
plt.scatter(x1, x2,c="blue",s=ps,marker='.',label="Merged Pareto Front");
plt.legend(loc='upper right')
plt.xlabel('x1')
plt.ylabel('x2')
    
plt.show()





