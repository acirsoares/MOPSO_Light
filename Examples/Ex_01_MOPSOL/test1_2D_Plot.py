import matplotlib.pyplot as plt
import numpy as np
import csv
from matplotlib import cm

def readMyFile(filename):  # Read y1,y2,x1,x2 from .csv file
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

# Plot y1 x y2 graph.

# Set the figure size
plt.rcParams["figure.figsize"] = [12.0, 5.0]
plt.rcParams["figure.autolayout"] = True

for i in range (1,11):

    # Set data
    y1,y2,x1,x2=readMyFile('MOPSOL_test1_%.2d.csv' %i)
    y1 = [float(i) for i in y1]
    y2 = [float(i) for i in y2]
    x1 = [float(i) for i in x1]
    x2 = [float(i) for i in x2]

    # Set point size and color
    datalen = len(y1)
    pointsize = 6
    ps = pointsize*np.ones(datalen) 

    fig1 = plt.figure("Figure 1")
    plt.subplot(2,5,i)
    plt.scatter(y1, y2,c="black",s=ps,marker='.',label="RUN_%.2d" %i);
    plt.legend(loc='upper right')
    plt.xlabel('OF_01')
    plt.ylabel('OF_02')

    # Set point size and color
    datalen = len(x1)
    pointsize = 6
    ps = pointsize*np.ones(datalen) 
    fig2 = plt.figure("Figure 2")

    plt.subplot(2,5,i)
    plt.scatter(x1, x2,c="blue",s=ps,marker='.',label="RUN_%.2d" %i);
    plt.legend(loc='upper right')
    plt.xlabel('x1')
    plt.ylabel('x2')


y1,y2,x1,x2=readMyFile('MOPSOLight_final.csv')
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


fig3 = plt.figure("Figure 3")
plt.scatter(y1, y2,c="black",s=ps,marker='.',label="Merged Pareto Front");
plt.legend(loc='upper right')
plt.xlabel('OF 1')
plt.ylabel('OF 2')


fig4 = plt.figure("Figure 4")
plt.scatter(x1, x2,c="blue",s=ps,marker='.',label="Merged Pareto Front");
plt.legend(loc='upper right')
plt.xlabel('x1')
plt.ylabel('x2')

    
plt.show()




