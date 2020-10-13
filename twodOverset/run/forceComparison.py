#!/usr/bin/python
import matplotlib.pyplot as pt
import numpy as np
from sklearn.linear_model import LinearRegression
import csv

result1 = np.genfromtxt('forces.c.dat',names=['time','error','ninterp'])
result2 = np.genfromtxt('forces.nc.dat',names=['time','error','ninterp'])
ax          =   pt.axes()
ax.plot(result1['time'],np.absolute(result1['error']), 'k+-',label='Conservative')
ax.plot(result2['time'],np.absolute(result2['error']), 'r*-',label='Standard')
ax.set_ylim(ymin=0,ymax=0.0125)


sax         =   ax.twinx()
sax.plot(result1['time'],np.absolute(result1['ninterp']), 'bo-',label= 'Interpolation Cell Count')
ax.set_xlabel(' Time (s)')
ax.set_ylabel('abs(Force Error)')
sax.set_ylabel('Interpolation Cell Count')
sax.set_ylim(ymax=1580)
ax.legend()

pt.legend(loc="upper left")
pt.show() 

