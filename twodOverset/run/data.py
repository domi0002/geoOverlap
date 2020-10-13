#!/usr/bin/python
import matplotlib.pyplot as pt
import numpy as np
from sklearn.linear_model import LinearRegression

meshCount=np.array([1352, 5408, 21632, 86528])
h=1/np.sqrt(meshCount)

#Create a linear fit( https://scipy-lectures.org/packages/scikit-learn/auto_examples/plot_linear_regression.html)
model       =   LinearRegression()

# Inverse Distance Results
l2_inv      =   np.array([0.069779, 0.02774, 0.01213, 0.00925])
model.fit(np.log10(h).reshape(-1,1), np.log10(l2_inv).reshape(-1,1))
l2_inv_lin  =   model.predict(np.log10(h).reshape(-1,1))
ax          =   pt.axes()
ax.scatter(np.log10(h),np.log10(l2_inv))
ax.plot(np.log10(h),l2_inv_lin, label= 'Slope (basic)= %.2f' %(model.coef_))
ax.set_xlabel('Log(Mesh spacing)')
ax.set_ylabel('Log(Error)')


#C_inv       =   np.array([0.0016214, 0.00070311, 0.0004368, 0.0002089, 9.0561e-5])
l2_inv_cons =   np.array([0.015296,  0.00696, 0.002179, 0.00147])
#C_inv_cons  =   np.array([1.3877e-17, 9.2898e-16, 1.51e-15, 8.21e-15,  3.7e-14])
model.fit(np.log10(h).reshape(-1,1), np.log10(l2_inv_cons).reshape(-1,1))
l2_inv_cons_lin  =   model.predict(np.log10(h).reshape(-1,1))
ax.scatter(np.log10(h),np.log10(l2_inv_cons))
ax.plot(np.log10(h),l2_inv_cons_lin, label= 'Slope (conservative)= %.2f' %(model.coef_))
pt.legend()
pt.show() 



# PLS Results
l2_pls      =   np.array([0.00649833, 0.0014804, 0.00103716, 0.000164582])
#C_pls       =   np.array([7.87e-5, 8.853e-6, 1.7023e-6, 1.595e-7,  1.738e-8])
model.fit(np.log10(h).reshape(-1,1), np.log10(l2_pls).reshape(-1,1))
l2_pls_lin  =   model.predict(np.log10(h).reshape(-1,1))
ax          =   pt.axes()
ax.scatter(np.log10(h),np.log10(l2_pls))
ax.plot(np.log10(h),l2_pls_lin, label= 'Slope (basic)= %.2f' %(model.coef_))
ax.set_xlabel('Log(Mesh spacing)')
ax.set_ylabel('Log(Error)')



l2_pls_cons =   np.array([0.002814, 0.000652, 0.0002284, 5.41886e-5])
#C_pls_cons  =   np.array([5.689e-16, 8.88e-16, 2.91e-16, 7.785e-15, 2.811e-14])
model.fit(np.log10(h).reshape(-1,1), np.log10(l2_pls_cons).reshape(-1,1))
l2_pls_cons_lin  =   model.predict(np.log10(h).reshape(-1,1))
ax.scatter(np.log10(h),np.log10(l2_pls_cons))
ax.plot(np.log10(h),l2_pls_cons_lin, label= 'Slope (conservative)= %.2f' %(model.coef_))
pt.legend()
pt.show() 

 

