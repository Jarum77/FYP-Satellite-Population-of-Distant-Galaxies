import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data/MW SFR/LMC-SFR.txt", skiprows=1,  delimiter=',') 

age = data[:,0]
SFH = data[:,1]


plt.semilogx(age,SFH)
plt.show()