import math
import numpy as np
from matplotlib import pyplot as plt

separation=0.711 #kpc
cutoff=0.250 #kpc

view_mag=4
theta_exact=math.asin(cutoff/separation)
print(theta_exact * 180 / np.pi)

rad=math.tan(theta_exact)*view_mag

r=view_mag*math.sin(theta_exact)

displacement=view_mag*math.cos(theta_exact)

r1=rad
r2=rad

coverage=0.5*math.tan(theta_exact)**2
print ("Fraction of incorrect identification is: %.5f" %(coverage))
interval=np.linspace(0,(2*math.pi),num=200)
z1=np.multiply(np.cos(interval),r1)
y1=np.multiply(np.sin(interval),r1)
zeros=np.zeros(len(z1))
zeros2=np.full(len(z1), separation)
dis=np.full(len(z1), displacement)
dis2=np.full(len(z1), -displacement)
plt.rcParams["font.family"] = "serif"
font = {'fontname':'serif'} 
plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True
fig = plt.figure()
ax = fig.add_subplot(projection="3d")
x, y, z = [0, separation], [0, 0], [0, 0]
ax.scatter(x, y, z, c='red', s=100)
ax.plot(x, y, z, color='black')
ax.plot(dis,z1,y1, c='cyan')
ax.plot(dis2,z1,y1, c='cyan')
ax.set_xlabel("X-Direction (Mpc)", **font)
ax.set_ylabel("Y-Direction (Mpc)", **font)
ax.set_zlabel("Z-Direction (Mpc)", **font)
ax.axis('equal')
i=0
while i <= (len(interval)-1):
    x,y,z=[dis[i],separation],[y1[i],0],[z1[i],0]
    ax.plot(x,y,z, alpha=0.05, c='cyan')
    x,y,z=[dis[i],0],[y1[i],0],[z1[i],0]
    ax.plot(x,y,z, alpha=0.05, c='cyan')
    x,y,z=[dis2[i],0],[y1[i],0],[z1[i],0]
    ax.plot(x,y,z, alpha=0.05, c='cyan')
    x,y,z=[dis2[i],separation],[y1[i],0],[z1[i],0]
    ax.plot(x,y,z, alpha=0.05, c='cyan')
    i+=1


 # draw sphere
u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:50j]
x = view_mag*np.cos(u)*np.sin(v)
y = view_mag*np.sin(u)*np.sin(v)
z = view_mag*np.cos(v)
ax.view_init(elev=12, azim=45)
ax.plot_surface(x, y, z, color='grey', alpha=0.1)
plt.savefig("plots/wrong_neighbour", dpi=300)
plt.show()



