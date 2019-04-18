import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from numpy import *
from pylab import *
from matplotlib.patches import Circle
import matplotlib.animation as manimation
import matplotlib.patches as patches
import matplotlib.colors as colors
import time
import math


start_time = time.time()

print("--- %s seconds ---" % (time.time() - start_time))

parameters = np.genfromtxt(r'./parameters.txt', unpack=True)

print(parameters)
print(parameters[0])

nx = int(parameters[0])
ny = int(parameters[1])
timesteps = int(parameters[2])
max_velocity = parameters[3]
max_vorticity = parameters[4]
min_v_y = parameters[5]
max_v_y = parameters[6]

nx_low = 0
nx_high = nx
ny_low = 0
ny_high = ny

print(nx)
print(ny)

plt.rcParams['animation.ffmpeg_path']='C:/Users/d-w-h/Desktop/Home/Programs/ffmpeg-20190318-15d016b-win32-static/bin/ffmpeg.exe'
writer=manimation.FFMpegWriter(bitrate=20000, fps=15)

fig = plt.figure(figsize=(8,4))
ax = plt.axes(xlim=(nx_low, nx_high-1), ylim=(ny_low, ny_high-1))  
plt.xlabel(r'x')
plt.ylabel(r'y')

x_x, x_y = np.genfromtxt(r'./object_location.txt', unpack=True)
    
def animate(i):

    x, y, fx, fy, vort_field = np.genfromtxt(r'./simdata/data_' + str(i) + '.txt', unpack=True)

    print("--- %s seconds ---" % (time.time() - start_time))

    X = np.array(np.reshape(x, (1, ny, nx)))
    Y = np.array(np.reshape(y, (1, ny, nx)))
    U = np.array(np.reshape(fx, (1, ny, nx)))
    V = np.array(np.reshape(fy, (1, ny, nx)))
    VortField = np.array(np.reshape(vort_field, (1, ny, nx)))
    
    xlow = 0
    xhigh = nx
    ylow = 0
    yhigh = ny 

    Vfieldnorm = np.log(1 + (VortField**2)**(0.5))
    Velocity = ((U**2 + V**2)**(0.5))
    Velocitynorm = 1/(Velocity + 1e-16)

    print(i)
    plt.gcf().clear()
    ax = plt.axes(xlim=(nx_low + 1/2, nx_high - 1/2), ylim=(ny_low + 1/2, ny_high - 1/2))
    VField = Vfieldnorm[0,:,xlow:xhigh]
    z = Velocity[0,:,xlow:xhigh]
    x = (X[0,:,xlow:xhigh])
    y = (Y[0,:,xlow:xhigh])
    u = (U[0,:,xlow:xhigh])
    v = (V[0,:,xlow:xhigh])
    Vnorm = Velocitynorm[0,:,xlow:xhigh]
    
    #cont = plt.contourf(x, y, u, v, pivot='mid', angles='xy', width=0.0005, color='w')    
    #cont = plt.contourf(x, y, VField, linspace(0, max_vorticity, 250), cmap='RdBu_r')        

    strm = ax.streamplot(x, y, u, v, linewidth=2)

    cont = plt.contourf(x, y, v, linspace(-0.28, 0.1, 250), cmap='viridis') 
    circle2 = plt.Circle((x_x[i], x_y[i]), 6.2)   
    ax.add_artist(circle2)  
    cbar = plt.colorbar(cont)
        
    return circle2


print("Starting Animation")

anim = manimation.FuncAnimation(fig, animate, frames=timesteps-1, repeat=False)

print("Done Animation, start saving")

anim.save('Sim_Results.mp4',
          writer=writer, dpi=200)

print("--- %s seconds ---" % (time.time() - start_time))