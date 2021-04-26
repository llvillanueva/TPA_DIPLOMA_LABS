import numpy as np
import matplotlib.pyplot as plt
from celluloid import Camera

# parameters
xsize=100 # domain in m
dx=1.0 # resolution in m 
dt=0.5 # timestep in s
uvel=1.0 # horizontal velocity m/s
time=40. # total time

nsteps=int(time/dt)

# x domain
xpt=np.linspace(0,xsize,int(xsize/dx)+1)

# initial function sine:
yinit=np.sin(2*np.pi*xpt/(xsize+1))

# squareform initial function:
#yinit=(x>xsize/2)*1
y0=yinit

# dy/dt+u dy/dx=0

fig0,ax0=plt.subplots()
camera = Camera(fig0)

k0=dt*uvel/(2.*dx) # constant in time

# solve with forward in time and centered in space
for it in range(nsteps):
    y1=y0-k0*(np.roll(y0,-1)-np.roll(y0,1))
    y0=y1
    ytrue=np.roll(yinit,int((it+1)*dt*uvel))
    ax0.plot(xpt,y1,color="blue",label="Numerical FTCS")
    ax0.plot(xpt,ytrue,color="red",linestyle="-",label="True Solution")
    camera.snap()
    
animation = camera.animate()
animation.save('celluloid_minimal.gif', writer = 'imagemagick')
    
# true answer
    
fig1,ax1=plt.subplots()
ax1.plot(xpt,yinit,label="initial conditions")
ax1.plot(xpt,y1,label="Numerical FTCS")
ax1.plot(xpt,ytrue,label="True Solution")
ax1.set_xlabel("x pts (m)")
ax1.set_ylabel("y")
ax1.legend()
plt.savefig("fig.pdf")


