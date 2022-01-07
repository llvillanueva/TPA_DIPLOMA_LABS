#!/usr/bin/env python
# coding: utf-8

# $\frac{Dq}{Dt}= \frac{dq}{dt}+u\frac{dq}{dx}=S=0$
# 

# In[1]:


import numpy as np 
import matplotlib.pyplot as plt


# In[2]:


# parameter
xsize=100. #(m)
dx=1  # (m)

dt=0.05 # timestep (sec)
uvel=1.0  # horizontal velocity (m/s)
time=40. # total time of the simulation

# set up x grid and initial profile of q 
xpt=np.linspace(0,xsize,int(xsize/dx)+1)
qinit=np.sin(2*np.pi*xpt/(xsize+1))

#qinit=(xpt>xsize/2)*2-1
#qinit=(qinit+np.roll(qinit,1)+np.roll(qinit,-1))/3
print(qinit)


# Now we want to discretize the equation
# 
# Using forward in time for $\frac{dq}{dt}$ and centered in space difference for $\frac{dq}{dx}$
# 
# $\frac{q_{t+1,i}-q_{t,i}}{\Delta t} + u \frac{q_{t,i+1}-q_{t,i-1}}{2 \Delta x}=0$
# 
# giving us an expression for the future value 
# 
# $q_{t+1,i} = q_{t,i} - u\Delta t \frac{q_{t,i+1}-q_{t,i-1}}{2 \Delta x}$
# 
# 

# In[3]:


nsteps=int(time/dt)
q0=qinit # current value
K0=uvel*dt/(2*dx)

# integrate equations
for it in range(nsteps):
    q1=q0-K0*(np.roll(q0,-1)-np.roll(q0,1))
    q0=q1
    qtrue=np.roll(qinit,int((it+1)*dt*uvel))

fig, ax = plt.subplots()
ax.set_title("this is the title")
ax.plot(xpt,qinit,label="Initial")       
ax.plot(xpt,q1,label="Final")
ax.plot(xpt,qtrue,label="True")

ax.legend() 
plt.show()


# In[4]:


xtest=np.array(range(10))
xtest


# In[5]:


np.roll(xtest,-1)


# In[ ]:




