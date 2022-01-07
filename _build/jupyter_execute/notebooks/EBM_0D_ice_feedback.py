#!/usr/bin/env python
# coding: utf-8

# # Ice albedo feedback in the 0D energy balance model
# 
# **Credits: This notebook is based on a class exercise provided by Simona Bordoni (University of Trento) with additional ideas taken from lesson 6 of Brian Rose's climate modelling course at Albany.**

# We will start our consideration of climate feedbacks with one that is not a particularly fast atmosphere process, but that is easy to understand to introduce the concept.  While glaciers have a long inertia in terms of their melting timescale, the change in snow cover on land and the growth and melting of sea ice, can occur on timescales that are seasonal and relatively fast compared to the forcing of climate change.  Take a look at this figure which shows the extent of artic sea ice cover in September (the annual relative minimum) in 1984, not long after the modern record of remotely sensed data began and 2012, which was a record mimimum until it was broken in 2020/21.
# 
# <div>
#     <img src="attachment:esm_ebm_Arctic_Sea_Ice_Minimum_Comparison.png" width="300"/>
# </div>
# Caption: Ice cover comparison 1984 and 2012 [credit: https://energyeducation.ca/]
# 
# 
# ## Question:  How do you think the shrinking ice cover with change the TOA radiative balance?

# Consider the globally averaged energy balance equation for a planet whose surface temperature is $T_s$, this time no equilibrium is assumed, so the difference between the incoming and outgoing radiation leads to a temperature tendancy:
#  
# $$
#   C \frac{dT_s}{dt} = \frac{S_0}{4} [ 1 - \alpha(T_s)]  - \epsilon_e \sigma T_s^4
# $$
#  
# where $C$ is the heat capacity of our system (units J m$^{-2}$ K$^{-1}$), related mainly to the deep ocean (we come back to this later), and $\epsilon_e$ is the effective surface emissivity of the Earth (refer to lecture notes). 
#  
# Here the albedo of the planet $\alpha$ is now a function of global mean temperature rather than a fixed constant.  
# Why is this?  Well, this is to mimic the feedback due to ice (although of course one could also try to represent other temperature sensitive albedo feedbacks such as land-surface cover, although difficult to do in such a simple model).  
# 
# As the temperature cools, more of the planet is assumed to be below the threshold for ice cover, increasing the global mean albedo.  This is represented by the following function:  
# 
# $$
#   \alpha = 0.45 - 0.25 \,\mathrm{tanh} \left( \frac{(T-272)}{23} \right)
# $$
#  
#  Note that we don't actually know what the ice cover is in this model, we only know the final impact on albedo.  
#  
#  First of all we will define this function of albedo:

# In[1]:


import matplotlib.pyplot as plt
import numpy as np

# function to calculate alfa
def alpha(T):
    """ function for albedo"""
    return 0.45-0.25*np.tanh((T-272)/23)


# Let's take a look at this function by plotting it.  We can do this by defining a vector of temperatures...

# In[2]:


# vector of surface temperatures
Ts=np.arange(200,350,0.5)
print (Ts)


# In[3]:


# vector of albedos
fig,ax=plt.subplots()
ax.plot(Ts,alpha(Ts))
ax.set_xlabel("Ts (K)")
ax.set_ylabel("Albedo")


# So we can see that for temperatures exceeding 300K roughly, the Earth becomes ice-free and the albedo tends to a fixed values of about 0.2.  As temperature cools the ice edge moves to the south, until for conditions below Ts=225K most of the globe is frozen and the albedo tends to 0.7.  
# 
# So now, rather than solve the equation, we are simply going to plot the RHS to get the sign of DT/Dt
# 

# In[4]:


# constants needed
sigma=5.67e-8 # Stefan-B constant 
eps_e=0.61    # equivalent emissivity of surface: 1-eps/2

# define the solar constant array
S0today=1370 

# define the energy balance equation in a function
def cdTdt(S,T):
    """Temperature tendancy=F(S,T), S=Solar constant, T=temperature"""
    return S*(1-alpha(T))/4 - eps_e*sigma*np.power(T,4)

# limits of temperature array we will need later:
tmin,tmax=np.min(Ts),np.max(Ts)

# plot
fig,ax=plt.subplots()
ax.plot(Ts,cdTdt(S0today,Ts))
ax.set_xlabel("Ts (K)")
ax.set_ylabel("C DT/Dt (W m$^{-2}$)")
ax.hlines(0,tmin,tmax,linestyles="dotted")


# ## Question: Where are the three equilibria temperatures and are they stable or unstable?
# 
# When you have discussed this, we will go on to plot a contour plot of $C \frac{DT}{Dt}$ for a range of values for S0 and Ts, with S0 ranging  
# 

# In[5]:


# define the solar constant array
S0=np.arange(1000,2000,10)

# make 2D arrays for contour plot:
T2d=np.tile(Ts,(len(S0),1)).transpose()
S2d=np.tile(S0,(len(Ts),1))


# In[6]:


# contour plot
fig,ax=plt.subplots()
X=ax.contour(S0,Ts,cdTdt(S2d,T2d),levels=np.arange(-200,400,25))
ax.clabel(X,fontsize=10)
ax.set_xlabel("S0 (Wm$^{-2}$)")
ax.set_ylabel("T (K)")
ax.vlines(S0today,tmin,tmax,linestyles="dotted")


# ## Question: 
# What happens to temperature if you reduce S0 from the present day value to 1000 Wm$^{-2}$?  If you then increase S0 to 1600 W m$^{-2}$, what happens?  Does the temperature trace the same path?
# 
# ## Exercise
# 
# Write a small python code that iterates the temperature until $\frac{Dt}{dt}=0$.  Now do the following...  Starting from $S_0=1340 W m^{-2}$ calculate $T_s$ (use a first guess $T_s$=320K.  Then reduce $S_0$ by 20 W m$^{-2}$, find the new equilibrium (starting the iteration from the previous equilibrium value), repeat until $S_0$=1000 W m$^{-2}$, storing the equilibrium . Then reverse the process by increasing $S_0$ in 20 $W m^{-2}$ steps back to the present day and beyond to 1600 W m$^{-2}$.  Plot $T_s$ against $S_0$.  
# 
# 

# 
# 
# ## Question: 
# We will now calculate the ice albedo feedback parameter.
# 
# Similarly as we did for the Planck feedback, we need to calculate the rate of change of flux imbalance with respect to temperature:
# 
# $$
# -\frac{\partial N}{\partial \alpha}\frac{\partial \alpha}{\partial T}=-\underbrace{\frac{-S0}{4}}_{\frac{\partial N}{\partial \alpha}}\underbrace{\frac{1}{4}\frac{1}{23}sech^2\frac{T-272}{23}}_{\frac{\partial \alpha}{\partial T}}
# $$
# 
# which uses the trig relation $\frac{d tanh(x)}{dx}=sech^2(x)$.
# 
# We can now code this up to calculate the feedback, we need the relationship for $sech^2(x)$:
# 
# $sech^2(x)=\left(\frac{2}{e^x+e^{-x}}\right)^2$
# 
# This gives
# 

# In[7]:


import math as mp
x=(Ts-272)/23

icefeedback=S0today*np.power(2.0/(np.exp(x)+np.exp(-x)),2)/(4*4*23)

fig,ax=plt.subplots()
ax.plot(Ts,icefeedback)
ax.set_ylabel("Albedo Feedback $\lambda_{ice}$ (Wm$^{-2}$K$^{-1}$)")
ax.set_xlabel("T (K)")


# Note that the magnitude of the feedback is strongly dependent on the temperature, i.e. the ice albedo feedback strength is not an absolute but depends on the current state.  
# 
# ## Question: Why is the feedback zero at cold or warm temperatures?
# 
# This makes sense if you think about it a little.  For example, at very cold temperatures when the world is ice covered, further decreases in temperature can not increase ice and do little to impact albedo, thus it is clear the feedback should be zero.  Likewise for warm temperatures when the world is ice free. 
# 
# ## Question: Is the overall feedback (ice+albedo) ever unstable?
# 
# Try to add the Planck feedback to the plot... 

# In[8]:


# add hline for magnitude of Planck
planckfeedback=-4.0*eps_e*sigma*np.power(Ts,3)

fig,ax=plt.subplots(nrows=2)
ax[0].plot(Ts,icefeedback,label="Ice feedback")
ax[0].plot(Ts,planckfeedback,label="Planck feedback")
ax[0].set_ylabel("$\lambda$ (Wm$^{-2}$K$^{-1}$)")
ax[0].set_xlabel("T (K)")
ax[0].legend()

feedback=icefeedback+planckfeedback
ax[1].plot(Ts,feedback)
ax[1].hlines(0.0,tmin,tmax,linestyles="dotted")
ax[1].set_xlabel("T (K)")
ax[1].set_ylabel("Total $\lambda$ (Wm$^{-2}$K$^{-1}$)")


# We can see that the system is indeed unstable between approximately 255K and 282 K, by finding out the lowest and highest entries where the ice albedo exceeds the Planck feedback, (although to be precise we should really solve with a root finding algorithm).

# In[9]:


#
# rough estimate (to nearest 2K) of unstable regime bounds
#
unstable_min=Ts[np.argwhere(feedback>0).min()]
unstable_max=Ts[np.argwhere(feedback>0).max()]
print("Earth system unstable between ",unstable_min," and ",unstable_max," K")


# ## Lesson take home messages:
# 
# <ul>
#   <li>An overall positive feedback indicates an unstable system</li>
#   <li>With only the Planck feedback operating, there is only one stable equilibrium possible in the system </li>
#   <li>With a nonlinear ice albedo feedback there is the possibility of hysteresis.  For some values of solar constant there are 3 equilibrium points, 2 stable and one unstable.  As S0 reduces the Earth can flip from the upper stable branch to the lower stable branch.  As S0 increases again the switch will occur at a difference (higher) solar constant value.</li>
#   <li>Thus ice-covered states are "sticky", once the Earth is frozen, difficult to get out</li>
#   <li> We note that when climate models are discussed, we often talk about them having a positive or negative feedback, this usually refers to all the feedbacks **except** the Planck feedback.  <u>Thus a climate model with positive feedback can still be stable as long as the magnitude does not exceed that of the Planck feedback.</u> </li>  
# </ul> 

# In[ ]:




