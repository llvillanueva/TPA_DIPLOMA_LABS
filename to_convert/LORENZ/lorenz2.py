from numpy import array, linspace
from scipy.integrate import solve_ivp
import pylab
from mpl_toolkits import mplot3d

def func(t, r):
    x, y, z = r 
    fx = 10 * (y - x)
    fy = 28 * x - y - x * z
    fz = x * y - (8.0 / 3.0) * z
    return array([fx, fy, fz], float)

# and I run it
r0 = [0, 1, 0]
sol = solve_ivp(func, [0, 50], r0, t_eval=linspace(0, 50, 5000))

# and plot it
fig = pylab.figure()
ax = pylab.axes(projection="3d")
ax.plot3D(sol.y[0,:], sol.y[1,:], sol.y[2,:], 'blue')
pylab.show()
