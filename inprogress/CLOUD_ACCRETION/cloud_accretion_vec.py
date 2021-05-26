import numpy as np
#import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
from mpl_toolkits.mplot3d import Axes3D

def terminal_w(r):
    """calculation of droplet terminal fallspeed"""
    w=np.where(r<30.e-6,X1*r**2,X2*r)
    w=np.where(r>1.e-3,X3*np.sqrt(r),w)
    return(w)

def r2mass(r):
    """converts radius to mass"""
    return((4.0*rho_l*np.pi*r**3)/3.0)

def mass2r(m):
    """converts radius to mass"""
    r=((3.0*m)/(4.0*rho_l*np.pi))**(1./3.)
    return(r)

class cloud():
    """class of a cloud with arrays of x,y,z,r and w"""
    def __init__(self,xpos,ypos,zpos,radius):

        # position vectors and radius (could use array)
        #self.x=ma.masked_array(xpos,mask=np.zeros(len(xpos)))
        #self.y=ma.masked_array(ypos,mask=np.zeros(len(ypos)))
        #self.z=ma.masked_array(zpos,mask=np.zeros(len(zpos)))
        #self.r=ma.masked_array(radius,mask=np.zeros(len(radius)))

        self.x=xpos
        self.y=ypos
        self.z=zpos
        self.r=radius
        self.dead=np.zeros(len(radius))
        
        # vertical velocity
        self.w=terminal_w(radius) # velocity  
        
    def collision(self):
        from scipy.spatial.distance import cdist

        ndrop=len(self.r) # how many drops?
        
        # need to get pairwise sum of radii to determine collisions sphere
        rsq=np.tile(self.r,(ndrop,1))
        rsum=rsq+np.transpose(rsq) # paired radius sum

        # now find collisions
        # check if pairs in same "x-y" space
        coords=np.column_stack((self.x,self.y))    
        xydist=cdist(coords,coords)
        
        # and if points have "swapped" their z-locations, i.e.
        # one cloud drop has overtaken another?
        zd1=np.array([i-j for i in self.z
                          for j in self.z]).reshape(ndrop,ndrop)
        zd2=np.array([i-j for i in self.znew
                          for j in self.znew]).reshape(ndrop,ndrop)
        
        # if zd1*zd2<0 then drop overtake has occurred
        ind1,ind2=np.where((xydist<rsum)&(zd1*zd2<0.0))
        
        # list of unique pairs neglecting diagonals (impacts itself!)
        unique=(ind1<ind2)
        ind1 = ind1[unique]
        ind2 = ind2[unique]

        # loops over lower triangle of all column-collisions
        for i1,i2 in zip(ind1,ind2):
            if (self.dead[i1]):
                continue
            if (self.dead[i2]):
                continue
            print ("check ",xydist[i1,i2],rsum[i1,i2],zd1[i1,i2],zd2[i2,i2],self.z[i1],self.z[i2],self.znew[i1],self.znew[i2],self.w[i1],self.w[i2])
            
            # new mass, radius and terminal velocity 
            mass1=r2mass(self.r[i1])
            mass2=r2mass(self.r[i2])
            print("mass ",mass1," bumped with ",mass2)
            
            self.r[i1]=mass2r(mass1+mass2)
            self.w[i1]=terminal_w(self.r[i1])
            self.znew[i1]=min(self.znew[i1],self.znew[i2])

            print("new radius,w ",self.r[i1]*1.e6,self.w[i1])
            # other drop is now dead:
            self.dead[i2]=True
            self.r[i2]=0.0
            
            #self.x.mask[i2]=self.y.mask[i2]=True
            #self.z.mask[i2]=True
            #self.r.mask[i2]=True
        

def main():
    # dimensions of the domain
    global xsize,ysize,zsize
    
    global X1,X2,X3
    global dt # timestep
    global rho_l
    global rsmall,rlarge

    # physical constants:
    #-------------------
    # terminal velocity coeffs
    X1,X2,X3=1.2e8,8e3,250. 
    # density of liquid water
    rho_l=1000.
    slice=1

    # domain size in meters
    xsize=ysize=100e-6
    zsize=2000

    # vertical extend of cloud
    zcloud1=1300
    zcloud2=2000

    #----------------------------------------------    
    # parameters that can be varied for experiment:
    #----------------------------------------------
    # timestep (seconds) 
    dt=50.0

    #
    # Student exercise 1: change radii and liquid water amount:
    #
    liq=0.5e-3 # cloud liq water content in kg/m**3    
    rsmall=5.e-6 # microns
    rlarge=20.e-6
    ratio_large=0.05 # proportion of large drops   
    #nlarge=200
    lanimate=True
    
    # drop density per m**3
    
    ndropdens=liq*3.0/(rho_l*4.0*np.pi*rsmall**3)
    
    # L=N*rho_l*4/3 pi r^3
    # ndrop=10 # total number of drops
    # initial proportion of large drops with rlarge radius
    # not used if distrbution of drop sizes assumed
    ndrop=int(ndropdens*xsize*ysize*(zcloud2-zcloud1))
    if (ndrop>7000):
        print("ndrop too large ",ndrop)
        exit()
    print("total number of drops",ndrop)

    nlarge=int(np.ceil(ratio_large*ndrop))
    print("total number of large drops",nlarge)
    
    # Student exercise 2: change distribution: 
    # can insert code here to simply set a distribution of radii:
    # set arrange from lognormal distribution, 
    # dropr=np.random.DIST(moments)

    # initial random positions:
    # could use a dictionary, but don't want to lose numpy advantage?
    dropx=np.random.uniform(low=0,high=xsize,size=ndrop)
    dropy=np.random.uniform(low=0,high=ysize,size=ndrop)
    dropz=np.random.uniform(low=zcloud1,high=zcloud2,size=ndrop)
    
    print("initialize")
    
    # one large drop falling through a cloud of small drops
    dropr=np.full(ndrop,fill_value=rsmall)
    dropr[-nlarge:]=rlarge

    # initial drop conditions
    drops=cloud(dropx,dropy,dropz,dropr)

    print("run")    
    # set up plot window
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.set_xlim(0,xsize)
    ax.set_ylim(0,ysize)
    ax.set_zlim(0,zsize)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_zlabel("z (m)")


    if (lanimate):
        sc3d=ax.scatter(drops.x[::slice],drops.y[::slice],drops.z[::slice],c='blue',marker='o')
        title = ax.set_title('Slab Cloud')

        def animate(i):
            timestep(drops)
            sc3d._offsets3d=(drops.x[::slice],drops.y[::slice],drops.z[::slice])  # update the data
            title.set_text('Slab Cloud, time={0:4.2f} seconds'.format(i*dt))
            return(sc3d)
    
        ani = animation.FuncAnimation(fig, animate, 
                                  interval=50, 
                                  blit=False)

        plt.show()
    else:
        for itime in range(1000):
            print(itime)
            timestep(drops)
            
def timestep(drops):
    drops.znew=drops.z-dt*drops.w
    drops.collision()
    drops.z=drops.znew # update
    #print(drops.z[1])
    
    #detect_rain()

    # drop falls out of bottom, make dead and add to rain stats
    #print(drops.z)


main()
    

    
