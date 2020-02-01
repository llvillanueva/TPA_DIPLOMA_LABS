#
# define grid and max number of iterations
#
res=1*pi/360.
lat=seq(0+res/2,pi/2,by=res)
sinlat=sin(lat)
coslat=cos(lat)
lat=lat*180/pi # convert to degrees for display
niter=200

#
# define my constants
#
s0=1370 # solar constant
olra=204.0
olrb=2.17
kt=3.81
icetemp=-10
icealb=0.3
freealb=0.3

#
# solar scale function
#
# sol=1-0.477*(0.5*(3.0*sinlat^2-1.0))  # parametrization 1
sol=0.25*(5-3*sinlat^2.0)               # parametrization 2 
sol=sol/mean(sol) # normalize

#
# initial conditions
#
temp=60*coslat-30
temp[]=-20

#
# set up plot 
#
pdf("zonal_model.pdf")
plot(lat,temp,type='l',col='blue',lwd=2,
  ylim=c(-60,30),xlim=c(0,90),
  ylab="Temperature (C)",xlab="Latitude (deg)")

#
# iteration loop
#
tbar0=-999 # memory of temperature for convergence

for (i in 1:niter){
   # albedo function
   alb=(temp<icetemp)*icealb+(temp>=icetemp)*freealb

   # averaged latitude
   tbar=sum(coslat*temp)/sum(coslat)
   tbar=mean(temp)

   # energy imbalance:
   delt=sol*s0/4*(1-alb)-(olra+olrb*temp) +kt*(tbar-temp)

   print(paste("solar",sol*s0/4*(1-alb)))
   print(paste("transport",kt*(tbar-temp)))
   print(paste("olr",olra+olrb*temp))

   # update T
   temp=temp+delt/20

   # uncomment if you want iteration lines:
   lines(lat,temp,type='l')

   # test for convergence:
   if (abs(tbar-tbar0)<1.e-9)
     {lines(lat,temp,type='l',col='red',lwd=3)
     lines(c(-90,90),c(icetemp,icetemp),lty=2,lwd=3,col='red')
     print(paste("convergence after interation ",i))
     print(paste("Final global mean temperature is ",tbar))
     q()}
   else
     {tbar0=tbar} 
}

print (paste("convergence failed after iterations: ",niter))

