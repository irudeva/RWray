from netCDF4 import Dataset
import numpy as np
import datetime as datetime  # Python standard library datetime  module
from  windspharm.standard import VectorWind
from windspharm.tools import prep_data, recover_data, order_latdim
from scipy import interpolate


# import cartopy.crs as ccrs
# from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
# from cartopy.util import add_cyclic_point


print "Calculating 2d ray paths"

# Parameters

pi = np.pi
dtr=pi/180
rtd=180/pi
radius=6.371e6 #radius of sphere having same volume as Earth (m)
e_omega=7.292e-5 #rotation rate of Earth (rad/s)


day=24*60*60 #in seconds
mins = 60
Periods = np.array([float('inf'), 50, 20 ])*day

freq = 2*pi/Periods
nfreq=freq.size
dt = 60 * mins   #time increment
int_time=10*day   #integration time
Nsteps = int_time/dt

k_wavenumbers=np.array([1, 2, 3, 4, 5, 6]) #  initial k wave number:

lon0 = np.array([0])
lat0 = np.array([30])

#set to 1 to do complex ray tracing
complex_tracing=False

print '---Parameters---------'
print "Wave periods: ", Periods/day, " days"
print "Wave numbers: ", k_wavenumbers
print "Periods: "
print "integration time ", int_time/day, " days"
print "time step ", dt/mins, " min"
print "Nsteps = ",Nsteps
print "Starting points: lon ",lon0,"E lat ",lat0,"N"
if complex_tracing is True :
    print "Complex tracing is on"
elif complex_tracing is False :
    print "Complex tracing is off"
else :
    print 'complex_tracing=',complex_tracing
    print "CHECK: What to do in case of complex solutions?"
    quit()
print '--------------------'

# Read data

fw = "../data/wnd300.mnth.erain.nc"
fsf = "../data/sf300.mnth.erain.nc"

dimnam=('longitude','latitude','time')
varnam=['longitude','latitude','time','u','v']

print "Reading wind from", fw
nc = Dataset(fw, 'r')
v=0
for var in varnam:
    if nc.variables[varnam[v]].name != var:
        print "Variables don't agree", var, nc.variables[varnam[v]].name, v
        exit()
    v += 1

lons = nc.variables[varnam[0]][:]
lats = nc.variables[varnam[1]][:]
time = nc.variables[varnam[2]][:]
uwnd = nc.variables[varnam[3]][:]
vwnd = nc.variables[varnam[4]][:]

print "Reading streamfunction from", fsf
nc1 = Dataset(fsf, 'r')
v=0
for var in varnam[0:2]:
    if nc1.variables[varnam[v]].name != var:
        print "ERROR reading ",fsf,":"
        print "Variables don't agree", var,"!=",nc.variables[varnam[v]].name
        exit()
    v += 1
lons1 = nc1.variables[varnam[0]][:]
lats1 = nc1.variables[varnam[1]][:]
time1 = nc1.variables[varnam[2]][:]
sf = nc.variables[varnam[3]][:]
if (lons1!=lons).any():
    print "ERROR wind.lons != streamfunction.lons"
    exit()
if (lats1!=lats).any():
    print "ERROR wind.lats != streamfunction.lats"
    exit()
if (lons1!=lons).any():
    print "ERROR wind.time != streamfunction.time"
    exit()


if(lats[0]<lats[-1]):
    print "ERROR: make sure that lat dim is N -> S"
    exit()

#  Time
dt_time = [datetime.date(1900, 1, 1) + datetime.timedelta(hours=int(t))\
           for t in time]

nt=np.array([0 for i in range(time.size)])
i =0
for yr in range(1980,1983) :
    for m in [12, 1, 2] :
        yr1 = yr
        if m == 12:
            yr1 = yr-1
        for t in dt_time :
            if t == datetime.date(yr1,m,1):
                print 'selected time: ', t
                ind = dt_time.index(t)
                nt[i] = ind
                i += 1


u = np.average(uwnd[nt[nt>0],:,:],axis=0)
v = np.average(vwnd[nt[nt>0],:,:],axis=0)
psi = np.average(sf[nt[nt>0],:,:],axis=0)

# Convert to Mercator projection

xm=lons*radius*dtr
xm360=360*radius*dtr
ym=lats+1  #array declaration
#ym[1:-2] = lats[1:-2]
ym[1:-2]=radius*np.log((1+np.sin(dtr*lats[1:-2]))/np.cos(dtr*lats[1:-2]));
ym[0]=float('inf')
ym[-1]=ym[0]

# dy = np.gradient(ym)
# dx = np.gradient(xm)

coslat=np.cos(dtr*lats)
#coslat[0]=0   # a very small number is used instead
#coslat[-1]=0  # ----"""----

# velocity in the Mercator projection
um=u/coslat[:,None]
vm=v/coslat[:,None]
# um checked!!!

# Create a VectorWind instance to handle the computations.
uwnd, uwnd_info = prep_data(uwnd, 'tyx')
vwnd, vwnd_info = prep_data(vwnd, 'tyx')

w = VectorWind(uwnd, vwnd)
# Compute absolute vorticity
q = w.absolutevorticity()

qbar = np.average(q[:,:,nt[nt>0]],axis=2)
print "qbar(4,0)=",qbar[4,0]
#qbar checked!!!


print "------------------------------"
print "gradients"
print np.version.version
print "  "
print "----- wind gradients ---------"

umx, umy = w.gradient(u)
vmx, vmy = w.gradient(v)
## umx, umy Checked!!!

#  alternatively  -  dy does not work!!!
#umy1, dum = np.gradient(um, dx)  ##ERROR:replace dx with dy!!!
#dum, umx1 = np.gradient(um, dx)

print "  "
print "----- q gradients ---------"

#Trick to calculate gradient for a scalar * cos(lat)

qx, qy = w.gradient(qbar*coslat[:,None])
#print "qmx[4,0]=",qmx[4,0]
#print "qmy[4,0]=",qmy[4,0]
#print "  "

#   alternatively -  dy does not work
#qy, dum = np.gradient(qbar, dx)  ##ERROR:replace dx with dy!!!
#dum, qx = np.gradient(qbar, dx)
#print "qx[9,4]=",qx[9,4]
#print "qmx[9,4]=",qmx[9,4]
#print "  "
# w.gradient is different from np.gradient at high lat >60N/S

print "  "
print "----- q second derivatives ---------"

#Trick to calculate gradient for a scalar * cos(lat)

qxx, qxy = w.gradient(qx*coslat[:,None])
qyx, qyy = w.gradient(qy*coslat[:,None])
# print "qmyy[4,0]=",qmyy[4,0]
# print "qmxx[30,0]=",qmxx[30,0]
# print "qmxy[30,5]=",qmxy[30,5]
# print "   "

# print "diff[4,5]: ", (qbar[4,4]-qbar[4,6])/(xm[4]-xm[6])
# print "qx[4,5]=",qx[4,5]
# print "qmx[4,5]=",qmx[4,5]
# print "  "
#
# print "diff[30,5]: ", (qbar[30,4]-qbar[30,6])/(xm[4]-xm[6])
# print "diff qbar[30,5]: ", (qbar[30,4]-qbar[30,6])
# print "diff xm[30,5]: ", (xm[4]-xm[6])
# print "qx[30,5]=",qx[30,5]
# print "qmx[30,5]=",qmx[30,5]
# print "  "


#----BetaM---------------------------------------------------------------------
print 'Calculate BetaM'

# BetaM == qy - checked !!!!

#coslat[0]=0
#coslat[-1]=cos1[0]
cos2=coslat*coslat

dum, cosuy=w.gradient(u*cos2[:,None])
dum, cosuyy = w.gradient(cosuy/coslat[:,None])

tmp = 2*e_omega *cos2/radius
BetaM=tmp[:,None]-cosuyy

#quit()

#---NetCDF write---------------------------------------------------------------
print("Start NetCDF writing")

ncvar = 'qx'
ftest = '../output/test/test.%s.nc' % (ncvar)
ncout = Dataset(ftest, 'w', format='NETCDF4')
ncout.description = "TEST %s" % (ftest)

# Using our previous dimension info, we can create the new time dimension
# Even though we know the size, we are going to set the size to unknown

dimnam=('longitude','latitude','time')
varnam=['longitude','latitude','time',ncvar]

ncout.createDimension(dimnam[0], lons.size)
ncout.createDimension(dimnam[1], lats.size)
#ncout.createDimension(dimnam[2], None)

for nv in range(0, 2) :
    ncout_var = ncout.createVariable(varnam[nv], nc.variables[varnam[nv]].dtype,dimnam[nv])
    for ncattr in nc.variables[varnam[nv]].ncattrs():
        ncout_var.setncattr(ncattr, nc.variables[varnam[nv]].getncattr(ncattr))
#print(nc.variables['latitude'].ncattrs())

ncout.variables[dimnam[0]][:] = lons
ncout.variables[dimnam[1]][:] = lats
#ncout.variables[dimnam[2]][:] = time
#ncout.variables[dimnam[2]][:] = 1

ncout_var = ncout.createVariable(ncvar, 'f',dimnam[1::-1])
#ncout_var.long_name = 'streamfunction'
var_scale = 1.e-12
var_add   = 0.
ncout_var.scale_factor = var_scale
ncout_var.add_offset   = var_add
#!!!automatically takes scale and offset into account
#!!! no need for: ncout_sf[:] = (sf-sf_add)/sf_scale
#ncout_var.units        = 'm**2 s**-1'
ncout_var.units        = 'not specified'

#print qx.shape
#print ncout_var.shape
ncout_var[:] = qx


nc.close()
ncout.close()
##---End NetCDF write---------------------------------------------------------------
print "All derivatives done"

##---Interpolation-----------------------------------------------------------------
print "  "
print "Interpolation"


umint = interpolate.interp2d(xm, ym[1:-1], um[1:-1,:], kind='cubic')
vmint = interpolate.interp2d(xm, ym[1:-1], vm[1:-1,:], kind='cubic')

umxint = interpolate.interp2d(xm, ym[1:-1], umx[1:-1,:], kind='cubic')
umyint = interpolate.interp2d(xm, ym[1:-1], umy[1:-1,:], kind='cubic')
vmxint = interpolate.interp2d(xm, ym[1:-1], vmx[1:-1,:], kind='cubic')
vmyint = interpolate.interp2d(xm, ym[1:-1], vmy[1:-1,:], kind='cubic')

qxint = interpolate.interp2d(xm, ym[1:-1], qx[1:-1,:], kind='cubic')
qyint = interpolate.interp2d(xm, ym[1:-1], qy[1:-1,:], kind='cubic')

qxxint = interpolate.interp2d(xm, ym[1:-1], qxx[1:-1,:], kind='cubic')
qyyint = interpolate.interp2d(xm, ym[1:-1], qyy[1:-1,:], kind='cubic')

qxyint = interpolate.interp2d(xm, ym[1:-1], qxy[1:-1,:], kind='cubic')


##---Derivatives(eq.9 and 10 in Karoly 1983)---------------------------------------

def ug(spotk,spotl,um,qx,qy) :
    Ks2=spotl*spotl+spotk*spotk
    Ks4=Ks2*Ks2
    return um+((spotk*spotk-spotl*spotl)*qy-2*spotk*spotl*qx)/Ks4

def vg(spotk,spotl,vm,qx,qy) :
    Ks2=spotl*spotl+spotk*spotk
    Ks4=Ks2*Ks2
    return vm+((spotk*spotk-spotl*spotl)*qx+2*spotk*spotl*qy)/Ks4

def kt(spotk,spotl,umx,vmx,qxy,qxx) :
    Ks2=spotl*spotl+spotk*spotk
    return -spotk*umx-spotl*vmx+(qxy*spotk-qxx*spotl)/Ks2

def lt(spotk,spotl,umy,vmy,qxy,qyy) :
    Ks2=spotl*spotl+spotk*spotk
    return -spotk*umy-spotl*vmy+(qyy*spotk-qxy*spotl)/Ks2

def rk(x,y,k,l):
    xt=ug(k,l,umint(x,y),qxint(x,y),qyint(x,y))
    yt=vg(k,l,vmint(x,y),qxint(x,y),qyint(x,y))
    dkdt=kt(k,l,umxint(x,y),vmxint(x,y),qxyint(x,y),qxxint(x,y))
    dldt=lt(k,l,umyint(x,y),vmyint(x,y),qxyint(x,y),qyyint(x,y))
    return xt,yt,dkdt,dldt

# Runge-Kutta method
def rk4(f, x0, y0, x1, n):
    vx = [0] * (n + 1)
    vy = [0] * (n + 1)
    h = (x1 - x0) / float(n)
    vx[0] = x = x0
    vy[0] = y = y0
    for i in range(1, n + 1):
        k1 = h * f(x, y)
        k2 = h * f(x + 0.5 * h, y + 0.5 * k1)
        k3 = h * f(x + 0.5 * h, y + 0.5 * k2)
        k4 = h * f(x + h, y + k3)
        vx[i] = x = x0 + i * h
        vy[i] = y = y + (k1 + k2 + k2 + k3 + k3 + k4) / 6
    return vx, vy

# Mercator to lonlat

def y2lat(a):
    return rtd*(2.0*np.arctan(np.exp(a/radius))-pi/2.0)

def lat2y(a) :
#lat in deg!!!
    return radius * np.log(np.tan(pi/4.0 + a * dtr/2.0))


###==================================================================================
print "  "
print "Start ray tracing:"
print "  "
##----------------------------------------------------------------------------------
# Solving for the ray path for different forcing sites (initial locations of rays):

Nloc = lon0.size
for iloc in range(0,Nloc) :
    print " Location #", iloc

    i = np.argmin(np.absolute(lons-lon0[iloc]))
    j = np.argmin(np.absolute(lats-lat0[iloc]))


    print "  Initial location of rays: "
    print "   Lon0: %6.2f corresponds to %6.2f" % (lon0[iloc],lons[i])
    print "   Lat0: %6.2f corresponds to %6.2f" % (lat0[iloc],lats[j])



    ##  Estimating the initial Ks from the forcing site for a specific BetaM and UbarM

    for fr in freq :
        #    period=round((2*pi/fr)/day);
        print "  Ray tracing: period", 2*pi/(fr*day)
        for k in k_wavenumbers :
            print "  initial k = ", k
            spotk = k/(radius*coslat[j])
            #print "spotk=", spotk

            ##  Calculate the initial l wave number from the initial omega
            ##  and k by solving the polynomial equation based on the
            ##  dispersion relation (equation 8 in Karoly 1983):
            ##  hange the following to have a non zero frequency:
            coeff = np.zeros(4)
            coeff[0]=vm[j,i]
            coeff[1]=um[j,i]*spotk-fr;
            coeff[2]=vm[j,i]*spotk*spotk+qx[j,i]
            coeff[3]=um[j,i]*np.power(spotk,3)-qy[j,i]*spotk-fr*spotk*spotk

            lroot = np.roots(coeff)
            print "  initial l = ", lroot*radius*coslat[j]
            for R in range(0,3) :
                spotl=lroot[R]
                print "  Root # ", R, "  spotl = ", spotl
                spotl=lroot[R]
                spotk = k/(radius*coslat[j]) #refresh!!!
                Ks=np.sqrt(spotl*spotl+spotk*spotk)

                if complex_tracing is False :
                    if np.not_equal(np.imag(spotl),0) :
                        print "   *** found complex initial l, not tracing. \n"
                        print " Location #", iloc
                        print "  Ray tracing: period", 2*pi/(fr*day)
                        print "  initial k ", k
                        print "  Root # ", R
                        quit()

                ## Starting the loop with the above initial k,l, and Ks

                lonn = np.zeros(Nsteps)
                latn = np.zeros(Nsteps)
                xn = np.zeros(Nsteps)
                yn = np.zeros(Nsteps)
                kn = np.zeros(Nsteps)
                ln = np.zeros(Nsteps)

#                for t in range(0,Nsteps) :


                for t in range(0,24) :
                    if np.equal(np.remainder(t,40),0) :
                       print "    t = ", t
                    # print "    t = ", t

                    if t==0 :
                        x0=xn[0]=xm[i]
                        y0=yn[0]=ym[j]
                        k0=kn[0]=spotk
                        l0=ln[0]=spotl
                        lonn[0]=lon0[iloc]
                        latn[0]=lat0[iloc]
                    else :
                        x0=xn[t-1]
                        y0=yn[t-1]
                        k0=kn[t-1]
                        l0=ln[t-1]

                    # # Runge-Kutta method

                    # RK step 1
                    kx0, ky0, kk0, kl0 = rk(x0,y0,k0,l0)

                    # RK step 2
                    x1 = x0+0.5*kx0*dt
                    y1 = y0+0.5*ky0*dt
                    k1=k0+0.5*kk0*dt
                    l1=l0+0.5*kl0*dt

                    kx1, ky1, kk1, kl1 = rk(x1,y1,k1,l1)

                    # RK step 3
                    x2 = x0+0.5*kx1*dt
                    y2 = y0+0.5*ky1*dt
                    k2=k0+0.5*kk1*dt
                    l2=l0+0.5*kl1*dt

                    kx2, ky2, kk2, kl2 = rk(x2,y2,k2,l2)

                    # RK step 4
                    x3 = x0+kx2*dt
                    y3 = y0+ky2*dt
                    k3=k0+kk2*dt
                    l3=l0+kl2*dt

                    kx3, ky3, kk3, kl3= rk(x3,y3,k3,l3)

                    #RK4 results
                    dx=dt*(kx0+2*kx1+2*kx2+kx3)/6;
                    dy=dt*(ky0+2*ky1+2*ky2+ky3)/6;
                    dk=dt*(kk0+2*kk1+2*kk2+kk3)/6;
                    dl=dt*(kl0+2*kl1+2*kl2+kl3)/6;

                    # print ' '
                    # print 'x0=',x0, ' kx0=', kx0
                    # print 'x1=',x1, ' kx1=', kx1
                    # print 'x2=',x2, ' kx2=', kx2
                    # print 'x3=',x3, ' kx0=', kx0
                    # print 'dx=',dx/dt
                    #
                    # print ' '
                    # print 'y0=',y0, ' ky0=', ky0
                    # print 'y1=',y1, ' ky1=', ky1
                    # print 'y2=',y2, ' ky2=', ky2
                    # print 'y3=',y3, ' ky0=', ky0
                    # print 'dy=',dy/dt
                    #
                    # print ' '
                    # print 'k0=',k0, ' kk0=', kk0
                    # print 'k1=',k1, ' kk1=', kk1
                    # print 'k2=',k2, ' kk2=', kk2
                    # print 'k3=',k3, ' kk0=', kk0
                    # print 'dk=',dk/dt
                    #
                    # print ' '
                    # print 'l0=',l0, ' kl0=', kl0
                    # print 'l1=',l1, ' kl1=', kl1
                    # print 'l2=',l2, ' kl2=', kl2
                    # print 'l3=',l3, ' kl0=', kl0
                    # print 'dl=',dl/dt

                    tn=t+1
                    xn[tn] = x0+dx
                    if xn[tn]>=xm360 :
                     xn[tn]=xn[tn]-xm360;
                    yn[tn] = y0+dy
                    kn[tn] = k0+dk
                    ln[tn] = l0+dl

                    # print ' '
                    # print 'x0+dx=',x0,'+',dx,'=',xn
                    # print 'y0+dy=',y0,'+',dy,'=',yn
                    # print 'k0+dk=',k0,'+',dk,'=',kn
                    # print 'l0+dl=',l0,'+',dl,'=',ln

                    # Ks =np.sqrt(kn*kn+ln*ln)

                    ## Finding the location

                    lonn[tn] = xn[tn]*rtd/radius
                    latn[tn] = y2lat(yn[tn])

                    print lonn[tn], latn[tn]


                fout = open('myfile','w')
                for t in range(0,25) :
                    fout.write("%f6.2 %f6.2 %f6.2 %f6.2" % (xn[t],yn[t],kn[t]/(radius*np.cos(yn[t]*dtr)),ln[t]/(radius*np.cos(yn[t]*dtr))))
                fout.close()






                print "All good to here"
                quit()
