from netCDF4 import Dataset
import numpy as np
import datetime as datetime  # Python standard library datetime  module
from  windspharm.standard import VectorWind
from windspharm.tools import prep_data, recover_data, order_latdim


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
    print "CHECK: What to do in case of caplex solutions?"
    exit()
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
ym=lats+1  #array declaration
#ym[1:-2] = lats[1:-2]
ym[1:-2]=radius*np.log((1+np.sin(dtr*lats[1:-2]))/np.cos(dtr*lats[1:-2]));
ym[0]=float('inf')
ym[-1]=ym[0]

dy = np.gradient(ym)
dx = np.gradient(xm)

um = u+1
vm = u+1
i = 0
print "before u=",u[1,1] #!!!!!!!!!!make sure u have not changed!!!!
for lat in lats:
    um[i,:]=u[i,:]/np.cos(lat*dtr)
    vm[i,:]=v[i,:]/np.cos(lat*dtr)
    i += 1
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

#  alternatively
umy1, dum = np.gradient(um, dx)  ##ERROR:replace dx with dy!!!
dum, umx1 = np.gradient(um, dx)
#umy1, umx1 = np.gradient(um, dx,dy)
# print "um[5,0]=",um[5,0]
# print "um[5,2]=",um[5,2]
# print "dx=",dx[1]
# print "umx1[5,1]=",umx1[5,1]
# print "umx1[100,100]=",umx1[100,100]


print "  "
print "----- q gradients ---------"

#Trick to calculate gradient for a scalar * cos(lat)
qbarnm = qbar+1
i=0
for lat in lats:
    qbarnm[i,:]=qbar[i,:]*np.cos(lat*dtr)
    i += 1

qmx, qmy = w.gradient(qbarnm)
print "qmx[4,0]=",qmx[4,0]
print "qmy[4,0]=",qmy[4,0]
print "  "

#   alternatively
qy, dum = np.gradient(qbar, dx)  ##ERROR:replace dx with dy!!!
dum, qx = np.gradient(qbar, dx)
print "qx[9,4]=",qx[9,4]
print "qmx[9,4]=",qmx[9,4]
print "  "

print "  "
print "----- q second derivatives ---------"

#Trick to calculate gradient for a scalar * cos(lat)
qmxm = qbar+1
qmym = qbar+1
i=0
for lat in lats:
    qmxm[i,:]=qmx[i,:]*np.cos(lat*dtr)
    qmym[i,:]=qmy[i,:]*np.cos(lat*dtr)
    i += 1

qmxx, qmxy = w.gradient(qmxm)
qmyx, qmyy = w.gradient(qmym)
print "qmyy[4,0]=",qmyy[4,0]
print "qmxx[30,0]=",qmxx[30,0]
print "qmxy[30,5]=",qmxy[30,5]
print "   "

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


cos2=lats+1
cos2[1:-2]=np.cos(dtr*lats[1:-2])
cos2[1:-2]=cos2[1:-2]*cos2[1:-2]
cos2[0]=float('inf')
cos2[-1]=cos2[0]

cosux, cosuy=w.gradient(u*cos2)
cosuyy = qbar+1
i=0
for lat in lats:
    cosuyy[i,:]=cosuy[i,:]*np.cos(lat*dtr)
    i += 1
dum, cosuyy = w.gradient(cosuy)

tmp = cos2+1
tmp[1:-2] = 2*e_omega *cos2[1:-2]/radius
tmp[0]=float('inf')
tmp[-1]=tmp[0]
BetaM = np.ones(241,480)
BetaM[1:-2,:] = tmp[-1:2]*ones

BetaM[1:-2]=2*e_omega *cos2[1:-2]/radius-cosuyy

BetaM[0]=float('inf')
BetaM[-1]=BetaM[0]



print BetaM[4]



#---NetCDF write---------------------------------------------------------------
print("Start NetCDF writing")

ncvar = 'qy'
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
var_scale = 1.e-18
var_add   = 0.
ncout_var.scale_factor = var_scale
ncout_var.add_offset   = var_add
#!!!automatically takes scale and offset into account
#!!! no need for: ncout_sf[:] = (sf-sf_add)/sf_scale
#ncout_var.units        = 'm**2 s**-1'
ncout_var.units        = 'not specified'

#print qx.shape
#print ncout_var.shape
ncout_var[:] = qy


nc.close()
ncout.close()
#---End NetCDF write---------------------------------------------------------------
