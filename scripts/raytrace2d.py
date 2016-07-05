from netCDF4 import Dataset
import numpy as np
import datetime as datetime  # Python standard library datetime  module

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

#  Time
dt_time = [datetime.date(1900, 1, 1) + datetime.timedelta(hours=int(t))\
           for t in time]

nt=np.array([0 for i in range(time.size)])
i =0
for yr in range(1980,1982) :
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

print psi.shape


# Convert to Mercator projection

xm=lons*radius*dtr
ym=lats
ym[1:-2]=radius*np.log((1+np.sin(dtr*lats[1:-2]))/np.cos(dtr*lats[1:-2]));
print ym[2]
ym[0]=float('inf')
ym[-1]=ym[0]
