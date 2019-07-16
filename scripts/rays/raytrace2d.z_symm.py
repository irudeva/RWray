from netCDF4 import Dataset
import numpy as np
import datetime as datetime  # Python standard library datetime  module
from  windspharm.standard import VectorWind
from windspharm.tools import prep_data, recover_data, order_latdim
from scipy import interpolate



# This is a variation of the RW ray tracing for zonally symmetric flow, i.e.
# u=u(y)  v=0
# ! take zonal average of u

# Then:
# omega = um*k-BetaM*k/K**2
# ug = w/k+2*BetaM*k**2/K**4
# vg = 2*BetaM*k*l/K**4
# dk/dt = 0
# dl/dt = 0


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


day2s=24*60*60 #in seconds
min2s = 60
#Periods = np.array([float('inf'), -50, -14 ])*day2s
#Periods = np.array([float('inf')])*day2s
Periods = np.array([-14])*day2s
#Periods = np.array([float('inf'), 14, -14])*day2s

freq = 2*pi/Periods
nfreq=freq.size
dt = 15 * min2s   #time increment in s
int_time=15*day2s  #integration time
Nsteps = int_time/dt

#k_wavenumbers=np.array([1, 2, 3, 4, 5, 6]) #  initial k wave number:
k_wavenumbers=np.array([4]) #  initial k wave number:

lon0 = np.array([180])
lat0 = np.array([20])
loc = np.array([0])  # the range of locations used

#set to 1 to do complex ray tracing
complex_tracing=False

print '---Parameters---------'
print "Wave periods: ", Periods/day2s, " day(s)"
print "Wave numbers: ", k_wavenumbers
print "Periods: "
print "integration time ", int_time/day2s, " day(s)"
print "time step ", dt/min2s, " min"
print "Nsteps = ",Nsteps
#print "Starting points: lon ",lon0[loc[0]:loc[1]+1],"E lat ",lat0[loc[0]:loc[1]+1],"N"
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

# Set a zonally symmetric flow!

uzon = np.mean(u,axis=1)
uwndzon = np.mean(uwnd,axis=2)
for lat in range(0,np.size(u,axis=0)) :
    for tt in range(0,np.size(uwnd,axis=0)) :
        uwnd[tt,:,:] = uwndzon[tt,lat]
    u[lat,:] = uzon[lat]
v = np.zeros_like(v)
vwnd = np.zeros_like(vwnd)



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

for j in range(0,np.size(um,axis=0)) :
    print lats[j], ym[j],u[j,i],um[j,i]


# Create a VectorWind instance to handle the computations.
uwnd, uwnd_info = prep_data(uwnd, 'tyx')
vwnd, vwnd_info = prep_data(vwnd, 'tyx')


#w = VectorWind(uwnd, vwnd)
# Compute absolute vorticity
w = VectorWind(u, v)
q = w.absolutevorticity()

qbar = q
#qbar = np.average(q[:,:,nt[nt>0]],axis=2)
print "qbar(4,0)=",qbar[4,0]
#qbar checked!!!


print "------------------------------"
print "gradients"
print np.version.version
print "  "
print "----- wind and q gradients ---------"

# umx, umy = w.gradient(u)
# vmx, vmy = w.gradient(v)

dx = np.gradient(xm)
umx = np.zeros_like(um)
vmx = np.zeros_like(vm)
qx = np.zeros_like(q)
qxx = np.zeros_like(q)
qxy = np.zeros_like(q)
for j in range(0,np.size(um,axis=0)) :
    umx[j,:] = np.gradient(um[j,:],dx)
    vmx[j,:] = np.gradient(vm[j,:],dx)
    qx[j,:] = np.gradient(q[j,:],dx)
    qxx[j,:] = np.gradient(qx[j,:],dx)

dy = np.gradient(ym)
umy = np.zeros_like(um)
vmy = np.zeros_like(vm)
qy = np.zeros_like(q)
qyy = np.zeros_like(q)
for i in range(0,np.size(um,axis=1)) :
    umy[:,i] = np.gradient(um[:,i],dy)
    vmy[:,i] = np.gradient(vm[:,i],dy)
    qy[:,i] = np.gradient(q[:,i],dy)
    qyy[:,i] = np.gradient(qy[:,i],dy)
    qxy[j,:] = np.gradient(qy[j,:],dx)

# print "  "
# print "----- q gradients ---------"
#
#Trick to calculate gradient for a scalar * cos(lat)

#qx, qy = w.gradient(qbar*coslat[:,None])

# print "  "
# print "----- q second derivatives ---------"

#Trick to calculate gradient for a scalar * cos(lat)

# qxx, qxy = w.gradient(qx*coslat[:,None])
# qyx, qyy = w.gradient(qy*coslat[:,None])

#----BetaM---------------------------------------------------------------------
print 'Calculate BetaM'

cos2=coslat*coslat

# dum, cosuy=w.gradient(u*cos2[:,None])
# dum, cosuyy = w.gradient(cosuy/coslat[:,None])
#
# tmp = 2*e_omega *cos2/radius
# BetaM=tmp[:,None]-cosuyy

#alternative BetaM

# for i in range(0,np.size(um,axis=1)) :
#  cosuy_np = np.gradient(um[:,i]*cos2,2*pi*radius/360)
#  cosuyy_np = np.gradient(um[:,i]/cos2,2*pi*radius/360)
cosuy_np = np.zeros_like(um)
cosuyy_np = np.zeros_like(um)
#print dy
#quit()
cosuy_np[:,0] = np.gradient(um[:,0]*cos2,dy)
cosuyy_np[:,0] = np.gradient(cosuy_np[:,0]/cos2,dy)
for j in range(0,np.size(um,axis=0)) :
    cosuy_np[j,:] = cosuy_np[j,0]
    cosuyy_np[j,:] = cosuyy_np[j,0]

tmp = 2*e_omega *cos2/radius
BetaM_np = tmp[:,None]-cosuyy_np
BetaM = BetaM_np


# for j in range(0,np.size(um,axis=0)) :
#     for i in range(0,np.size(um,axis=1)) :
#        print "          ",cosuy[j,i],cosuy_np[j,i]
#        print i,j,lons[i]," ",lats[j]," ",BetaM[j,i],BetaM_np[j,i]," ", BetaM[j,i]-BetaM_np[j,i]
#
#----Ks------------------------------------------------------------------------
# print 'Calculate Ks'
#
# Ks=np.zeros_like(u)
#
# for j in range(0,np.size(um,axis=0)) :
#     for i in range(0,np.size(um,axis=1)) :
#         Ks[j,i] = np.sqrt(BetaM[j,i]/um[j,i])
#
# Ksy = np.zeros_like(Ks)
# for i in range(0,np.size(Ks,axis=1)) :
#     Ksy[:,i] = np.gradient(Ks[:,i],dy)


#---NetCDF write---------------------------------------------------------------
print("Start NetCDF writing")

varlist = np.zeros(15, dtype = {'names': ['name', 'outname', 'data', 'scale'],
                                'formats': ['a5', 'a5', '(241,480)f4', 'f4']} )


varlist[0] = ("u","u",u,1)
varlist[1] = ("v","v",v,1)
varlist[2] = ("um","um",um,1)
varlist[3] = ("vm","vm",vm,1)
varlist[4] = ("umx","umx",umx,1.e-6)
varlist[5] = ("umy","umy",umy,1.e-6)
varlist[6] = ("vmx","vmx",vmx,1.e-6)
varlist[7] = ("vmy","vmy",vmy,1.e-6)
varlist[8] = ("qbar","q",qbar,1.e-4)
varlist[9] = ("qx","qx",qx,1.e-12)
varlist[10] = ("qy","qy",qy,1.e-11)
varlist[11] = ("qxx","qxx",qxx,1.e-18)
varlist[12] = ("qyy","qyy",qyy,1.e-18)
varlist[13] = ("qxy","qxy",qxy,1.e-18)
varlist[14] = ("BetaM","BetaM",BetaM,1.e-11)


#Arr= [("u",100),("v",200)]
#for Key,Value in Arr:
#    print Key,"=",Value

print 'u=',u[10,10]
print 'v=',v[10,10]
print 'um=',um[10,10]
print 'vm=',vm[10,10]
print 'umx=',umx[10,10]
print 'umy=',umy[10,10]
print 'vmx=',vmx[10,10]
print 'vmy=',vmy[10,10]
print 'q=',qbar[10,10]
print 'qx=',qx[10,10]
print 'qy=',qy[10,10]
print 'qxx=',qxx[10,10]
print 'qyy=',qyy[10,10]
print 'qxy=',qxy[10,10]
print 'Beta=',BetaM[10,10]


for iv in range(varlist['name'].size) :


    ncvar = varlist["outname"][iv]
    print 'ncvar=',ncvar
    ftest = '../output/zon_symm/test.%s.nc' % (varlist["outname"][iv])
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
    # var_scale = varlist["scale"][iv]
    # var_add   = 0.
    ncout_var.scale_factor = varlist["scale"][iv]
    ncout_var.add_offset   = 0.
    #!!!automatically takes scale and offset into account
    #!!! no need for: ncout_sf[:] = (sf-sf_add)/sf_scale
    #ncout_var.units        = 'm**2 s**-1'
    #ncout_var.units        = 'not specified'
    ncout_var.units        = 'scale   %s' % varlist["scale"][iv]

    #print qx.shape
    #print ncout_var.shape
    ncout_var[:] = varlist["data"][iv]


    ncout.close()
nc.close()

##---End NetCDF write---------------------------------------------------------------
print "All derivatives done"


##---Interpolation-----------------------------------------------------------------
print "  "
print "Interpolation"


uint = interpolate.interp2d(xm, ym[1:-1], u[1:-1,:], kind='cubic')
vint = interpolate.interp2d(xm, ym[1:-1], v[1:-1,:], kind='cubic')

umint = interpolate.interp2d(xm, ym[1:-1], um[1:-1,:], kind='cubic')
vmint = interpolate.interp2d(xm, ym[1:-1], vm[1:-1,:], kind='cubic')

umxint = interpolate.interp2d(xm, ym[1:-1], umx[1:-1,:], kind='cubic')
umyint = interpolate.interp2d(xm, ym[1:-1], umy[1:-1,:], kind='cubic')
vmxint = interpolate.interp2d(xm, ym[1:-1], vmx[1:-1,:], kind='cubic')
vmyint = interpolate.interp2d(xm, ym[1:-1], vmy[1:-1,:], kind='cubic')

qint = interpolate.interp2d(xm, ym[1:-1], qbar[1:-1,:], kind='cubic')

qxint = interpolate.interp2d(xm, ym[1:-1], qx[1:-1,:], kind='cubic')
qyint = interpolate.interp2d(xm, ym[1:-1], qy[1:-1,:], kind='cubic')

qxxint = interpolate.interp2d(xm, ym[1:-1], qxx[1:-1,:], kind='cubic')
qyyint = interpolate.interp2d(xm, ym[1:-1], qyy[1:-1,:], kind='cubic')

qxyint = interpolate.interp2d(xm, ym[1:-1], qxy[1:-1,:], kind='cubic')

BetaMint = interpolate.interp2d(xm, ym[1:-1], BetaM[1:-1,:], kind='cubic')

##---Derivatives(eq.9 and 10 in Karoly 1983)---------------------------------------

def Kt(k,um,fr,BetaM):
    Kt = np.nan
    if BetaM > 0 and um-fr/k > 0:
        Kt = np.sqrt(BetaM/(um-fr/k))
    return Kt

def ug(k,l,fr,BetaM) :
    Ks2=l*l+k*k
    Ks4=Ks2*Ks2
    return fr/k+2*BetaM*k*k/Ks4

def vg(k,l,BetaM) :
    Ks2=l*l+k*k
    Ks4=Ks2*Ks2
    return 2*BetaM*k*l/Ks4

# def kt(k,l,umx,vmx,qxy,qxx) :
    # return 0

#def lt(k,l,umy,vmy,qxy,qyy) :
def lt(k,fr,um,vg,l,BetaM,ydel) :
     print 'lt:'
     print ' um', um
     print ' BetaM', BetaM
     print ' ydel',ydel
     K_lt = np.empty(2)
     K_lt[:] = np.nan
     Ky_lt = np.zeros_like(K_lt)
     if all(BetaM >= 0) and all(um-fr/k > 0) :
         for i in range(2) :
             K_lt[i] = Kt(k,um[i],fr,BetaM[i])
         print K_lt
         Ky_lt = (K_lt[1]-K_lt[0])/ydel
         print ' K=',K_lt
         print ' K*=',K_lt*radius
         print ' Ky=', Ky_lt
         lt = vg*K_lt[0]*Ky_lt/l
         print ' lt=',lt
         print ' '
     else :
         print 'ERROR: Either BetaM or um < 0 '
         print ' BetaM=', BetaM
         print ' um=', um
         lt = -1
     return lt

# Runge-Kutta method
def rk(x,y,k,l,fr):
    xt=ug(k,l,fr,BetaMint(x,y))
    yt=vg(k,l,BetaMint(x,y))
    dkdt=0
    #dldt:
    y1=y
    y2=y+yt*dt
    print ' x=',x,' y=',y1, ' u(x,y)=',uint(x,y1), ' um(x,y)=',umint(x,y1)
    print ' x=',x,' y=',y2, ' u(x,y)=',uint(x,y2), ' um(x,y)=',umint(x,y1)
    yrange = np.linspace(y1,y2,num=2)
    dldt=lt(k,fr,np.array([umint(x,y1),umint(x,y2)]),yt,l,np.array([BetaMint(x,y1),BetaMint(x,y2)]),y2-y1)
    return xt,yt,dkdt,dldt


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

#Nloc = lon0.size
for iloc in loc :
    print " Location #", iloc

    i = np.argmin(np.absolute(lons-lon0[iloc]))
    j = np.argmin(np.absolute(lats-lat0[iloc]))


    print "  Initial location of rays: "
    print "   Lon0: %6.2f corresponds to %6.2f" % (lon0[iloc],lons[i])
    print "   Lat0: %6.2f corresponds to %6.2f" % (lat0[iloc],lats[j])



    ##  Estimating the initial Ks from the forcing site for a specific BetaM and UbarM

    for fr in freq :
        #    period=round((2*pi/fr)/day);
        print "  Ray tracing: period", 2*pi/(fr*day2s)
        for k in k_wavenumbers :
            print "  initial k = ", k
            spotk = k/(radius*coslat[j])
            #print "spotk=", spotk

            ##  Calculate the initial l wave number from the initial omega
            ##  and k by solving the polynomial equation based on the
            ##  dispersion relation (equation 8 in Karoly 1983):
            ##  hange the following to have a non zero frequency:
            # coeff = np.zeros(2)
            # coeff[0]=um[j,i]*spotk-fr;
            # coeff[1]=0
            # coeff[2]=um[j,i]*np.power(spotk,3)+fr*spotk*spotk-BetaM[j,i]*spotk
            #
            # lroot = np.roots(coeff)
            # print "  initial l = ", lroot*radius

            linit = np.square(Kt(spotk,um[j,i],fr,BetaM[j,i]))-spotk*spotk
            if linit >= 0 :
                linit = np.sqrt(linit)
            else:
                linit = np.nan
            lroot = np.array([linit,-linit])

            print 'linit ', linit
            print np.isnan(linit)

            if not np.isnan(linit) :

                print "  initial l = ", lroot*radius

                for R in range(0,2) :
                    spotl=lroot[R]
                    print "  Root # ", R, "  spotl = ", spotl
                    #spotk = k/(radius*coslat[j]) #refresh!!!
                    #Ks=np.sqrt(spotl*spotl+spotk*spotk)

                    # testomega=um[j,i]*spotk+vm[j,i]*spotl+(qx[j,i]*spotl-qy[j,i]*spotk)/(Ks*Ks)
                    # print 'um=',um[j,i]
                    # print 'vm=',vm[j,i]
                    #
                    # print 'spotk=',spotk
                    # print 'spotl=',spotl
                    # print 'Ks2=',Ks*Ks
                    #
                    # print 'qx=',qx[j,i]
                    # print 'qy=',qy[j,i]
                    #
                    #
                    # print 'testomega=', testomega
                    #

                    if complex_tracing is False :
                        if np.not_equal(np.imag(spotl),0) :
                            print "   *** found complex initial l, not tracing. "
                            print "   *** Location #", iloc
                            print "   *** Ray tracing: period", 2*pi/(fr*day)
                            print "   *** initial k ", k
                            print "   *** Root # \n", R
                            continue

                    ## Starting the loop with the above initial k,l, and Ks

                    lonn = np.empty(Nsteps+1)
                    latn = np.empty(Nsteps+1)
                    xn = np.empty(Nsteps+1)
                    yn = np.empty(Nsteps+1)
                    kn = np.empty(Nsteps+1)
                    ln = np.empty(Nsteps+1)

                    lonn[:] = np.nan
                    latn[:] = np.nan
                    xn[:] = np.nan
                    yn[:] = np.nan
                    kn[:] = np.nan
                    ln[:] = np.nan

                    for t in range(0,Nsteps) :
                    #for t in range(0,24) :
                        print '  t=',t
                        # if np.equal(np.remainder(t,40),0) :
                        #    print "    t = ", t

                        if t==0 :
                            x0=xn[0]=xm[i]
                            y0=yn[0]=ym[j]
                            k0=kn[0]=spotk
                            l0=ln[0]=np.real(spotl)
                            lonn[0]=lon0[iloc]
                            latn[0]=lat0[iloc]
                        else :
                            x0=xn[t]
                            y0=yn[t]
                            k0=kn[t]
                            l0 = np.square(Kt(k0,umint(x0,y0),fr,BetaMint(x0,y0)))-k0*k0
                            if l0 >= 0:
                                l0 = np.sqrt(l0)
                                if ln[t] >= 0 :
                                    l0 = l0
                                else :
                                    l0 = -l0
                            else:
                                l0 = np.nan
                            # if R == 0 :
                            print ' l0=',l0,'ln=',ln[t]
    # zs_vln
                            l0=ln[t]

                            # print x0,y0
                            print 'umint=',umint(x0,y0),'BetaM=',BetaMint(x0,y0)
                            # print Kt(k0,umint(x0,y0),fr,BetaMint(x0,y0))
                            # print k0,l0
                            # print ' '

                        if np.isnan(l0) :
                            break

                        # # Runge-Kutta method

                        # RK step 1
                        print 'RK step 1'
                        # print 'x0=', x0
                        # print 'y0=', y0
                        kx0, ky0, kk0, kl0 = rk(x0,y0,k0,l0,fr)
                        print ' l0, kl0', l0, kl0
                        if kl0 == -1 :
                            print ' RAY Terminated: dldt = nan'
                            break

                        # RK step 2
                        print 'RK step 2'
                        x1 = x0+0.5*kx0*dt
                        y1 = y0+0.5*ky0*dt
                        k1=k0+0.5*kk0*dt
                        # print 'k1:'
                        # print ' k0=', k0
                        # # print ' kk0=',kk0
                        # # print ' dt=',dt
                        # # print ' kk0*dt',kk0*dt
                        # print ' 0.5*kk0*dt',0.5*kk0*dt
                        # print ' k1', k1
                        # print ' k1*',k1*(radius*coslat[j])
                        # print 'endk1'
                        l1=l0+0.5*kl0*dt
                        # print ' l0, l1, lkl0', l0, l1, kl0
                        # print ' l1*',l1*radius*coslat[j]

                        kx1, ky1, kk1, kl1 = rk(x1,y1,k1,l1,fr)
                        if kl1 == -1 :
                            print ' RAY Terminated: dldt = nan'
                            break

                        # RK step 3
                        print 'RK step 3'
                        x2 = x0+0.5*kx1*dt
                        y2 = y0+0.5*ky1*dt
                        k2=k0+0.5*kk1*dt
                        # print 'k2:'
                        # print ' k0=', k0
                        # print ' 0.5*kk1*dt',0.5*kk1*dt
                        # print ' k2', k2
                        # print ' k2*',k2*(radius*coslat[j])
                        # print 'endk2'
                        l2=l0+0.5*kl1*dt
                        # print ' l1, l2, kl1', l1, l2, kl1
                        # print ' l2*',l2*radius

                        kx2, ky2, kk2, kl2 = rk(x2,y2,k2,l2,fr)
                        if kl2 == -1 :
                            print ' RAY Terminated: dldt = nan'
                            break

                        # RK step 4
                        print 'RK step 4'
                        x3 = x0+kx2*dt
                        y3 = y0+ky2*dt
                        k3=k0+kk2*dt
                        # print 'k3:'
                        # print ' k0=', k0
                        # print ' 0.5*kk2*dt',0.5*kk2*dt
                        # print ' k3', k3
                        # print ' k3*',k3*(radius*coslat[j])
                        # print 'endk3'
                        l3=l0+kl2*dt
                        # print ' l2, l3, kl2', l2, l3, kl2
                        # print ' l3*',l3*(radius)

                        kx3, ky3, kk3, kl3= rk(x3,y3,k3,l3,fr)
                        if kl3 == -1 :
                            print ' RAY Terminated: dldt = nan'
                            break

                        #RK4 results
                        dx=dt*(kx0+2*kx1+2*kx2+kx3)/6;
                        dy=dt*(ky0+2*ky1+2*ky2+ky3)/6;
                        dk=dt*(kk0+2*kk1+2*kk2+kk3)/6;
                        dl=dt*(kl0+2*kl1+2*kl2+kl3)/6;

                        tn=t+1
                        xn[tn] = x0+dx
                        if xn[tn]>=xm360 :
                         xn[tn]=xn[tn]-xm360;
                        yn[tn] = y0+dy

                        kn[tn] = k0+dk
                        ln[tn] = l0+dl

                        # print ' '
                        # print ' ug',dx/dt
                        # print ' vg',dy/dt
                        # if np.absolute(dy/dt)<0.5 or np.isnan(dl):
                        if np.isnan(dl):
                            print 'Ray terminated: vg =', dy/dt
                            break

                        # print 'x0+dx=',x0,'+',dx,'=',xn[tn]
                        # print 'y0+dy=',y0,'+',dy,'=',yn[tn]
                        # print 'k0+dk=',k0,'+',dk,'=',kn[tn]
                        print ' kl0,kl1,kl2,kl3=',kl0,kl1,kl2,kl3
                        print ' dl/dt', dl/dt
                        print 'l0+dl=',l0,'+',dl,'=',ln[tn]

                        # Ks =np.sqrt(kn*kn+ln*ln)

                        ## Finding the location

                        lonn[tn] = xn[tn]*rtd/radius
                        latn[tn] = y2lat(yn[tn])


                    if fr==0 :
                        fout = open('../output/zon_symm/raypath_zs_vln_loc{:d}N_{:d}E_period{}_k{:d}_root{:d}'.format(lat0[iloc],lon0[iloc],'_inf',k,R),'w')
                    else :
                        fout = open('../output/zon_symm/raypath_zs_vln_loc{:d}N_{:d}E_period{:0.0f}_k{:d}_root{:d}'.format(lat0[iloc],lon0[iloc],2*pi/(fr*day2s),k,R),'w')
                    frmt = "{:>5} {:>3} {:>4}"+" {:>6}"*2+(" {:>6}"+" {:>9}")*3+" {:>9}"*3+" {:>7}"*4+" {:>9}"*11+" \n"
                    fout.write(frmt.
                        format('t','hr','day','lon','lat','k*rad','k','l*rad','l','l0*rad','l0','K','KK','Kom','u','v','um','vm','umx','vmx','umy','vmy','q','qx','qy','BetaM','qxx','qyy','qxy'))
                    frmt = "{:>5d} {:>3d} {:>4.1f}"+" {:>6.2f}"*2+(" {:>6.2f}"+" {:>9.2e}")*3+" {:>9.2f}"+" {:>9.2e}"*2+" {:>7.2f}"*4+" {:>9.2e}"*11+" \n"
                    for t in range(0,Nsteps+1,12*3600/dt) :
                        x=xn[t]
                        y=yn[t]

                        KK=np.sqrt(kn[t]*kn[t]*np.square(np.cos(latn[t]*dtr))+ln[t]*ln[t])*radius
                        KK1 = np.sqrt(kn[t]*kn[t] + ln[t]*ln[t])
                        KKom = Kt(kn[t],umint(x,y),fr,BetaMint(x,y))
                        if np.isnan(KKom) :
                            KKom=np.array([np.nan])

                        l0 = np.square(Kt(kn[t],umint(x,y),fr,BetaMint(x,y)))-kn[t]*kn[t]
                        if l0 >= 0 :
                            l0 = np.sqrt(l0)
                            if ln[t] >= 0 :
                                l0 = l0
                            else :
                                l0 = -l0
                        else:
                            l0=np.array([np.nan])
                        #
                        # print t,t*dt/3600,t*dt/(3600*24.),lonn[t],latn[t]
                        #
                        # print kn[t]*radius*np.cos(latn[t]*dtr)
                        # print kn[t]
                        # print ln[t]*radius
                        # print ln[t]
                        # print l0[0]*radius
                        # print l0[0]
                        # print KK,KK1,KKom[0]
                        # print uint(x,y)[0],vint(x,y)[0],umint(x,y)[0],vmint(x,y)[0],umxint(x,y)[0],vmxint(x,y)[0],umyint(x,y)[0],vmyint(x,y)[0]
                        # print qint(x,y)[0],qxint(x,y)[0],qyint(x,y)[0],BetaMint(x,y)[0],qxxint(x,y)[0],qyyint(x,y)[0],qxyint(x,y)[0]

                        fout.write(frmt.format(t,t*dt/3600,t*dt/(3600*24.),lonn[t],latn[t],
                        kn[t]*radius*np.cos(latn[t]*dtr),kn[t],ln[t]*radius,ln[t],l0[0]*radius,l0[0],
                        KK,KK1,KKom[0],
                        uint(x,y)[0],vint(x,y)[0],umint(x,y)[0],vmint(x,y)[0],umxint(x,y)[0],vmxint(x,y)[0],umyint(x,y)[0],vmyint(x,y)[0],
                        qint(x,y)[0],qxint(x,y)[0],qyint(x,y)[0],BetaMint(x,y)[0],qxxint(x,y)[0],qyyint(x,y)[0],qxyint(x,y)[0]))
                    fout.close()
                    print 'loc{:d}N_{:d}E_period{}_k{:d}_root{:d}'.format(lat0[iloc],lon0[iloc],'_inf',k,R)


                print "All good to here"
                #quit()
