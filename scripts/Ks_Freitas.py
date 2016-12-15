#import cartopy.crs as ccrs
#from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
#from cartopy.util import add_cyclic_point
#import matplotlib as mpl
#mpl.rcParams['mathtext.default'] = 'regular'
#import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from windspharm.standard import VectorWind
from windspharm.tools import prep_data, recover_data, order_latdim
import datetime as datetime  # Python standard library datetime  module


pi = np.pi
dtr=pi/180
rtd=180/pi
radius=6.371e6 #radius of sphere having same volume as Earth (m)
e_omega=7.292e-5 #rotation rate of Earth (rad/s)

fin = "../data/wnd300.mnth.erain.nc"

dimnam=('longitude','latitude','time')
varnam=['longitude','latitude','time','u','v']

# Read zonal and meridional wind components from file using the netCDF4
# module. The components are defined on pressure levels and are in separate
# files.
nc = Dataset(fin, 'r')
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



#uwnd = nc.variables[varnam[3]][:]
#vwnd = nc.variables[varnam[4]][:]
#lons = nc.variables[dimnam[0]][:]
#lats = nc.variables[dimnam[1]][:]
#time = nc.variables[dimnam[2]][:]
#ncv = Dataset('wnd.mnth.eraint.nc), 'r')
#vwnd = ncv.variables['vwnd'][:]
#ncv.close()
print("Data uploaded")

#bgs   = (['DJF','JJA','Sep'])
bgs   = (['JJA','DJF'])
for ssn in bgs :
    if ssn == 'DJF':
        bgmon = np.array([12, 1, 2])
    elif ssn == 'JJA':
        bgmon = np.array([6,7,8])
    elif ssn == 'Sep':
        bgmon = np.array([9])



    #----time average---------------------------------------------------------------------
    print 'Take u,v average'
    #  Time
    dt_time = [datetime.date(1900, 1, 1) + datetime.timedelta(hours=int(t))\
               for t in time]

    nt=np.array([0 for i in range(time.size)])
    i =0
    for yrs in range(1980,2011) :
        fout = "../output/Ks/Ks_Fr1.{:s}{%d}.erain.nc".format(ssn,yrs)

        for yr in range(yrs,yrs+1) :
            for m in bgmon :
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


        #----Mercator---------------------------------------------------------------------
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



        #----BetaM---------------------------------------------------------------------
        print 'Calculate BetaM'

        cos2=coslat*coslat

        # for i in range(0,np.size(um,axis=1)) :
        #  cosuy_np = np.gradient(um[:,i]*cos2,2*pi*radius/360)
        #  cosuyy_np = np.gradient(um[:,i]/cos2,2*pi*radius/360)
        cosuy_np = np.zeros_like(um)
        cosuyy_np = np.zeros_like(um)
        cosulat_np = np.zeros_like(um)
        cosulat2_np = np.zeros_like(um)
        #print dy
        #quit()
        dy = np.gradient(ym)

        cosuy_np[:,0] = np.gradient(um[:,0]*cos2,dy)
        cosuyy_np[:,0] = np.gradient(cosuy_np[:,0]/cos2,dy)

        dlat = np.gradient(lats)
        for i,count in enumerate(um[0,:]) :
            cosulat_np[:,i] = np.gradient(um[:,i]*cos2,dlat)
            cosulat2_np[:,i] = np.gradient(cosulat_np[:,i],dlat)


        # tmp = 2*e_omega *cos2/radius
        tmp = radius * 2*e_omega *cos2
        BetaM_np = tmp[:,None]-cosulat2_np
        BetaM = BetaM_np
        #BetaM = BetaM_np + cosulat2_np*radius*radius

        Ks = np.sqrt(BetaM/um)



        print("Ks done")

        #---NetCDF write---------------------------------------------------------------
        print("Start NetCDF writing")

        ncout = Dataset(fout, 'w', format='NETCDF4')
        ncout.description = "Ks %s" % (fout)

        # Using our previous dimension info, we can create the new time dimension
        # Even though we know the size, we are going to set the size to unknown

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

        ncout_Ks = ncout.createVariable('Ks', 'f',(dimnam[1],dimnam[0]))
        ncout_Ks.long_name = 'Total stationary wavenumber'
        #Ks_scale = 1.e-7
        #Ks_add   = 0.
        #ncout_Ks.scale_factor = Ks_scale
        #ncout_Ks.add_offset   = Ks_add
        #ncout_sf.units        = 'm**2 s**-1'
        ncout_um = ncout.createVariable('um', 'f',(dimnam[1],dimnam[0]))
        ncout_um.long_name = 'meridional wind component'
        ncout_BetaM = ncout.createVariable('BetaM', 'f',(dimnam[1],dimnam[0]))
        ncout_BetaM.long_name = 'cos(lat) times poleward gradient of abs vorticity'
        #ncout_BetaM.scale_factor = 1.e-7
        #ncout_BetaM.add_offset   = 0.

        #!!!automatically takes scale and offset into account
        #!!! no need for: ncout_sf[:] = (sf-sf_add)/sf_scale
        ncout_Ks[:] = Ks
        ncout_um[:] = um
        ncout_BetaM[:] = BetaM

        print ncout


        nc.close()
        ncout.close()
