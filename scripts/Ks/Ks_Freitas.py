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
bgs   = (['DJF'])
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

    yr1 = 1980
    yr2 = 1985
    fout = "../output/Ks/Ks_Fr1.{:s}{:d}-{:d}.erain.nc".format(ssn,yr1,yr2)

    i =0
    for yrs in range(yr1,yr2+1) :

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

    b=(1+np.cos(2*dtr*lats))/2;
    trm1=2*e_omega * b/radius
    c=um*b[:,None]
    d=1./b

    cdy = np.zeros_like(um)
    cdy2= np.zeros_like(um)
    cdy[:,:] = np.nan
    cdy2[:,:] = np.nan
    for i,count in enumerate(um[2:-2,:]) :
        # print i
        j = i+2
        cdy[j,:]=(np.mean(c[j-1:j+1,:])-np.mean(c[j:j+2,:]))/(np.mean(ym[j-1:j+1])-np.mean(ym[j:j+2]));
        # print dy[j,0]
        # print c[j-1:j+1,0]
        # print c[j:j+2,0]
        # print ym[j-1:j+1]
        # print ym[j:j+2],np.mean(ym[j:j+2])
    cdy=cdy*d[:,None]

    for i,count in enumerate(um[3:-3,:]) :
        # print i
        j = i+3
        cdy2[j,:]=(np.mean(cdy[j-1:j+1,:])-np.mean(cdy[j:j+2,:]))/(np.mean(ym[j-1:j+1])-np.mean(ym[j:j+2]));

    BetaM=trm1[:,None]-cdy2;

    print cdy2[10:40,6]
    print trm1[10:6]
    print BetaM[:,6]
    print trm1[100] - cdy2[100,6]
    print ym[10:40]
    print xm[6]
    print um[10:40,6]
    # exit()
    #
    #
    #
    #
    #
    # #--------------------------------------------------------------------------------
    #
    # cos2=coslat*coslat
    #
    # # for i in range(0,np.size(um,axis=1)) :
    # #  cosuy_np = np.gradient(um[:,i]*cos2,2*pi*radius/360)
    # #  cosuyy_np = np.gradient(um[:,i]/cos2,2*pi*radius/360)
    # cosuy_np = np.zeros_like(um)
    # cosuyy_np = np.zeros_like(um)
    # cosulat_np = np.zeros_like(um)
    # cosulat2_np = np.zeros_like(um)
    # #print dy
    # #quit()
    # dy = np.gradient(ym)
    #
    # # cosuy_np[:,0] = np.gradient(um[:,0]*cos2,dy)
    # # cosuyy_np[:,0] = np.gradient(cosuy_np[:,0]/cos2,dy)
    # for i,count in enumerate(um[0,:]) :
    #     # cosulat_np[:,i] = np.gradient(um[:,i]*cos2,dlat)
    #     cosuy_np[:,i] = np.gradient(um[:,i]*cos2,dy)
    #     cosuyy_np[:,i] = np.gradient(cosulat_np[:,i]/cos2,dy)
    #
    # dlat = np.gradient(lats)
    # for i,count in enumerate(um[0,:]) :
    #     # cosulat_np[:,i] = np.gradient(um[:,i]*cos2,dlat)
    #     cosulat_np[:,i] = np.gradient(um[:,i],dlat)
    #     cosulat2_np[:,i] = np.gradient(cosulat_np[:,i],dlat)
    #
    #
    # # tmp = 2*e_omega *cos2/radius
    # tmp =  2*e_omega *cos2 / radius
    # # BetaM_np = tmp[:,None]-cosulat2_np
    # BetaM_np = tmp[:,None]-cosuyy_np
    # BetaM = BetaM_np
    # #BetaM = BetaM_np + cosulat2_np*radius*radius

    Ks =  radius* np.sqrt(BetaM/um)

    # #  umasked = np.ma.masked_less(u, 5)
    #
    Ksm1 = np.ma.masked_where(um<5, Ks)
    Ksm2 = np.ma.masked_where(BetaM<0, Ksm1)
    print Ksm1[:,100]
    print Ksm2[:,100]

    Ks = Ksm2.filled()

    # print Ks.count()
    # print np.ma.count_masked(Ks)
    # print Ks.size
    # Ksf = Ks.filled()
    print Ks[:,150]

    print("Ks done")

    #---NetCDF write---------------------------------------------------------------
    print("Start NetCDF writing")

    ncout = Dataset(fout, 'w', format='NETCDF4')
    ncout.description = "Ks from %s" % (fin)

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
    ncout.close()


nc.close()
