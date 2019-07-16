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

def leapyear(year):
    if year/400 == 0 :
        return False
    if year%100 == 0 :
        return False
    if year%4 == 0 :
        return True

XXwnd = 22  # critical wind

yr1 = 1983
yr2 = 2017

# clmon = np.array([3,9])
# chmon = np.array(['Mar','Sep'], dtype=str)
# clmon = np.array([3])
# chmon = np.array(['Mar'], dtype=str)
clmon = np.array([1])
chmon = np.array(['Jan'], dtype=str)


for iyr in range(yr1,yr2+1):
    print iyr
    fout = "../output/wind/frwnd6h_gt%d.%s%d.erain.nc"%(XXwnd,chmon[0],iyr)


    fin = "/Users/irudeva/work/DATA/ERAint/Plev/erain.hgt_air_wind.6h.%d.nc"%(iyr)
    # fout = "../data/sf300.mnth.erain.nc"
    print "read file %s"%(fin)


    dimnam=('longitude','latitude','level','time')
    varnam=['longitude','latitude','level','time','z','t','pv','u','v']

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
    lev  = nc.variables[varnam[2]][:]
    # if iyr == yr1:
    #     # tmp = nc.variables[varnam[3]][:]
    #     time = np.zeros((nyr,1464))
    #     uwnd = np.zeros((nyr,1464,lev.size,lats.size,lons.size))
    #     vwnd = np.zeros((nyr,1464,lev.size,lats.size,lons.size))

    time = nc.variables[varnam[3]][:]
    uwnd = nc.variables[varnam[7]][:]
    vwnd = nc.variables[varnam[8]][:]
    # del tmp

    nc.close()

# print("Data uploaded")

# quit()
#----selecting timesteps-------------------------------------------------------------
    print "creating time mask..."
    btime = np.zeros_like(time,dtype = bool)
    bwnd = np.zeros_like(uwnd,dtype = bool)

# for iyr in range(yr1,yr2+1):
#     print iyr

    # tsize = 1460
    # if leapyear(iyr):
    #     tsize =1464
    # print tsize
    # dt_time = [datetime.datetime(1900, 1, 1) + datetime.timedelta(hours=int(t))\
    #        for t in time[iyr-yr1,:tsize]]
    dt_time = [datetime.datetime(1900, 1, 1) + datetime.timedelta(hours=int(t))\
           for t in time]


    # print dt_time
    # quit()
    # datetime.date(iyr,m,1)

    for m in clmon :
        # yr1 = iyr
        for t in dt_time :
            if t.month == m :
                print 'selected time: ', t
                ind = dt_time.index(t)
                # print ind
                btime[ind] = 1
                bwnd[ind,:,:,:] = 1

    del dt_time


    # mx = ma.masked_array(x, mask=[0, 0, 0, 1, 0])
    # mask = (m1 == 1) & (m2 == 1)
    #
    # mask_wnd = wnd > 25
    #

    #-------masking wind ---------------------------------------------------
    print "calculating wnd..."
    uwnd2 = np.square(uwnd)
    print "uwnd2 done"
    vwnd2 = np.square(vwnd)
    print "vwnd2 done"

    wnd = np.sqrt(uwnd2 + vwnd2)
    wnd = uwnd
    print "wnd done"

    print "masking wind..."
    mwnd = np.ma.array(wnd, mask = ~bwnd)
    mwndXX = np.ma.array(mwnd, mask = wnd<XXwnd)
    mwndXX_vpos = np.ma.array(mwndXX, mask = vwnd<0)
    mwndXX_vneg = np.ma.array(mwndXX, mask = vwnd>0)

    print "counting wnd points..."
    frwndXX = mwndXX.count(axis=0)
    frwndXX_vpos = mwndXX_vpos.count(axis=0)
    frwndXX_vneg = mwndXX_vneg.count(axis=0)
    # print numwnd22.shape

    # u = np.average(uwnd[nt[nt>0],:,:],axis=0)
    # v = np.average(vwnd[nt[nt>0],:,:],axis=0)
    #------compressing----------------------------------------------------------------

    # msize  = np.count_nonzero(btime)
    # cltime = np.zeros(msize)
    # clyear = np.zeros(msize)
    # clmon  = np.zeros(msize)
    # cldate = np.zeros(msize)
    # clhr   = np.zeros(msize)
    # cluwnd = np.zeros((msize,lev.size,lats.size,lons.size))
    # clvwnd = np.zeros((msize,lev.size,lats.size,lons.size))
    #
    # i = 0
    #
    # for iyr in range(yr1,yr2+1):
    #     for bt in btime[iyr - yr1,:]:
    #         if bt :
    #             cltime[i] = time[iyr - yr1,t]
    #             cdate = datetime.datetime(1900, 1, 1) + datetime.timedelta(hours=int(time[iyr - yr1,t]))
    #             clyr[i] = cdate.year
    #             clmon[i] = cdate.month
    #             cldate[i] = cdate.date
    #             clhr[i] = cdate.hour
    #             cluwnd[i,:,:,:] = uwnd[iyr,t,:,:,:]
    #             clvwnd[i,:,:,:] = vwnd[iyr,t,:,:,:]
    #             i += 1




    #-------------------------------------------------------------------------------

    # # The standard interface requires that latitude and longitude be the leading
    # # dimensions of the input wind components, and that wind components must be
    # # either 2D or 3D arrays. The data read in is 3D and has latitude and
    # # longitude as the last dimensions. The bundled tools can make the process of
    # # re-shaping the data a lot easier to manage.
    # uwnd, uwnd_info = prep_data(uwnd, 'tyx')
    # vwnd, vwnd_info = prep_data(vwnd, 'tyx')
    # #print uwnd.shape
    #
    # # It is also required that the latitude dimension is north-to-south. Again the
    # # bundled tools make this easy.
    # lats, uwnd, vwnd = order_latdim(lats, uwnd, vwnd)
    #
    # # Create a VectorWind instance to handle the computation of streamfunction and
    # # velocity potential.
    # w = VectorWind(uwnd, vwnd)
    #
    # # Compute the streamfunction and velocity potential. Also use the bundled
    # # tools to re-shape the outputs to the 4D shape of the wind components as they
    # # were read off files.
    # sf, vp = w.sfvp()
    # sf = recover_data(sf, uwnd_info)
    # vp = recover_data(vp, uwnd_info)
    #
    # #print sf.shape
    # #print sf[1,1,1]
    # print("Streamfunction done")

    #---NetCDF write---------------------------------------------------------------
    print("Start NetCDF writing")

    ncout = Dataset(fout, 'w', format='NETCDF4')
    ncout.description = "jet streams form %s" % (fin)

    nc = Dataset(fin, 'r')  # open last fin to copy attributes

    # Using our previous dimension info, we can create the new time dimension
    # Even though we know the size, we are going to set the size to unknown

    ncout.createDimension(dimnam[0], lons.size)
    ncout.createDimension(dimnam[1], lats.size)
    ncout.createDimension(dimnam[2], lev.size)
    ncout.createDimension(dimnam[3], None)

    for nv in range(0, 3) :
        ncout_var = ncout.createVariable(varnam[nv], nc.variables[varnam[nv]].dtype,dimnam[nv])
        for ncattr in nc.variables[varnam[nv]].ncattrs():
            ncout_var.setncattr(ncattr, nc.variables[varnam[nv]].getncattr(ncattr))
    #print(nc.variables['latitude'].ncattrs())

    ncout.variables[dimnam[0]][:] = lons
    ncout.variables[dimnam[1]][:] = lats
    ncout.variables[dimnam[2]][:] = lev
    #  for time:
    nv = 3  #
    ncout_var = ncout.createVariable(varnam[nv], 'i2',dimnam[nv])
    ncout_var.long_name = 'years'
    ncout_var.months = chmon[0]
    ncout_var[:] = np.arange(iyr,iyr+1)

    ncout_var1 = ncout.createVariable('frwnd', 'i2',dimnam[2::-1])
    ncout_var1.long_name = 'frequency_wind_greater_%d'%(XXwnd)
    # sf_scale = 1.e+7
    # sf_add   = 0.
    # ncout_sf.scale_factor = sf_scale
    # ncout_sf.add_offset   = sf_add
    ncout_var1.units        = 'm s**-1'

    print dimnam[2::-1]
    print frwndXX.shape

    #!!!automatically takes scale and offset into account
    #!!! no need for: ncout_sf[:] = (sf-sf_add)/sf_scale
    ncout_var1[:] = frwndXX

    ncout_var1 = ncout.createVariable('frwnd_vpos', 'i2',dimnam[2::-1])
    ncout_var1.long_name = 'frequency_wind_greater_%d_v_positive'%(XXwnd)
    ncout_var1.units        = 'm s**-1'
    ncout_var1[:] = frwndXX_vpos

    ncout_var1 = ncout.createVariable('frwnd_vneg', 'i2',dimnam[2::-1])
    ncout_var1.long_name = 'frequency_wind_greater_%d_v_negative'%(XXwnd)
    ncout_var1.units        = 'm s**-1'
    ncout_var1[:] = frwndXX_vneg


    nc.close()
    ncout.close()

    del uwnd
    del vwnd
    del uwnd2
    del vwnd2
    del wnd

    del mwnd
    del mwndXX
    del mwndXX_vpos
    del mwndXX_vneg

    del frwndXX
    del frwndXX_vpos
    del frwndXX_vneg

    del btime
    del bwnd

    del time
    del lons
    del lats
    del lev
