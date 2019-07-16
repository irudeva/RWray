"""Compute Rossby wave source from the long-term mean flow.

This example uses the standard interface.

Additional requirements for this example:

* netCDF4 (http://unidata.github.io/netcdf4-python/)
* matplotlib (http://matplotlib.org/)
* cartopy (http://scitools.org.uk/cartopy/)

"""
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.util import add_cyclic_point
import matplotlib as mpl
mpl.rcParams['mathtext.default'] = 'regular'
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import datetime as datetime  # Python standard library datetime  module

from windspharm.standard import VectorWind
from windspharm.tools import prep_data, recover_data, order_latdim
from windspharm.examples import example_data_path

gravacc =  9.80665

# The first date of calculation
fyear = 1979
fmon = 1

# The last date of calculation
lyear = 2017
lmon = 12

# direcotry name
diri = "/Users/irudeva/work/DATA/ERAint/Plev/"

# bgs   = 'DJF'
# if bgs == 'DJF':
#     bgmon = np.array([12, 1, 2])
# elif bgs == 'JJA':
#     bgmon = np.array([6,7,8])


lev = "200"  #hgt level for wind
lev1 = "300"  # hgt level

dimnam=('longitude','latitude','level','time')
varnam=['longitude','latitude','level','time',"z","t",'u','v']

for iyr in range(1979,lyear+1) :
    fw = "{}erain.hgt_air_wind.monmean.{:d}.nc".format(diri,iyr)
    print "Reading wind from", fw
    nc = Dataset(fw, 'r')
    v=0
    for var in varnam:
        if nc.variables[varnam[v]].name != var:
            print "Variables don't agree", var, nc.variables[varnam[v]].name, v
            exit()
        v += 1

    if (iyr == fyear) :
        lons = nc.variables[varnam[0]][:]
        lats = nc.variables[varnam[1]][:]
        levs = nc.variables[varnam[2]][:]
        time = nc.variables[varnam[3]][:]
        z    = nc.variables[varnam[4]][:]
        uwnd = nc.variables[varnam[6]][:]
        vwnd = nc.variables[varnam[7]][:]
    else:
        if (nc.variables[varnam[0]][:]!=lons).any():
            print "ERROR lons != old lons in {:d}".format(iyr)
            exit()
        if (nc.variables[varnam[1]][:]!=lats).any():
            print "ERROR lats != old lats in {:d}".format(iyr)
            exit()
        if (nc.variables[varnam[2]][:]!=levs).any():
            print "ERROR levs != old levs in {:d}".format(iyr)
            exit()
        time = np.append(time,nc.variables[varnam[3]][:],0)
        z    = np.append(   z,nc.variables[varnam[4]][:],0)
        uwnd = np.append(uwnd,nc.variables[varnam[6]][:],0)
        vwnd = np.append(vwnd,nc.variables[varnam[7]][:],0)

z = z/gravacc

if(lats[0]<lats[-1]):
    print "ERROR: make sure that lat dim is N -> S"
    exit()


# print "Reading hgt from", fhgt
# nc1 = Dataset(fhgt, 'r')
# v=0
# for var in varnam[0:2]:
#     if nc1.variables[varnam[v]].name != var:
#         print "ERROR reading ",fhgt,":"
#         print "Variables don't agree", var,"!=",nc.variables[varnam[v]].name
#         exit()
#     v += 1
# lons1 = nc1.variables[varnam[0]][:]
# lats1 = nc1.variables[varnam[1]][:]
# time1 = nc1.variables[varnam[2]][:]
# z = nc1.variables['z'][:]
#
# if (lons1!=lons).any():
#     print "ERROR wind.lons != hgt.lons"
#     exit()
# if (lats1!=lats).any():
#     print "ERROR wind.lats != hgt.lats"
#     exit()
# if (lons1!=lons).any():
#     print "ERROR wind.time != hgt.time"
#     exit()

#  Time
dt_time = [datetime.date(1900, 1, 1) + datetime.timedelta(hours=int(t))\
           for t in time]

# climatological averaging
# nt=np.array([0 for i in range(time.size)])
#
# for iy in range(1980,1981) :
#     nt = nt*0
#     i =0
#     # for yr in range(iy,iy+1) :
#     for yr in range(iy,iy+31) :
#         print yr
#         for m in bgmon :
#             yr1 = yr
#             if m == 12:
#                 yr1 = yr-1
#             for t in dt_time :
#                 if t == datetime.date(yr1,m,1):
#                     print 'selected time: ', t
#                     ind = dt_time.index(t)
#                     nt[i] = ind
#                     i += 1

    # u = np.average(uwnd[nt[nt>0],:,:],axis=0)
    # v = np.average(vwnd[nt[nt>0],:,:],axis=0)
    # hgt = np.average(z[nt[nt>0],:,:],axis=0)/gravacc

# no climatogical average

for il, ilev in enumerate(levs):
    if ilev==300 :
        u = uwnd[:,il,:,:]
        v = vwnd[:,il,:,:]
        hgt = z[:,il,:,:]

            # print u.shape
            # print uwnd.shape
            #
            #
            # u, u_info = prep_data(u, 'tyx')
            # v, v_info = prep_data(v, 'tyx')
            #
            # print u_info
            #
            # quit()
            #
            #
            # # It is also required that the latitude dimension is north-to-south. Again the
            # # bundled tools make this easy.
            # latsO, uwnd, vwnd = order_latdim(lats, uwnd, vwnd)
        #
        # Create a VectorWind instance to handle the computations.
        w = VectorWind(u, v)


        # Compute components of rossby wave source: absolute vorticity, divergence,
        # irrotational (divergent) wind components, gradients of absolute vorticity.
        eta = w.absolutevorticity()
        div = w.divergence()
        uchi, vchi = w.irrotationalcomponent()
        etax, etay = w.gradient(eta)

        # Combine the components to form the Rossby wave source term. Re-shape the
        # Rossby wave source array to the 4D shape of the wind components as they were
        # read off files.
        S = -eta * div - (uchi * etax + vchi * etay)
        # S = recover_data(S, uwnd_info)



            # Pick out the field for December and add a cyclic point (the cyclic point is
            # for plotting purposes).
            # S_plt, lons_c = add_cyclic_point(S, lons)


            # Plot Rossby wave source.
            # ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
            # clevs = [-30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30]
            # # fill = ax.contourf(lons_c, lats, S_plt * 1e11, clevs, cmap=plt.cm.RdBu_r,
            # #                    transform=ccrs.PlateCarree(), extend='both')
            # # ax.coastlines()
            # # ax.gridlines()
            # print ax
            # ax.set_xticks([0, 60, 120, 180, 240, 300, 359.99], crs=ccrs.PlateCarree())
            # ax.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
            # lon_formatter = LongitudeFormatter(zero_direction_label=True,
            #                                    number_format='.0f')
            # lat_formatter = LatitudeFormatter()
            # ax.xaxis.set_major_formatter(lon_formatter)
            # ax.yaxis.set_major_formatter(lat_formatter)
            # # plt.colorbar(fill, orientation='horizontal')
            # plt.title('Rossby Wave Source ($10^{-11}$s$^{-1}$)', fontsize=16)
            # print 'start to plot'
            # plt.show()

        ncvar = ['u', 'v', 'RWS', 'hgt']
        # print 'ncvar=',ncvar
        ftest = '../output/RWS/RWS.hgt%s.mon%d_%d.nc' % (ilev,fyear,lyear)
        ncout = Dataset(ftest, 'w', format='NETCDF4')
        ncout.description = "lev = %s hPa" % (lev)

        # Using our previous dimension info, we can create the new time dimension
        # Even though we know the size, we are going to set the size to unknown

        dimnam=('longitude','latitude','time')
        varnam=['longitude','latitude','time',ncvar]

        ncout.createDimension(dimnam[0], lons.size)
        ncout.createDimension(dimnam[1], lats.size)
        ncout.createDimension(dimnam[2], None)

        for nv in range(0, 3) :
            ncout_var = ncout.createVariable(varnam[nv], nc.variables[varnam[nv]].dtype,dimnam[nv])
            for ncattr in nc.variables[varnam[nv]].ncattrs():
                ncout_var.setncattr(ncattr, nc.variables[varnam[nv]].getncattr(ncattr))
        #print(nc.variables['latitude'].ncattrs())

        ncout.variables[dimnam[0]][:] = lons
        ncout.variables[dimnam[1]][:] = lats
        ncout.variables[dimnam[2]][:] = time
        #ncout.variables[dimnam[2]][:] = 1

        varform =  "({:d},{:d},{:d})f4".format(time.size,lats.size,lons.size)
        varlist = np.zeros(16, dtype = {'names': ['name', 'outname', 'data', 'scale', 'units'],
                                        'formats': ['a5', 'a5', varform, 'f4', 'a13']} )

        varlist[0] = ("u","u",u,1,"m s**-1")
        varlist[1] = ("v","v",v,1,"m s**-1")
        varlist[2] = ("S","RWS",S, 1e-11,"10**-11 s**-1")
        varlist[3] = ("hgt","hgt",hgt,1,"m")
        varlist[4] = ("eta","eta",eta, 1e-5,"10**-5 s**-1")
        varlist[5] = ("div","div",div, 1e-6,"10**-6 s**-1")


        for iv in range(6) :

            var = varlist["outname"][iv]
            print var,dimnam[2::-1]
            ncout_var = ncout.createVariable(var, 'f',dimnam[2::-1])
            #ncout_var.long_name = 'streamfunction'
            # var_scale = varlist["scale"][iv]
            # var_add   = 0.
            ncout_var.scale_factor = varlist["scale"][iv]
            ncout_var.add_offset   = 0.
            #!!!automatically takes scale and offset into account
            #!!! no need for: ncout_sf[:] = (sf-sf_add)/sf_scale
            #ncout_var.units        = 'not specified'
            # ncout_var.units        = 'scale   %s' % varlist["scale"][iv]
            ncout_var.units        = varlist["units"][iv]

            #print qx.shape
            # print "shape",ncout_var.shape
            ncout_var[:] = varlist["data"][iv]


        ncout.close()


        nc.close()
