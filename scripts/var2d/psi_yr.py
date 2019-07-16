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

for iyr in range(1979,2018):

    fin = "/Users/irudeva/work/DATA/ERAint/Plev/erain.hgt_air_wind.monmean.%d.nc"%iyr
    fout = "../data/erain.sf300.monmean.%d.nc"%iyr


    dimnam=('longitude','latitude','level','time')
    varnam=['longitude','latitude','level','time','z','t','u','v']

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
    levs = nc.variables[varnam[2]][:]
    time = nc.variables[varnam[3]][:]
    uwnd = nc.variables[varnam[6]][:]
    vwnd = nc.variables[varnam[7]][:]


    #uwnd = nc.variables[varnam[3]][:]
    #vwnd = nc.variables[varnam[4]][:]
    #lons = nc.variables[dimnam[0]][:]
    #lats = nc.variables[dimnam[1]][:]
    #time = nc.variables[dimnam[2]][:]
    #ncv = Dataset('wnd.mnth.eraint.nc), 'r')
    #vwnd = ncv.variables['vwnd'][:]
    #ncv.close()

    # find 300 hPa
    ilev = np.where(levs==300)[0]

    print("Data uploaded for %d"%iyr)

    #print uwnd.shape
    #print uwnd[1,1,1]

    # The standard interface requires that latitude and longitude be the leading
    # dimensions of the input wind components, and that wind components must be
    # either 2D or 3D arrays. The data read in is 3D and has latitude and
    # longitude as the last dimensions. The bundled tools can make the process of
    # re-shaping the data a lot easier to manage.
    uwnd1, uwnd_info = prep_data(np.squeeze(uwnd[:,ilev,:,:]), 'tyx')
    vwnd1, vwnd_info = prep_data(np.squeeze(vwnd[:,ilev,:,:]), 'tyx')

    # It is also required that the latitude dimension is north-to-south. Again the
    # bundled tools make this easy.
    lats, uwnd1, vwnd1 = order_latdim(lats, uwnd1, vwnd1)

    # Create a VectorWind instance to handle the computation of streamfunction and
    # velocity potential.
    w = VectorWind(uwnd1, vwnd1)

    # Compute the streamfunction and velocity potential. Also use the bundled
    # tools to re-shape the outputs to the 4D shape of the wind components as they
    # were read off files.
    sf, vp = w.sfvp()
    sf = recover_data(sf, uwnd_info)
    vp = recover_data(vp, uwnd_info)

    #print sf.shape
    #print sf[1,1,1]
    print("Streamfunction done")

    #---NetCDF write---------------------------------------------------------------
    print("Start NetCDF writing")

    ncout = Dataset(fout, 'w', format='NETCDF4')
    ncout.description = "Streamfunction form %s" % (fout)

    # Using our previous dimension info, we can create the new time dimension
    # Even though we know the size, we are going to set the size to unknown
    dimnamout=('longitude','latitude','time')
    varnamout=['longitude','latitude','time','u','v']

    ncout.createDimension(dimnamout[0], lons.size)
    ncout.createDimension(dimnamout[1], lats.size)
    ncout.createDimension(dimnamout[2], None)

    for nv in range(0, 2) :
        ncout_var = ncout.createVariable(varnamout[nv], nc.variables[varnamout[nv]].dtype,dimnamout[nv])
        for ncattr in nc.variables[varnamout[nv]].ncattrs():
            ncout_var.setncattr(ncattr, nc.variables[varnamout[nv]].getncattr(ncattr))
    for nv in range(2, 3) :

        ncout_var = ncout.createVariable(varnamout[nv], nc.variables[varnamout[nv]].dtype,dimnamout[nv])
        for ncattr in nc.variables[varnamout[nv]].ncattrs():
            ncout_var.setncattr(ncattr, nc.variables[varnamout[nv]].getncattr(ncattr))
    #print(nc.variables['latitude'].ncattrs())
    for nv in range(3, 5) :
        ncout_var = ncout.createVariable(varnamout[nv], nc.variables[varnamout[nv]].dtype,dimnamout[::-1])
        for ncattr in nc.variables[varnamout[nv]].ncattrs():
            ncout_var.setncattr(ncattr, nc.variables[varnamout[nv]].getncattr(ncattr))

    ncout.variables[dimnamout[0]][:] = lons
    ncout.variables[dimnamout[1]][:] = lats
    ncout.variables[dimnamout[2]][:] = time
    ncout.variables[varnamout[3]][:] = np.squeeze(uwnd[:,ilev,:,:])
    ncout.variables[varnamout[4]][:] = np.squeeze(vwnd[:,ilev,:,:])

    ncout_sf = ncout.createVariable('sf', 'f',dimnamout[::-1])
    ncout_sf.long_name = 'streamfunction'
    sf_scale = 1.e+7
    sf_add   = 0.
    ncout_sf.scale_factor = sf_scale
    ncout_sf.add_offset   = sf_add
    ncout_sf.units        = 'm**2 s**-1'

    #!!!automatically takes scale and offset into account
    #!!! no need for: ncout_sf[:] = (sf-sf_add)/sf_scale
    ncout_sf[:] = sf


    nc.close()
    ncout.close()
