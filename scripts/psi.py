#import cartopy.crs as ccrs
#from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
#from cartopy.util import add_cyclic_point
#import matplotlib as mpl
#mpl.rcParams['mathtext.default'] = 'regular'
#import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np

# test me

from windspharm.standard import VectorWind
from windspharm.tools import prep_data, recover_data, order_latdim

fin = "wnd300.mnth.erain.nc"
fout = "sf300.mnth.erain.nc"


dimnam=('longitude','latitude','time')

# Read zonal and meridional wind components from file using the netCDF4
# module. The components are defined on pressure levels and are in separate
# files.
nc = Dataset(fin, 'r')
uwnd = nc.variables['u'][:]
vwnd = nc.variables['v'][:]
lons = nc.variables[dimnam[0]][:]
lats = nc.variables[dimnam[1]][:]
nc.close()
#ncv = Dataset('wnd.mnth.eraint.nc), 'r')
#vwnd = ncv.variables['vwnd'][:]
#ncv.close()
print("Data uploaded")

#print uwnd.shape
#print uwnd[1,1,1]

# The standard interface requires that latitude and longitude be the leading
# dimensions of the input wind components, and that wind components must be
# either 2D or 3D arrays. The data read in is 3D and has latitude and
# longitude as the last dimensions. The bundled tools can make the process of
# re-shaping the data a lot easier to manage.
uwnd, uwnd_info = prep_data(uwnd, 'tyx')
vwnd, vwnd_info = prep_data(vwnd, 'tyx')
#print uwnd.shape

# It is also required that the latitude dimension is north-to-south. Again the
# bundled tools make this easy.
#lats, uwnd, vwnd = order_latdim(lats, uwnd, vwnd)

# Create a VectorWind instance to handle the computation of streamfunction and
# velocity potential.
#w = VectorWind(uwnd, vwnd)

# Compute the streamfunction and velocity potential. Also use the bundled
# tools to re-shape the outputs to the 4D shape of the wind components as they
# were read off files.
#sf, vp = w.sfvp()
#sf = recover_data(sf, uwnd_info)
#vp = recover_data(vp, uwnd_info)

#print sf.shape
#print sf[1,1,1]
print("Streamfunction done")

#---NetCDF write---------------------------------------------------------------
print("Start NetCDF writing")

ncout = Dataset(fout, 'w', format='NETCDF4')
ncout.description = "Streamfunction form %s" % (fout)

# Using our previous dimension info, we can create the new time dimension
# Even though we know the size, we are going to set the size to unknown

ncout.createDimension(dimnam[0], lons.size)
ncout.createDimension(dimnam[1], lats.size)
ncout.createDimension(dimnam[2], None)


ncout_var = ncout.createVariable(dimnam[0], nc.variables[dimnam[0]].dtype,dimnam[0])
for ncattr in nc.variables[dimnam[0]].ncattrs():
    ncout_var.setncattr(ncattr, nc.variables[dimnam[0]].getncattr(ncattr))
ncout.variables[dimnam[0]][:] = lons
#nc_dim = fid.createVariable('time', fout_id.variables['time'].dtype,\
#                                   ('time',))

#ncatt = nc.ncattrs()
#print(ncatt)

for name in nc.ncattrs():
     print "Global attr", name, "=", getattr(nc,name)

print("spat'")
#print(nc.variables['latitude'].ncattrs())

ncout.close()
