import os.path
import numpy as np
import matplotlib.pyplot as plt
# from mpl_toolkits.basemap import Basemap
# plt.plot([1,2,3,4])
# plt.ylabel('some numbers')
# plt.show()


# m = Basemap(resolution='f',projection='merc',
#             lon_0=160,
#             llcrnrlat=-30.0,
#             urcrnrlat=30.0,
#             llcrnrlon=100.,
#             urcrnrlon=270.0,
#             lat_ts=0.0)


# lats = np.array([ 14.375,  14.125,  14.125,   9.375,  13.625,  14.375,   5.625,
#     13.875,  14.625,   5.875,   8.875,   5.625,   8.875,  13.375,
#     11.125,   8.375,  12.375,   6.125,   5.375,   8.375,   7.375,
#      7.875,  14.375,  14.875,   9.875,  11.125,  14.875,   7.875,
#      9.125,  11.625,   5.125,  10.875,   5.125,  12.125,  12.625,
#      8.625,   5.125,   8.375,  11.625,  11.375,  12.875,  14.375,
#      8.875,   8.375,   6.375,   8.625,   5.875,   8.125,   9.375,
#      5.875,   8.125,   8.875,   5.375,   8.875,   5.625,  11.875,
#      9.875,   9.875,  10.875,  11.375,   9.875,   9.375,  13.125,
#     14.125,   8.125,  14.875,   9.875,   9.625,  10.625,  12.125,
#      9.375,   5.625,  13.625,   6.375,  10.125,  14.875,   8.875,
#      5.125,   6.125,   8.625,   6.875,   9.375,   9.125,   9.375,
#     13.625,   6.125,   6.125,   6.875,  11.375,  13.375,  10.375,
#      6.875,   7.625,   7.625,  13.875,   5.125,   6.125,  14.125,
#      7.375,   5.375])
# lons = np.array([ 122.125,  125.875,  122.375,  132.125,  124.375,  122.625,
#     206.125,  122.625,  120.375,  187.625,  243.625,  161.625,
#     220.375,  121.875,  125.625,  130.125,  219.625,  223.875,
#     223.375,  190.125,  132.875,  185.875,  122.625,  120.875,
#     259.375,  125.875,  204.875,  131.875,  130.875,  125.375,
#     234.375,  241.375,  188.125,  124.375,  124.375,  140.375,
#     205.625,  130.875,  256.875,  239.375,  124.125,  123.125,
#     131.625,  126.375,  241.125,  130.875,  233.625,  184.375,
#     205.125,  169.375,  150.375,  183.375,  168.875,  239.625,
#     180.125,  241.375,  131.375,  244.875,  244.375,  125.625,
#     131.375,  150.125,  203.625,  122.125,  243.125,  159.125,
#     125.625,  243.375,  125.875,  127.375,  130.125,  146.625,
#     203.125,  185.125,  204.625,  120.625,  130.125,  233.875,
#     131.875,  258.625,  130.125,  173.375,  258.125,  129.125,
#     137.875,  215.625,  214.125,  234.625,  125.375,  136.625,
#     231.125,  145.625,  128.625,  232.875,  125.125,  145.125,
#     239.875,  138.625,  241.375,  169.625])


#read rays
season = np.array(["DJF"])
period =  np.array(["Inf"])
lons = np.arange(60,360,90)
lats = np.arange(80,-85,-5)
krange = np.arange(1,7)
years = np.arange(1980,1986)
roots = np.arange(1,4)

# number of dots per ray
raylength = 21
nlines = np.arange(raylength)


#create arrays
miss = -9999.
hr = np.full((nlines.size,roots.size),miss)

day = np.full_like(hr,miss)
x = np.full_like(hr,miss)
y = np.full_like(hr,miss)
rlon = np.full_like(hr,miss)
rlat = np.full_like(hr,miss)
rk = np.full_like(hr,miss)
rl = np.full_like(hr,miss)


for ssn in season :
    for ip in period :
        # print ip, ssn
        for lat in lats :
            for lon in lons :
                for k in krange :
                    fout = "../output/matlab/yearly/raydens_{}{:d}_{:d}_{:d}N_{:d}E_period{}_k{:d}".format(ssn,years(0),years(-1),lat,lon,ip,k)
                    for yr in years :
                        for root in roots :
                            fin = "../output/matlab/yearly/ray_{}{:d}_{:d}N_{:d}E_period{}_k{:d}_root{:d}".format(ssn,yr,lat,lon,ip,k,root)
                            if not os.path.isfile(fin) :
                                # print fin, ": File is missing "
                                continue
                            else :
                                # print fin, ": File exists"
                                print fin
                                fray = open(fin,'r')

                                for l in nlines :
                                    print l
                                    line = fray.readline()
                                    if line == '':
                                     print ' fray is at the eof'
                                     break
                                    else :
                                     print ' line: ', line.split()
                                     columns = line.split(',')

                                     hr[l,root]=columns[1]
                                     day[l,root]=columns[2]
                                     x[l,root]=columns[3]
                                     y[l,root]=columns[4]
                                     rlon[l,root]=columns[5]
                                     rlat[l,root]=columns[6]
                                     rk[l,root]=columns[7]
                                     rl[l,root]=columns[8]



                                quit()





# nx, ny = 10, 3
#
# # compute appropriate bins to histogram the data into
# lon_bins = np.linspace(lons.min(), lons.max(), nx+1)
# lat_bins = np.linspace(lats.min(), lats.max(), ny+1)
#
# # Histogram the lats and lons to produce an array of frequencies in each box.
# # Because histogram2d does not follow the cartesian convention
# # (as documented in the numpy.histogram2d docs)
# # we need to provide lats and lons rather than lons and lats
# density, _, _ = np.histogram2d(lats, lons, [lat_bins, lon_bins])
#
# # Turn the lon/lat bins into 2 dimensional arrays ready
# # for conversion into projected coordinates
# lon_bins_2d, lat_bins_2d = np.meshgrid(lon_bins, lat_bins)
#
# print lon_bins
# print lon_bins_2d
# print lat_bins_2d
# print lat_bins
# print density
#
#
# # # convert the xs and ys to map coordinates
# # xs, ys = m(lon_bins_2d, lat_bins_2d)
#
# # plt.pcolormesh(xs, ys, density)
# # plt.colorbar(orientation='horizontal')
# #
# # # overlay the scatter points to see that the density
# # # is working as expected
# # plt.scatter(*m(lons, lats))
# #
# # plt.show()
