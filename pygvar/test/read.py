from pygvar.reader import parse_gvar
from pygvar.calibration import vis_calibrate, ir_calibrate
from pygvar.geos_nav import GeosNav
from pygvar.navigation import gvarchannels2geos
import math
import os


import matplotlib.pyplot as plt
import pylab


import cartopy.crs as ccrs


def view_gvar(gvar_file, channel=None, calibrate_to=None, reproject2geos=False):
    assert calibrate_to in (None, 'radiance', 'reflectance', 'temperature')
    satellite = os.path.basename(gvar_file).split('.')[0].upper()

    b0, data = parse_gvar(gvar_file, channels_to_extract=[channel])

    chn_data = data[channel]


    if calibrate_to is not None:
        if channel ==  1:
            calib_data = vis_calibrate(counts=chn_data, channel_detectors=b0.detectors[channel],
                                       detector_lines=b0.detector_lines[channel],
                                       satellite=satellite, calibrate_to=calibrate_to)
        else:
            side = 2 if b0.iscan.side2 == 1 else 1
            calib_data = ir_calibrate(counts=chn_data, detector_lines = b0.detector_lines[channel],
                                      satellite=satellite, channel=channel, side=side, calibrate_to=calibrate_to)
    else:
        calib_data = chn_data



    if reproject2geos is True:
        # Create a GeosNav object representing  the image space in normalised geostationary projection
        # located and the channel original subsatellite point
        #full disk field of view in degrees, what the satellite sees from space, part of the Earth disk
        fd_fov_deg = 17.326243724756093  # this is chosen on the basis that 1 pixel will have roughly 1000 m
        # (1002.014) and the FD image 10800 nl/nc
        satellite_height = 42164.365 * 1e3 - 6378.137 * 1e3  # from GVAR NAV  nominal radial distance of instrument (km) -  earth equatorial radius (km)
        orig_nl, orig_nc, yres_km, xres_km = b0.channels_shape[channel]
        res_km = min(yres_km, xres_km)

        # compute the res in microradians by dividing the res in km obtained from block 0 by the satellite height and multiplying
        # by 9 zeros , three for getting meters from km and 6 for getting microrads from rads
        # note the res_mrad is rounded intentionally to the nearest int as per OGE. This means the res in meters will not be a nice integer
        res_mrad = int(round((res_km) / satellite_height * 1e9))
        # compute res in meters and FD lines
        res_metres = abs(res_mrad * satellite_height) / 1e6
        fd_lines = int(round(math.radians(fd_fov_deg) / (res_mrad * 1e-6)))
        lon_0 = b0.subpoint[1]
        navigator = GeosNav(
            lon_0=lon_0, sweep='y', ellipsoid='GRS80', sat_height=satellite_height,
            res_mrad=res_mrad, fd_lines=fd_lines
        )


        calib_data = gvarchannels2geos(raw_channel_data=calib_data,channel=channel,
                                       block0=b0, navigator=navigator,nchunks_per_dim=5)
        geos_proj = ccrs.Geostationary(satellite_height=satellite_height,
                                          central_longitude=lon_0,sweep_axis='y')

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection=geos_proj)
        ax.coastlines(resolution='10m', color='red', linestyle='-', alpha=1)

        ax.set_global()

        ax.imshow(calib_data,extent=navigator.extent, origin='upper', cmap='viridis', vmin=.0, vmax=.4)
        plt.show()


    else:

        pylab.imshow(calib_data, interpolation='nearest')
        pylab.title('JUSSI')
        pylab.show()


if __name__ == '__main__':
    import logging
    logging.basicConfig()
#    import pylab
    logger = logging.getLogger()
    logger.setLevel('INFO')

    src_file = '/data/sample_gvar/goes12.2003.152.223144'
    #src_file = '/data/sample_gvar/goes08.1996.335.234513'
    # view reprojected calibrated reflectance
    d = view_gvar(src_file, channel=1, calibrate_to='reflectance', reproject2geos=True)




