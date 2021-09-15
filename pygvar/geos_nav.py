"""
Navigation functions for  normalised geostationary projection.



The Geostationary satellite view is a  vertical  perspective projection with a fixed latitude (0 degrees).
This represents a view or the Earth  from space or a vantage point (satellite) where the projection  surface is a plane
tangent to the earth sphere at the surface position directly below the instrument.
The view contains the whole earth and the instruments produces an image associated with each viewing, a process that can
 also be called scanning.

Scanning

The satellite instrument scans the Earth in a given predefined pattern. The scanning process happens
 in  blocks  or chunks of space  on the vertical axis (lines) with entire horizontal chunks being scanned
  in a very short interval of time (miliseconds).
This constitutes a scan swath.
A scan swath generates a chunk of image usually a couple of lines  where each line consists  usually of all pixels the
instrument can produce in  the horizontal direction.
The internal sensors or detectors pan horizontally across in small steps . These steps are expressed in radians or
microradians because they represent increments or offsets in the  horizontal field of view of the instrument which is
 defined in radians. Similarly, on the vertical axis while performing a scan swath the data points are collected or
 sampled at specific offsets. These offsets are usually the same size as the horizontal ones depicted above.

The sum of data points collected in one scan instance comprises the whole image. This image can be perceived as a
virtual fixed grid with dimensions equal to the instrument field of view and an identical spatial separation measured in
microradians in both axes.



This fixed grid is usually rectified to an ellipsoid viewed from the idealized geostationary projection. This ellipsoid
has to be used when geo-referencing data points from the fixed grid (transforming from radians to geographic coordinates)


The instrument field of view and the stepping in radians defines the number of pixels in the acquired image.
fd_lines = field of view in radians / stepping

These relationship is very important as it allows to convert easily between the integer indices on the fixed grid and the radian values.

Furthermore the scanning angles from the fixed grid are related to the projections coordinates of the normalised geostatioanary
satellite projection and the height of the instrument above ground :

scanning angle (radians) = projection_coordinate(meters) / satellite height


Note that from a satellite perspective all pixels produced  are valid (they have been scanned after all), but not all of
them are on the disk of the earth.
Some pixels are outside the earth disk and have only radians or meters but no geographic corresponding coordinates.

Navigation

Generally speaking the navigation or the process of navigating an image means the process of converting the
radians corresponding to a given  image pixels to their earth locations.
The fixed grid coordinates, N/S elevation angle and E/W scanning angle, coupled
with the location of the satellite and the parameters associated with the selected earth model (GRS80 or WGS84 or custom) are
used to determine the geodetic latitude/longitude coordinates.




The transformation employs an Earth-Centered Fixed coordinate system  with its origin located at the centre of the Earth
and a satellite coordinate frame with it's origin at the mass center of the satellite

Detailed formulas for these transformation are available in many places across the internet, for example
for GOES R at
https://www.goes-r.gov/products/docs/PUG-L2+-vol5.pdf
page  21


Due to the fact that different satellite are build by different vendors the images produced by them come in various formats
including various forms of navigation. Certain data vendors provide users with navigation coefficients (CFAC, LFAC, COFF, LOFF)
that relate the image pixels to fixed grid locations.



l = LOFF + nint( y * 2 **-16 * LFAC
c = COFF + nint( x * 2 ** -16 * CFAC)


where CFAC and LFAC are the scaling factors and COFF LOFF are the offsets.


the scaling factors are related to the fixed grid stepping or resolution using following relationship

L/CFAC = -2**16 / resolution in microrad.

thus

resolution in microrad.= -2 ** 16 / float(np.rad2deg(L/CFAC))
and
L/COFF = fd_lines / 2


From here it is easy to calculate  the full disk field of view and navigate these images.



"""
import numpy as np
import pyproj
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
logger = logging.getLogger(__name__)


def pj4_setup(r_eq=None, r_pol=None, sat_height=None):
    """
    Utility function used by the forward projection pj4_ll2yx and the inverse projection pj4_yx2ll
    Takes input the defintioo of the ellipsoid and satellite height and returns a series of parameters used in the
    navigation process
    :param r_eq: number, equatorial radius of earth in meters or km
    :param r_pol: number, polar radus of Earth in meters or numbers
    :param sat_height: number, satellite height above  ground in meters or km
    :return: sa tuple of coeffients used by the forward and inverse fucntions
    """

    rad_factor = r_eq / sat_height  # this is a factor that
    es = (r_eq ** 2 - r_pol ** 2) / r_eq ** 2
    radius_g_1 = sat_height / r_eq
    radius_g = 1 + radius_g_1
    C = radius_g * radius_g - 1.
    if r_eq != r_pol:  # eccentricity condition, geos navigation does not use the spheroid functions
        one_es = 1 - es
        rone_es = 1. / one_es
        radius_p = np.sqrt(one_es)
        radius_p2 = one_es
        radius_p_inv2 = rone_es
    else:
        radius_p = radius_p2 = radius_p_inv2 = 1.

    return rad_factor, radius_p, radius_p2, radius_g, radius_g_1, radius_p_inv2, C


def pj4_ll2yx(lat=None, lon=None, lon_0=None, sweep=None, r_eq=None, r_pol=None, sat_height=None):
    """
        Given a geostationary satellite located at sat_height and lon_0 longitude converts the  latitude and longitude
        coordinates represented on an ellipsoid defined by r_eq and  r_pol to  satellite intermediary coordinates (radians)
        Proj4 CPP based. Ported for optimisation and control reasons
        :param lat: number or numpy array
        :param lon: number or numpy array
        :param lon_0: number
        :param sweep: x or y
        :param r_eq: number
        :param r_pol: number
        :param sat_height: number
        :return: y and x (number or numpy array depending on the input)

        NB: The PROJ4 scales the intermediary coordinates for computational reasons or other reasons
        unknown to us. I have disabled the scaling so the y and x coords are identical to the coords used in CGMS
        and other official docs

    """

    lat = np.asarray(lat)
    lon = np.asarray(lon)
    assert lat.ndim == lon.ndim, 'Input lat and lon have to have same dimensions'
    # store original dimension of inut data
    ondim = lat.ndim

    if ondim == 1:
        # use  broadcasting
        lon = lon[np.newaxis, :]
        lat = lat[:, np.newaxis]
    else:
        lat = np.atleast_2d(lat)
        lon = np.atleast_2d(lon)
    if (lat.ndim > 2) or (lon.ndim > 2):
        raise Exception('Invalid dimensions for lat %s or lon %s. Max allowed is 2' % (lat.ndim, lon.ndim))

    check_nav_constants(r_eq=r_eq, r_pol=r_pol, sat_height=sat_height, sweep=sweep, lon_0=lon_0)
    rad_factor, radius_p, radius_p2, radius_g, radius_g_1, radius_p_inv2, C = pj4_setup(r_eq=r_eq, r_pol=r_pol,sat_height=sat_height)
    m = ~(np.isnan(lat) | np.isnan(lon))
    input_masked =  m[~m].size > 0
    if input_masked:
        lat = lat[m]
        lon = lon[m]


    # compute relative lon from lon_0
    rel_lon = lon-lon_0
    # connvert to radians
    lam = np.radians(rel_lon)
    phi = np.radians(lat)
    # geocentric latitiude
    phi = np.arctan(radius_p2 * np.tan(phi))
    # Calculation of the three components of the vector from satellite to position on earth surface (lat, lon)
    r = radius_p / np.hypot(radius_p * np.cos(phi), np.sin(phi))
    Vx = r * np.cos(lam) * np.cos(phi)
    Vy = r * np.sin(lam) * np.cos(phi)
    Vz = r * np.sin(phi)
    visible = ((radius_g - Vx) * Vx - Vy * Vy - Vz * Vz * radius_p_inv2) > 0.
    # Calculation based on view angles from satellite
    tmp = radius_g - Vx

    # compared to the original proj4 equivalent function "e_forward from geos.cpp in proj4 source" this function
    # does not scale(multiply) the radians with radius_g_1
    # in this case if the radians are multiplied with sat_height we get meters directkly.
    # in the proj4 original code the radians had to be multiplied by the inverse of radius_g1

    if sweep == 'x':
        # x = radius_g_1 * np.arctan(Vy / np.hypot(Vz, tmp))
        x = np.arctan(Vy / np.hypot(Vz, tmp))
        # y = radius_g_1 * np.arctan(Vz / tmp)
        y = np.arctan(Vz / tmp)
    else:
        # x = radius_g_1 * np.arctan(Vy / tmp)
        x = np.arctan(Vy / tmp)
        # y = radius_g_1 * np.arctan(Vz / np.hypot(Vy, tmp))
        y = np.arctan(Vz / np.hypot(Vy, tmp))

    # return is complicated in order to honour the dimensionality of inputs combined with the visible masking
    if ondim == 0:
        if not visible: return np.nan, np.nan
        return y.item(), x.item()
    else:
        y[~visible] = np.nan
        x[~visible] = np.nan
        if input_masked:
            ry = np.empty(m.shape, dtype=y.dtype)
            ry[:] = np.nan
            ry[m] = y
            rx = np.empty(m.shape, dtype=x.dtype)
            rx[:] = np.nan
            rx[m] = x
        else:
            ry = y
            rx = x
        if ondim == 1:
            ry = np.atleast_1d(np.squeeze(ry))
            rx = np.atleast_1d(np.squeeze(rx))

        return ry, rx


def pj4_yx2ll(y=None, x=None, lon_0=None, sweep=None, r_eq=None, r_pol=None, sat_height=None):
    """
        Given a geostationary satellite located at sat-height and lon_0 longitude convert the input intermediary
        coordinates (radians) to latitude and longitude  represented on an ellipsoid defined by r_eq  and  r_pol
        :param y: number or numpy array
        :param x: number or numpy array
        :param lon_0: number
        :param r_eq: number
        :param r_pol: number
        :param sat_height: number
        :return: lat and lon (number or numpy array depending on the input)

    NB: The PROJ4 scales the intermediary coordinates for computational reasons or other reasons
        unknown to us. I have disabled the scaling so the y and x coords are identical to the coords used in CGMS
        and other official docs
    """

    x = np.asarray(x)
    y = np.asarray(y)
    ondim = x.ndim
    x = np.atleast_2d(x)
    y = np.atleast_2d(y)
    if (x.ndim > 2) or (y.ndim > 2):
        raise Exception('Invalid dimensions for x %s or y %s. Max allowed is 2' % (x.ndim, y.ndim))

    check_nav_constants(r_eq=r_eq, r_pol=r_pol, sat_height=sat_height, sweep=sweep, lon_0=lon_0)
    rad_factor, radius_p, radius_p2, radius_g, radius_g_1, radius_p_inv2, C = pj4_setup(r_eq=r_eq, r_pol=r_pol,
                                                                                        sat_height=sat_height)

    # Setting three components of vector from satellite to position
    # compared to the original proj4 equivalent function "e_inverse from geos.c" this function
    # does not divide the radians with radius_g_1 because also the inverse function pj_ll2yx did not scale them

    Vx = -1.0
    if sweep == 'x':
        # Vz = np.tan(y / radius_g_1)
        Vz = np.tan(y)
        # Vy = np.tan(x / radius_g_1 ) * np.hypot(1., Vz)
        Vy = np.tan(x) * np.hypot(1., Vz)
    else:
        # Vy = np.tan(x / radius_g_1)
        Vy = np.tan(x)
        # Vz = np.tan(y / radius_g_1 ) * np.hypot(1.0, Vy)
        Vz = np.tan(y) * np.hypot(1.0, Vy)
    # Calculation of terms in cubic equation and determinant
    a = Vz / radius_p
    a = Vy * Vy + a * a + Vx * Vx
    b = 2 * radius_g * Vx
    det = b * b - 4 * a * C

    # handle visbility in numpy style
    visible = det > 0
    if np.any(~visible):
        k = np.zeros_like(det)
        k[:] = np.nan
        a = a[visible]
        det = det[visible]
        k[visible] = (-b - np.sqrt(det)) / (2 * a)
    else:
        k = (-b - np.sqrt(det)) / (2 * a)

    # Calculation of three components of vector from satellite to position
    Vx = radius_g + k * Vx
    Vy *= k
    Vz *= k
    # Calculation of longitude and latitude
    lam = np.arctan2(Vy, Vx)
    phi = np.arctan(Vz * np.cos(lam) / Vx)
    phi = np.arctan(radius_p_inv2 * np.tan(phi))
    lat = np.degrees(phi)
    # adjust relative longitude with lon_0
    lon = np.degrees(lam) + lon_0

    # return is complicated in order to honour the dimensionality of inputs combined with the visible masking

    if ondim == 0:
        if not visible: return np.nan, np.nan
        return lat.item(), lon.item()
    else:
        # lat[notvisible] = np.nan
        # lon[notvisible] = np.nan
        if ondim == 1:
            return np.atleast_1d(lat.squeeze()), np.atleast_1d(lon.squeeze())
        else:
            return lat, lon


def proj4_ll2yx(lat=None, lon=None, lon_0=None, sweep=None, r_eq=None, r_pol=None, sat_height=None):
    """
    PROJ4 API geostationary navigation. It could be used but it was created for testing the correctness of
    the navigation
    :param lat:
    :param lon:
    :param lon_0:
    :param sweep:
    :param r_eq:
    :param r_pol:
    :param sat_height:
    :return:
    """


    proj4_string = f'+proj=geos +a={r_eq:.1f} +b={r_pol:.1f} +h={sat_height:.1f} +lon_0={lon_0:.2f} +units=m +sweep={sweep}'
    geosp = pyproj.Proj(proj4_string)
    try:  # fix the bug in Proj4 that Tomas run into. The bug is related to not
        shp = lon.shape
        X, Y = geosp(list(lon.flatten()), list(lat.flatten()))
        X = np.array(X).reshape(shp)
        Y = np.array(Y).reshape(shp)
        return Y, X
    except AttributeError:
        X, Y = geosp(lon, lat)
        return Y / sat_height, X / sat_height


def proj4_yx2ll(y=None, x=None, lon_0=None, sweep=None, r_eq=None, r_pol=None, sat_height=None):
    """
    PROJ4 API geostationary navigation. It could be used but it was created for testing the correctness of
    the navigation
    :param y:
    :param x:
    :param lon_0:
    :param sweep:
    :param r_eq:
    :param r_pol:
    :param sat_height:
    :return:
    """

    Y = y * sat_height
    X = x * sat_height
    proj4_string = f'+proj=geos +a={r_eq:.1f} +b={r_pol:.1f} +h={sat_height:.1f} +lon_0={lon_0:.2f} +units=m +sweep={sweep}'
    geosp = pyproj.Proj(proj4_string)
    lon, lat = geosp(X, Y, inverse=True)
    return lat, lon


def check_nav_constants(r_eq=None, r_pol=None, sat_height=None, sweep=None, lon_0=None):
    for arg_name, arg_value in locals().items():
        assert arg_value is not None, f'Invalid {arg_name} with value {arg_value}'


def yx2YX(y=None, x=None, geos_obj=None, sat_height=None ):
    """
    Convert native geostationary coordinates (radians) to  geostationary coordinates (meters)

    :param y:
    :param x:
    :param sat_height:
    :param geos_obj:
    :return:
    """
    satheight = geos_obj.sat_height if geos_obj is not None else sat_height
    assert satheight is not None, 'sat_height or a valid geos_obj is required'
    return y*satheight, x*satheight

def YX2yx(Y=None, X=None, geos_obj=None, sat_height=None  ):
    """
    Convert geostationary coordinates (meters) to native geostationary coordinates (radians)

    :param Y:
    :param X:
    :param sat_height:
    :param geos_obj:
    :return:
    """
    satheight = geos_obj.sat_height if geos_obj is not None else sat_height
    assert satheight is not None, 'sat_height or a valid geos_obj is required'
    return Y/satheight, X/satheight



def YX2lc(Y=None, X=None, geos_obj=None, sat_height=None, res_mrad=None, fd_lines=None, fd_fov=None,
          fd_geostransform=None, line_adj=0, col_adj=0):
    """
    Convert projected geostationary coordinates (meters) to line/column indices either using native radian coordinates
    of using the GDAL fd_geotransform
    :param Y:
    :param X:
    :param geos_obj:
    :param res_mrad:
    :param fd_lines:
    :param fd_fov:
    :param geostransform
    :return:
    """

    if geos_obj is not None:
        sat_height = geos_obj.sat_height
        res_mrad = geos_obj.res_mrad
        fd_lines = geos_obj.fd_lines
        fd_fov = geos_obj.fd_fov
        line_adj = geos_obj.line_adj
        col_adj = geos_obj.line_adj

    if fd_geostransform is not None:
        xmin, xres, _, ymax, __, yres = fd_geostransform
        yres = float(yres)
        xres = float(xres)
        ty = ((Y - ymax) / yres) - line_adj
        tx = ((X - xmin) / xres) - col_adj
        return np.floor(ty).astype(np.uint32), np.floor(tx).astype(np.uint32)

    y, x = YX2yx(Y=Y, X=X, sat_height=sat_height)

    return yx2lc(y=y, x=x, res_mrad=res_mrad, fd_lines=fd_lines, fd_fov=fd_fov, line_adj=line_adj, col_adj=col_adj)


def ll2yx(lat=None, lon=None, lon_0=None, sweep=None, r_eq=None, r_pol=None, sat_height=None):
    """
    Given a geostationry satellite located at SAT_HEIGHT and lon_0 longitude convert the  latitude and longitude coordinates
     represended on an ellipsoid defined by EQUATORIAL_RADIUS and  POLAR_RADIUS to  satellite intermediary coordinates (radians)

    :param lat: number or numpy array
    :param lon: number or numpy array
    :param lon_0: number
    :param sweep: x or y
    :param r_eq: number
    :param r_pol: number
    :param sat_height: number
    :return: y and x (number or numpy array depending on the input)
    """
    return pj4_ll2yx(lat=lat, lon=lon, lon_0=lon_0, sweep=sweep, r_eq=r_eq, r_pol=r_pol,sat_height=sat_height)

def ll2YX(lat=None, lon=None, lon_0=None, sweep=None, r_eq=None, r_pol=None, sat_height=None):
    y, x = pj4_ll2yx(lat=lat, lon=lon, lon_0=lon_0, sweep=sweep, r_eq=r_eq, r_pol=r_pol, sat_height=sat_height)
    return yx2YX(y=y, x=x,sat_height=sat_height)

def yx2ll(y=None, x=None, lon_0=None, sweep=None, r_eq=None, r_pol=None, sat_height=None):
    """
        Given a geostationary satellite located at SAT_HEIGHT and lon_0 longitude convert the input intermediary coordinates (radians)
        to latitude and longitude  represented on an ellipsoid defined by EQUATORIAL_RADIUS and  POLAR_RADIUS
        :param y: number or numpy array
        :param x: number or numpy array
        :param lon_0: number
        :param sweep: x or y
        :param r_eq: number
        :param r_pol: number
        :param sat_height: number
        :return: lat and lon (number or numpy array depending on the input)
        """

    check_nav_constants(lon_0=lon_0, sweep=sweep, r_eq=r_eq, r_pol=r_pol, sat_height=sat_height)
    return pj4_yx2ll(y=y, x=x, lon_0=lon_0, sweep=sweep, r_eq=r_eq, r_pol=r_pol, sat_height=sat_height)



def yx2lc(y=None, x=None, res_mrad=None, fd_lines=None, fd_fov=None, line_adj=0, col_adj=0):
    """
    Convert  radians or intermediary coordinates of an ideal geostationary projection associated with and image
    generated by a geostationary instrument into  line and columns indices
    Some parameters are complementary, for example the res_mrad is the most important. If not supplied the
    fd_fov and fd_lines become mandatory and are used to compute the res_mrad
    :param y: numpy array or number
    :param x: numpy array or number
    :param res_mrad: number the resoluation of the instrument in microradians
    :param fd_lines: number, int, the number of pixels in a full disk image in one dimension.
    :param fd_fov: number, the instrument full disk field of view in radians
    :param line_adj, number, image N-S shift in pixels wrt to the idealised geostationary projection
    :param col_adj, number, image W-E shift in pixels wrt to the idealised geostationary projection

    :return: tuple of int numbers or int32 numpy arrays representing the indices of the input radians
    """
    if res_mrad is None:
        assert fd_fov not in (None, np.nan), f'Invalid fd_fov {fd_fov}'
        assert fd_lines not in (None, np.nan), f'Invalid fd_lines {fd_lines}'
        res_mrad = fd_fov / fd_lines * 1e6

    if fd_lines is None:
        assert fd_fov not in (None, np.nan), f'Invalid fd_fov {fd_fov}'
        assert res_mrad not in (None, np.nan), f'Invalid res_mrad {res_mrad}'
        fd_lines = fd_fov / res_mrad * 1e-6

    xres_rad = res_mrad * 1e-6
    yres_rad = -xres_rad

    offset = fd_lines // 2 #+ .5 #account for rounding, normally I shouod have added .5
    #to each index but is is easier to add it below
    # Taking the absolute value before flooring is necessary to avoid negative indices. In some cases the values are
    # near the machine precision of a float which can cause them to become negative.
    # l_ind and c_ind always have to be positive integers

    # the line_adj and col_adj relate to adjustment applied to l,c indices when being converted to coordinates
    # this is the opposite case, so we need to subtract the adj values

    l_ind = np.floor(np.abs(offset - line_adj + y / yres_rad)).astype(np.int32)
    c_ind = np.floor(np.abs(offset - col_adj + x / xres_rad)).astype(np.int32)

    return l_ind, c_ind


def lc2yx(l=None, c=None, res_mrad=None, fd_fov=None, fd_lines=None, line_adj=0, col_adj=0):
    """
    Convert line and column indices associated with and image generated by a geostationary instrument
    into  radians of an ideal normalised geostationary projection.

    Some parameters are complementary, for example the res_mrad is the most important. If not supplied the
    fd_fov and fd_lines become mandatory and are used to compute the res_mrad
    :param l: numpy array or number, int
    :param c: numpy array or number, int
    :param res_mrad: number the resoluation of the instrument in microradians
    :param fd_lines: number, int, the number of pixels in a full disk image in one dimension.
    :param fd_fov: number, the instrument full disk field of view in radians
    :param line_adj, number, image N-S shift  in pixels wrt to the idealised geostationary projection
    :param col_adj, number, image W-E shift in pixels wrt to the idealised geostationary projection

    :return: tuple of float numbers or float numpy arrays representing the radians of the input indices
    """

    if res_mrad is None:
        assert fd_fov not in (None, np.nan), f'Invalid fd_fov {fd_fov}'
        assert fd_lines not in (None, np.nan), f'Invalid fd_lines {fd_lines}'
        res_mrad = fd_fov / fd_lines * 1e6

    if fd_fov is None:
        assert fd_lines not in (None, np.nan), f'Invalid fd_lines {fd_lines}'
        assert res_mrad not in (None, np.nan), f'Invalid res_mrad {res_mrad}'

    xres_rad = res_mrad * 1e-6
    yres_rad = -xres_rad
    offset = fd_lines // 2
    # fact y space  goes  from + to  -, x space goes from - to +
    # yres is negative, xres is positive
    # half yres has to be added (substracted in reality)
    # half x res has to be added
    # this is all because we are converting lower precision coords (l, c) to high precision coords (y, x)
    # whenever this kind of transformation is involved this is the correct approach

    ty = (l - offset + line_adj) * yres_rad + yres_rad / 2
    tx = (c - offset + col_adj) * xres_rad + xres_rad / 2

    return ty, tx


def ll2lc(lat=None, lon=None, geos_obj=None, lon_0=None, sweep=None, r_eq=None, r_pol=None, sat_height=None,
          res_mrad=None, fd_lines=None, fd_fov=None, line_adj=0, col_adj=0):
    """
    Convenience method to convert geographic coordinates (lat, lon)
    to indices of an image in an ideal normalised geostationary projection

    :param lat: numebr pr numy array
    :param lon: number or numpy array
    :param geos_obj: an instance of GeosNav class containing all required satellite scene parameters
    :param lon_0: number the instrument SSP longitude
    :param sweep: string, y or x
    :param r_eq: number, equatorial radius
    :param r_pol: number, polar radius
    :param sat_height: number, satellite height above ground in meters
    :param res_mrad: number, float, the resolution of the image in microradians
    :param fd_lines: number int,
    :param fd_fov: number, the full disk field of view in radians
    :param: line_adj, number that is added to line indices computed inside yx2lc function
    :param: col_adj, number that is added to column indices computed inside yx2lc function


    :return:
    """

    if geos_obj is None:
        y, x = ll2yx(lat=lat, lon=lon, lon_0=lon_0, sweep=sweep, r_eq=r_eq,
                     r_pol=r_pol, sat_height=sat_height)
        return yx2lc(y=y, x=x, res_mrad=res_mrad, fd_lines=fd_lines, fd_fov=fd_fov, line_adj=line_adj, col_adj=col_adj)
    else:
        y, x = ll2yx(lat=lat, lon=lon, lon_0=geos_obj.lon_0, sweep=geos_obj.sweep, r_eq=geos_obj.r_eq,
                     r_pol=geos_obj.r_pol, sat_height=geos_obj.sat_height)
        return yx2lc(y=y, x=x, res_mrad=geos_obj.res_mrad, fd_lines=geos_obj.fd_lines, fd_fov=geos_obj.fd_fov,
                     line_adj=geos_obj.line_adj, col_adj=geos_obj.col_adj)


def lc2ll(l=None, c=None, geos_obj=None, lon_0=None, sweep=None, r_eq=None, r_pol=None, sat_height=None,
          res_mrad=None, fd_lines=None, fd_fov=None, line_adj=0, col_adj=0):
    """
    Convenience method to convert indices of an image in an ideal normalised geostationary projection to
    geographic coordinates (lat, lon)

    :param l: number or numpy array, int, image line indices
    :param c: number or numpy array, int image columns indices
    :param geos_obj: an instance of GeosNav class containing all required satellite scene parameters
    :param lon_0: number the instrument SSP longitude
    :param sweep: string, y or x
    :param r_eq: number, equatorial radius
    :param r_pol: number, polar radius
    :param sat_height:   number, satellite height above ground in meters
    :param res_mrad: number, float, the resolution of the image in microradians
    :param fd_lines: number int,
    :param fd_fov: number, the full disk field of view in radians
    :return:
    """

    if geos_obj is None:
        y, x = lc2yx(l=l, c=c, res_mrad=res_mrad, fd_fov=fd_fov, fd_lines=fd_lines, line_adj=line_adj, col_adj=col_adj)
        return yx2ll(y=y, x=x, lon_0=lon_0, sweep=sweep, r_eq=r_eq, r_pol=r_pol, sat_height=sat_height)
    else:
        y, x = lc2yx(l=l, c=c, res_mrad=geos_obj.res_mrad, fd_fov=geos_obj.fd_fov, fd_lines=geos_obj.fd_lines,
                     line_adj=geos_obj.line_adj, col_adj=geos_obj.col_adj)
        return yx2ll(y=y, x=x, lon_0=geos_obj.lon_0, sweep=geos_obj.sweep, r_eq=geos_obj.r_eq,
                     r_pol=geos_obj.r_pol, sat_height=geos_obj.sat_height)


def lc2YX(l=None, c=None, res_mrad=None, fd_fov=None, fd_lines=None, sat_height=None, line_adj=0, col_adj=0):

    y, x = lc2yx(l=l, c=c, res_mrad=res_mrad, fd_fov=fd_fov, fd_lines=fd_lines, line_adj=line_adj, col_adj=col_adj)
    return yx2YX(y=y, x=x, geos_obj=None, sat_height=sat_height)


def reproject_geos2geos_chunk(sl=None, el=None, sc=None, ec=None,
                              src_geos_nav=None, dst_geos_nav=None,
                              ):
    """
    Reproject indices of a geostationary projection defined by src_geos_nav  indices  of a geostationary projection
    defined by dst_geos_nav  in the interval sl-el, sc-ec

    :param sl: number, start line
    :param el: number, end line
    :param sc: number, start column
    :param ec: number, end column
    :param src_geos_nav: an instance of GeosNav class containing all required satellite scene parameters (source)
    :param dst_geos_nav: an instance of GeosNav class containing all required satellite scene parameters (destination)
    :return: numpy array


    """

    # create destination lines and cols
    dst_l, dst_c = np.mgrid[sl:el:1, sc:ec:1]
    # # accouont for subset
    dst_l += src_geos_nav.start_line
    dst_c += src_geos_nav.start_col
    lats, lons = lc2ll(l=dst_l, c=dst_c, geos_obj=dst_geos_nav)
    #convert ll to lc in source prj
    src_l, src_c = ll2lc(lat=lats, lon=lons, geos_obj=src_geos_nav)
    #mask
    lm = (src_l > src_geos_nav.start_line) & (src_l < src_geos_nav.end_line)
    cm = (src_c > src_geos_nav.start_col) & (src_c < src_geos_nav.end_col)
    m = lm & cm

    # account for subset
    src_l -= src_geos_nav.start_line
    src_c -= src_geos_nav.start_col


    return sl, el, sc, ec, m, src_l[m], src_c[m]





def reproject_geos2geos_mp(src_geos_array=None, src_geos_nav=None,
                           dst_geos_nav=None, dst_fill_value=None,
                           dst_dtype=None, return_indices=False):
    """
    Reproject an image from a geostationary projection defined by src_geos_nav  to an image in geostationary projection
    defined by dst_geos_nav arguments in parallel

    :param src_geos_array: numpy array
    :param src_geos_nav: an instance of GeosNav class containing all required satellite scene parameters (source)
    :param dst_geos_nav: an instance of GeosNav class containing all required satellite scene parameters (destination)
    :param dst_fill_value, number representing the vallue the dst array will be filled prior to be reprojected
    :param dst_dtype: numpy dtype fo the destination image
    :param return_indices, bool, defaults to False. Specify whether the indices used in reprojection should be returned as well
    :return: numpy array
    """

    src_nl, src_nc = src_geos_array.shape
    dest_nl = src_nl
    dest_nc = src_nc
    dest_dtype = dst_dtype or src_geos_array.dtype

    # allocate memory for dest
    dest_geos_data = np.full_like(src_geos_array, dtype=dest_dtype, fill_value=dst_fill_value)

    if return_indices:
        l = np.zeros((dest_nl, dest_nc), dtype='i2')
        c = np.zeros((dest_nl, dest_nc), dtype='i2')
        m = np.zeros((dest_nl, dest_nc), dtype=np.bool)
    # split the full disk space into several chunks
    nsegs = 5 if src_nl < 5000 else 8
    lseg_len = int(src_nl / nsegs)
    cseg_len = int(src_nc / nsegs)

    with ProcessPoolExecutor(max_workers=nsegs) as pexec:
        futures = list()
        for i in range(nsegs):
            sl = i * lseg_len
            el = sl + lseg_len
            if el > src_nl: el = src_nl
            for j in range(nsegs):
                sc = j * cseg_len
                ec = sc + cseg_len
                if ec > src_nc: ec = src_nc
                fut = pexec.submit(reproject_geos2geos_chunk,
                                   sl=sl, el=el, sc=sc, ec=ec,
                                   src_geos_nav=src_geos_nav,
                                   dst_geos_nav=dst_geos_nav,

                                   )
                futures.append(fut)
        try:

            for f in as_completed(futures,timeout=60*3):
                result = f.result()
                sl, el, sc, ec, cm, ml, mc = result
                dest_geos_data[sl:el, sc:ec][cm] = src_geos_array[ml, mc]
                if return_indices:
                    l[sl:el, sc:ec][cm] = ml
                    c[sl:el, sc:ec][cm] = mc
                    m[sl:el, sc:ec][cm] = cm[cm]

        except KeyboardInterrupt:
            #be aware it is NOT possible to stop running processes nicely.
            #the script will usually terminate after some time. This is an improvement over python2 where the script
            # could just hand for days
            #in some exceptional circumstances
            pexec.shutdown(wait=False)
            for f in futures:
                done = f.done()
                if not done:
                    r = f.cancelled()
                    if not r:
                        logger.info(f'{f} is still reprojecting a chunk of image and will take some time before it finishes')

            raise
    if not return_indices:
        return dest_geos_data
    else:
        return dest_geos_data, m, l, c




def reproject_geos2geos(
        src_geos_array=None, src_geos_nav=None,
        dst_geos_nav=None, dst_dtype=None, dst_fill_value=None,
        return_indices=False
    ):
    """
    Reproject an image from a geostationary projection defined by src_geos_obj  to an image in geostationary projection
    defined by dst_geosz_obj

    :param src_geos_array: numpy array
    :param src_geos_nav: an instance of GeosNav class containing all required satellite scene parameters (source)
    :param dst_geos_nav: an instance of GeosNav class containing all required satellite scene parameters (destination)
    :param dst_dtype: numpy dtype of the destination image
    :param dst_fill_value, number representing the value the dst array will be filled prior to be reprojected
    :param return_indices, bool defaults to false, if True also return the  indices computed during re-projection
    :return: numpy array
    """

    src_nl, src_nc = src_geos_array.shape
    dest_nl = src_nl
    dest_nc = src_nc
    dest_dtype = dst_dtype or src_geos_array.dtype

    # allocate memory for dest
    dest_geos_data = np.full_like(src_geos_array, dtype=dest_dtype, fill_value=dst_fill_value)

    # the transformation needed is lc from destination projection to ll
    # the ll are then converted to lc in source projection and the src image is indexed with this src lc

    # create destination lines and cols
    dst_l, dst_c = np.mgrid[0:dest_nl:1, 0:dest_nc:1]
    #accouont for subset
    dst_l += src_geos_nav.start_line
    dst_c += src_geos_nav.start_col

    # convert the indices to latlon
    lats, lons = lc2ll(l=dst_l, c=dst_c, geos_obj=dst_geos_nav)

    #convert not latlon to srource indices defined by source lon_0
    src_l, src_c = ll2lc(lat=lats, lon=lons, geos_obj=src_geos_nav)


    #account for subset
    src_l -= src_geos_nav.start_line
    src_c -= src_geos_nav.start_col

    # mask to be sure
    m = (src_l < src_nl) & (src_l >= 0) & (src_c < src_nc) & (src_c >= 0)
    #reproject
    dest_geos_data[m] = src_geos_array[src_l[m], src_c[m]]

    if return_indices:
        return dest_geos_data, m, src_l, src_c
    return dest_geos_data


def reproject2ll(geos_array=None, bbox=None, geos_obj=None, **kwargs):
    """

    :param geos_array: an input 2d array of data in a geostationary projection
    :param bbox: a bounding box object (from generalutils.latlon) defining the target extent and resolution
    :param geos_obj: an instance of GeosNav class containing all required satellite scene parameters
    :return:
    """
    # lats and lons can be 2D or 1D.
    # we prefer 1D = use broadcasting because of these results of ll2rad
    # 2D arrays 'll2rad' avg run time was 0.501402807 sec in 10 runs
    # 1D 'll2rad' avg run time was 0.135642862 sec in 10 runs
    # no questions about this!

    lat = bbox.latitudes()
    lon = bbox.longitudes()

    return reproject_geos2ll(geos_array=geos_array, lat=lat, lon=lon, geos_obj=geos_obj, **kwargs)


def reproject_geos2ll(geos_array=None, lat=None, lon=None, geos_obj=None, lon_0=None, sweep=None, r_eq=None, r_pol=None,
                      sat_height=None, res_mrad=None, fd_lines=None, fd_fov=None):

    """
    Given a 2D array of geostationary satellite data data for a specific channel reproject (nearest neighbour)
    this data to latlon using the coordinates supplied in latitude, longitude arrays

    :param geos_array: an input 2d array of data in a geostationary projection
    :param lat: 1D or 2D array of latitudes
    :param lon: 1D or 2D array of longitudes
    :param geos_obj: an instance of GeosNav class containing all required satellite scene parameters
    :param lon_0: number the instrument SSP longitude
    :param sweep: string, y or x
    :param r_eq: number, equatorial radius
    :param r_pol: number, polar radius
    :param sat_height: number, satellite height above ground in meters
    :param res_mrad: number, float, the resolution of the image in microradians
    :param fd_lines: number int,
    :param fd_fov: number, the full disk field of view in radians
    :return:
    """

    nl, nc = geos_array.shape
    # assert nl == nc, 'The supplied array does not have an equal number of lines %s and columns %s' % (nl, nc)

    if geos_obj is None:
        l, c = ll2lc(lat=lat, lon=lon, lon_0=lon_0, sweep=sweep, r_eq=r_eq, r_pol=r_pol, sat_height=sat_height,
                     res_mrad=res_mrad, fd_lines=fd_lines, fd_fov=fd_fov)
    else:
        l, c = ll2lc(lat=lat, lon=lon, geos_obj=geos_obj)

    return geos_array[l, c]


class GeosNav:

    def __init__(self, lon_0=None, sweep=None, r_eq=None, r_pol=None, inv_flatten=None, ellipsoid=None, sat_height=None,
                 proj4_str=None, res_mrad=None, res_meter=None, fd_lines=None, fd_cols=None, fd_fov=None,
                 fd_geotransform=None, loff=None, coff=None, lfac=None, cfac=None, subset=None,
                 line_adj=0, col_adj=0):

        """
        :param lon_0: number - Longitude of the sub-satellite point.
        :param sweep: string - Sweep angle axis. Accepted values: 'y' (for most of the satellites) or 'x' for GOES 16 and 17.
        :param r_eq: number -  Ellipsoid equatorial radius in metres (a.k.a. semi major axis).
        :param r_pol: number - Ellipsoid polar radius in metres (a.k.a. semi minor axis).
        :param inv_flatten: float - Ellipsoid inverse flattening.
        :param ellipsoid: string - Ellipsoid name. Typical values are 'GRS80' or 'WGS84'.
        :param sat_height: - number - Satellite altitude above sub-satellite point in metres.
        :param proj4_str: string - Definition of the satellite's geostationary projection in proj4 format.
        :param res_mrad: float - Instrument angular resolution in micro-radians (μrad).
        :param res_meter: number - Instrument resolution in metres (pixel size).
        :param fd_lines: int - Number of lines in a full disk image.
        :param fd_cols: int - Number of lines in a full disk image.
        :param fd_fov: float - Field of view angle in radians (angular size of a full disk image).
        :param fd_geotransform: iterable (top left x, w-e meter resolution, 0, top left y, 0, n-s meter resolution)
        :param loff: number - Navigation coefficient 'line offset' (as used in HRIT format navigation).
        :param coff: number - Navigation coefficient 'column offset' (as used in HRIT format navigation).
        :param lfac: number - Navigation coefficient 'line factor' (as used in HRIT format navigation).
        :param cfac: number - Navigation coefficient 'column factor' (as used in HRIT format navigation).
        :param subset: iterable of coordinates (meters) defining an area of the full disk where the navigation applies to.
        :param: line_adj, number that is added to line indices computed inside ll2lc function
        :param: col_adj, number that is added to column indices computed inside ll2lc function

        NB: the line_adj and col_adj are so far used only my MSG HRIT reader (search for MSG .5 pixel shift story)

        """


        # assertions regarding satellite position and ellipsoid
        if proj4_str is None:
            assert lon_0 is not None, 'Satellite longitude is unknown. Specify "lon_0" parameter or a valid "proj4_str".'
            assert sat_height is not None, 'Satellite height is unknown. Specify "sat_height" parameter or a valid "proj4_str".'
            assert sweep is not None, 'Sweep angle axis unknown. Specify "sweep" parameter or a valid "proj4_str".'
            assert ((r_eq is not None and (r_pol is not None or inv_flatten is not None)) or ellipsoid is not None), \
                'Ellipsoid is unknown. Specify r_eq and r_pol (or inv_flatten), ' \
                'or ellipsoid parameters, or a valid "proj4_str".'
        else:
            assert sweep is not None or '+sweep=' in proj4_str, \
                'Sweep angle axis unknown. Specify "sweep" parameter or a valid "proj4_str".'
            assert sat_height is not None or '+h=' in proj4_str, \
                'Satellite height is unknown. Specify "sat_height" parameter or a valid "proj4_str".'
            assert lon_0 is not None or '+lon_0=' in proj4_str, \
                'Satellite longitude is unknown. Specify "lon_0" parameter or a valid "proj4_str".'
            assert ((r_eq is not None and (r_pol is not None or inv_flatten is not None)) or ellipsoid is not None) \
                or (('+a=' in proj4_str and ('+b=' in proj4_str or '+rf=' in proj4_str or '+f=' in proj4_str)) or
                    '+ellps=' in proj4_str), 'Ellipsoid is unknown. Specify r_eq and r_pol (or inv_flatten), ' \
                                             'or ellipsoid parameters, or a valid "proj4_str".'

        # assertions regarding satellite scene geometry
        assert all([loff is not None, coff is not None, lfac is not None, cfac is not None]) or \
            [(res_mrad is not None or res_meter is not None), fd_lines is not None, fd_fov is not None].count(True) >= 2 \
            or fd_geotransform is not None, \
            'Specify any two out of three instrument constants:\n \
            1. resolution (res_mrad or res_meter), 2. fulldisk size (fd_lines), 3. field of view angle (fd_fov).\n \
            Otherwise, fd_geotransform or four navigation coefficients need to be specified \n \
            line offset (loff), column offset (coff), line factor (lfac) and column factor (cfac).'

        if fd_fov is not None:
            assert fd_fov < 3.15, ' Expected units for the field of view angle (fd_fov) are radians.'

        self.lon_0 = lon_0
        self.sweep = sweep
        self.r_eq = r_eq
        self.r_pol = r_pol
        self.inv_flatten = inv_flatten
        self.ellipsoid = ellipsoid
        self.sat_height = sat_height
        self.proj4_str = proj4_str
        self.res_mrad = res_mrad
        self.res_meter = res_meter
        self.fd_fov = fd_fov
        self.fd_geotransform = fd_geotransform
        self.fd_lines = fd_lines
        self.fd_cols = fd_cols
        self.loff = loff
        self.coff = coff
        self.lfac = lfac
        self.cfac = cfac
        self.line_adj = line_adj
        self.col_adj = col_adj

        if self.proj4_str is not None:
            # parse the proj4 string into dict
            proj4_dict = dict(p.split('=') for p in self.proj4_str.replace('+no_defs', '').split())
            assert proj4_dict['+proj'] == 'geos', 'A valid geostationary projection string is required.'

            # solve sat_height in presence of a proj4 string
            if self.sat_height is None:
                self.sat_height = float(proj4_dict['+h'])
            elif '+h' in proj4_dict.keys():
                assert abs(self.sat_height - float(proj4_dict['+h'])) < 1, 'The sat_height and proj4 +h parameters conflict.'

            # solve lon_0 in presence of a proj4 string
            if self.lon_0 is None:
                self.lon_0 = float(proj4_dict['+lon_0'])
            elif '+lon_0' in proj4_dict.keys():
                assert abs(self.lon_0 - float(proj4_dict['+lon_0'])) < 0.001, 'The lon_0 and proj4 +lon_0 parameters conflict.'

            # solve ellipsoid in presence of a proj4 string
            if self.ellipsoid is None:
                if '+ellps' in proj4_dict.keys():
                    self.ellipsoid = proj4_dict['+ellps']
            elif '+ellps' in proj4_dict.keys():
                assert self.ellipsoid == proj4_dict['+ellps'], 'The ellipsoid and proj4 +ellps parameters conflict.'

            # solve semi major axis in presence of a proj4 string
            if self.r_eq is None:
                if '+a' in proj4_dict.keys():
                    self.r_eq = float(proj4_dict['+a'])
                elif '+ellps' in proj4_dict.keys():
                    self.r_eq = pyproj.pj_ellps[proj4_dict['+ellps']]['a']
            elif '+a' in proj4_dict.keys():
                assert abs(self.r_eq - float(proj4_dict['+a'])) < 1, 'The r_eq and proj4 +a parameters conflict.'

            # solve semi minor axis in presence of a proj4 string
            if self.r_pol is None:
                if '+b' in proj4_dict.keys():
                    self.r_pol = float(proj4_dict['+b'])
                elif '+rf' in proj4_dict.keys() and self.r_eq:
                    self.r_pol = (self.r_eq * 1 / float(proj4_dict['+rf']) - self.r_eq) * -1
                elif '+f' in proj4_dict.keys() and self.r_eq:
                    self.r_pol = (self.r_eq * float(proj4_dict['+f']) - self.r_eq) * -1
                elif '+ellps' in proj4_dict.keys() and self.r_eq:
                    self.r_pol = (self.r_eq * 1 / pyproj.pj_ellps[proj4_dict['+ellps']]['rf'] - self.r_eq) * -1
            elif '+b' in proj4_dict.keys():
                assert abs(self.r_pol - float(proj4_dict['+b'])) < 1, 'The r_pol and proj4 +b parameters differ.'

            # solve inverse flattening in presence of a proj4 string
            if self.inv_flatten is None:
                if '+rf' in proj4_dict.keys():
                    self.inv_flatten = float(proj4_dict['+rf'])
                elif '+f' in proj4_dict.keys():
                    self.inv_flatten = 1 / float(proj4_dict['+f'])
                elif '+ellps' in proj4_dict.keys():
                    self.inv_flatten = pyproj.pj_ellps[proj4_dict['+ellps']]['rf']
            elif '+rf' in proj4_dict.keys():
                assert abs(self.inv_flatten - float(proj4_dict['+rf'])) < 0.001, 'The inv_flatten and proj4 +rf parameters conflict.'
            elif '+f' in proj4_dict.keys() and self.r_eq:
                assert abs(self.inv_flatten - 1 / float(proj4_dict['+f'])) < 0.001, 'The inv_flatten and proj4 +rf parameters conflict.'
            elif '+a' in proj4_dict.keys() and '+b' in proj4_dict.keys():
                i_f = float(proj4_dict['+a']) / (float(proj4_dict['+a']) - float(proj4_dict['+b']))
                assert abs(self.inv_flatten - i_f) < 0.001, 'The inv_flatten and proj4 +a, +b parameters conflict.'

            # solve sweep axis in presence of a proj4 string
            if self.sweep is None:
                self.sweep = proj4_dict['+sweep']

        else:
            # solve ellipsoid parameters in absence of a proj4 string
            if self.r_eq is None and self.ellipsoid is not None:
                self.r_eq = pyproj.pj_ellps[self.ellipsoid]['a']
            if self.r_pol is None and self.ellipsoid is not None:
                self.r_pol = (self.r_eq * 1 / pyproj.pj_ellps[self.ellipsoid]['rf'] - self.r_eq) * -1
            if self.inv_flatten is None and self.ellipsoid is not None:
                self.inv_flatten = pyproj.pj_ellps[self.ellipsoid]['rf']
            if self.r_pol is None and self.inv_flatten is not None:
                self.r_pol = (self.r_eq * 1 / self.inv_flatten - self.r_eq) * -1

            # solve proj4 string
            self.proj4_str = f'+proj=geos +a={self.r_eq:.2f} +b={self.r_pol:.2f} +h={self.sat_height:.2f} +lon_0={self.lon_0} +sweep={self.sweep} +no_defs'

        # some assertions regarding scene geometry parameters
        if self.res_meter is not None and self.res_mrad is not None:
            assert abs(self.res_meter - self.res_mrad * self.sat_height / 1e6) < 0.1, 'The res_meter and res_mrad resolution differ.'
        if self.lfac is not None:
            assert abs(self.lfac) == abs(self.cfac), 'The lfac and cfac navigation coefficients are expected to be equal.'
            # the following condition is not the case e.g. for half-disk images. However, also partial images will be represented as fulldisk by convention.
            assert abs(self.loff) == abs(self.coff), 'The loff and coff navigation coefficients are expected to be equal.'
        if self.fd_geotransform is not None:
            assert abs(self.fd_geotransform[0]) == abs(self.fd_geotransform[3]) and abs(self.fd_geotransform[1]) == abs(self.fd_geotransform[5]) \
                   and (self.fd_geotransform[2], self.fd_geotransform[4]) == (0, 0), 'The fd_geotransform does not represent a regular geostationary satellite scene.'

        # solve res_mrad
        if self.res_mrad is None:
            if self.lfac is not None:
                self.res_mrad = abs(-2 ** 16 / float(np.rad2deg(self.lfac))) * 1e6
            elif self.cfac is not None:
                self.res_mrad = abs(-2 ** 16 / float(np.rad2deg(self.cfac))) * 1e6
            elif self.fd_lines is not None and self.fd_fov is not None:
                self.res_mrad = self.fd_fov / self.fd_lines * 1e6
            elif self.fd_cols is not None and self.fd_fov is not None:
                self.res_mrad = self.fd_fov / self.fd_cols * 1e6
            elif self.res_meter is not None:
                self.res_mrad = abs(self.res_meter / self.sat_height) * 1e6
            elif self.fd_geotransform is not None:
                self.res_mrad = abs(self.fd_geotransform[5] / self.sat_height) * 1e6 # JANO isn't this going to make negatovge res_mrad allways?

        # solve fd_lines
        if self.fd_lines is None:
            if self.loff is not None:
                self.fd_lines = int(abs(self.loff * 2))
            elif self.res_mrad is not None and self.fd_fov is not None:
                self.fd_lines = int(round(self.fd_fov / self.res_mrad * 1e6))
            elif self.fd_geotransform is not None:
                self.fd_lines = int(abs(round(self.fd_geotransform[3] / self.fd_geotransform[5]) * 2))

        # solve fd_fov
        if self.fd_fov is None:
            if self.fd_lines is not None:
                self.fd_fov = self.fd_lines * self.res_mrad / 1e6
            elif self.fd_cols is not None:
                self.fd_fov = self.fd_cols * self.res_mrad / 1e6

        # solve ncols
        if self.fd_cols is None:
            if self.coff is not None:
                self.fd_cols = int(abs(self.coff * 2))
            elif self.res_mrad is not None and self.fd_fov is not None:
                self.fd_cols = int(round(self.fd_fov / self.res_mrad * 1e6))
            elif self.fd_geotransform is not None:
                self.fd_cols = int(abs(round(self.fd_geotransform[0] / self.fd_geotransform[1]) * 2))

        # solve res_meter
        if self.res_meter is None:
            self.res_mrad = abs(self.res_mrad)
            self.res_meter = abs(self.res_mrad * self.sat_height) / 1e6

        # solve lfac, cfac, loff, coff
        if self.lfac is None:

            self.lfac = int(np.round(np.deg2rad(2 ** 16 / (self.res_mrad*1e-6))))
            self.cfac = int(np.round(np.deg2rad(2 ** 16 / (self.res_mrad*1e-6))))
            self.loff = self.fd_lines // 2
            self.coff = self.fd_cols // 2

        # solve fd_geotransform and extent
        left = self.fd_lines * self.res_meter / -2
        top = self.fd_lines * self.res_meter / 2
        if self.fd_geotransform is None:
            self.fd_geotransform = left, self.res_meter, 0, top, 0, -self.res_meter
        self.fd_extent = (left, -left, -top, top)

        # add rounded resolution in km
        self.res_km = np.round(self.res_meter / 1e3, 1)

        if subset is not None:
            assert len(subset) == 4, f'Invalid subset={subset}. It needs to be an iterable of coordinates ' \
                                          f'representing a chunk of geostationary space (xmin, xmax, ymin, ymax) '
            assert subset[0] < subset[1] and subset[2] < subset[3], \
                                          f'Invalid subset={subset} has invalid values. It needs to be an iterable of coordinates ' \
                                          f'representing a chunk of geostationary space (xmin, xmax, ymin, ymax) '

            # expect allways  xmin, xmax, ymin, ymax
            xmin, xmax, ymin, ymax = subset
            self.start_line, self.start_col = YX2lc(Y=ymax, X=xmin, geos_obj=self)
            self.end_line, self.end_col = YX2lc(Y=ymin, X=xmax, geos_obj=self)
            self.nlines = self.end_line - self.start_line
            self.ncols = self.end_col - self.start_col
            self.geotransform = xmin, self.res_meter, 0, ymax, 0, -self.res_meter
            self.extent = subset
        else:
            self.start_line = self.start_col = 0
            self.end_line, self.end_col = self.fd_lines, self.fd_cols
            self.nlines, self.ncols = self.fd_lines, self.fd_cols
            self.geotransform = self.fd_geotransform
            self.extent = self.fd_extent
        self.image_extent = self.start_line, self.end_line, self.start_col, self.end_col


    def __str__(self):
        return f'Geostationary navigation object summary:' \
            f'\n\tproj4 string: {self.proj4_str}' \
            f'\n\tgeotransform: {self.geotransform}' \
            f'\n\textent in meters : {self.extent}' \
            f'\n\textent in lines/columns : {self.image_extent}' \
            f'\n\timage size : {self.nlines, self.ncols}' \
            f'\n\tnative resolution in microradians: {self.res_mrad} ~ {self.res_meter} meters' \
            f'\n\tfull disk field of view: {self.fd_fov} radians ~ {np.rad2deg(self.fd_fov)} degrees' \
            f'\n\tfull disk image size {self.fd_lines, self.fd_cols}' \
            f'\n\tscaling and offset coefficients: LFAC: {self.lfac} LOFF: {self.loff} CFAC: {self.cfac} COFF: {self.coff}'

    @property
    def is_subset(self):
        """
        Compute if a navigation was subsetted. Return true if it was else False
        :return:
        """
        return self.nlines != self.fd_lines or self.ncols != self.fd_cols

    def to_lon_0(self, src_geos_array=None, dst_geos_nav=None, dst_fill_value=None, run_in_parallel=False,
                 return_indices=False):
        """
        Reproject the src_geos_array associated with this navigation (self) to dst_geos_nav.

        :param src_geos_array: numpy 2D array representing source data
        :param dst_geos_nav: instance of the destination projection GesoNav class
        :param dst_fill_value: number, the fill value for destination or reprojected data
        :param run_in_parallel: bool, if True the reprojection is distributed  onto multiple cores
        :param return_indices: bool, if True the repr indices (lines and col) are also returend
        :return:
        """

        if not run_in_parallel:

            return reproject_geos2geos(
                                        src_geos_array=src_geos_array, src_geos_nav=self,
                                        dst_geos_nav=dst_geos_nav, dst_fill_value=dst_fill_value,
                                        return_indices=return_indices)
        else:
            return reproject_geos2geos_mp(src_geos_array=src_geos_array, src_geos_nav=self,
                                          dst_geos_nav=dst_geos_nav, dst_fill_value=dst_fill_value,
                                          return_indices=return_indices)




    @property
    def geometric_mask(self):
        """
        Compute a geometric mask of the earth inside  the native
        coordinate space of the image
        :return:
        """
        h = self.sat_height + self.r_eq
        aeq = 1 - self.r_eq ** 2 / (h ** 2)
        ap = 1 - self.r_pol ** 2 / (h ** 2)
        xmax = np.arccos(np.sqrt(aeq))
        ymax = np.arccos(np.sqrt(ap))

        # if self.sweep == 'y':
        #     xmax = np.arccos(np.sqrt(aeq))
        #     ymax = np.arccos(np.sqrt(ap))
        # else:
        #     xmax = np.arccos(np.sqrt(ap))
        #     y max = np.arccos(np.sqrt(aeq))

        yr, xr = self.get_radian_arrays()

        yr /= ymax
        xr /= xmax
        m = xr ** 2 + yr ** 2 <= 1.

        return m

    def get_radian_arrays(self, one_dim=False):
        """
        Creates radians of array coordinates for this navigation object
        :param one_dim, bool, defaults to False, if True the reaturned array are one dimensional
        :return: a tuple of 2D array representing the native radian coordinates of this navigation
        """
        xres_rad = self.res_mrad * 1e-6
        yres_rad = -self.res_mrad * 1e-6
        yoffset = self.fd_lines // 2
        xoffset = self.fd_cols // 2
        y_rad = (np.arange(self.start_line, self.end_line) - yoffset) * yres_rad + yres_rad/2
        x_rad = (np.arange(self.start_col, self.end_col) - xoffset) * xres_rad + xres_rad/2

        if one_dim: return y_rad, x_rad
        y_rad = np.lib.stride_tricks.as_strided(y_rad, shape=(self.nlines, self.ncols),
                                                strides=(y_rad.itemsize, 0))
        x_rad = np.lib.stride_tricks.as_strided(x_rad, shape=(self.nlines, self.ncols),
                                                strides=(0, x_rad.itemsize,))

        return y_rad, x_rad

    def get_latlon_1d(self, mask=None):
        """
        Calulates geographic coordinates for each non-masked pixel in the geostationary image.
        Returns the result as two 1d arrays (latitudes, longitudes)

        :param mask: Allows user to specify desired mask of the full disk image. Geometric mask is used if not specified.
        :return: latitudes, longitudes (numpy 1d arrays)
        """

        line_indices = np.arange(self.start_line, self.end_line + 1)
        column_indices = np.arange(self.start_col, self.end_col + 1)

        lines2d = np.lib.stride_tricks.as_strided(line_indices, shape=(self.nlines, self.ncols), strides=(line_indices.itemsize, 0))  # col duplicate
        columns2d = np.lib.stride_tricks.as_strided(column_indices, shape=(self.nlines, self.ncols), strides=(0, column_indices.itemsize))  # row duplicate

        if mask is None:
            mask = ~self.geometric_mask
        lines2d_masked = np.ma.masked_array(lines2d, mask=mask)
        columns2d_masked = np.ma.masked_array(columns2d, mask=mask)

        lines1d = lines2d_masked.compressed()
        columns1d = columns2d_masked.compressed()

        return lc2ll(l=lines1d, c=columns1d, geos_obj=self)
