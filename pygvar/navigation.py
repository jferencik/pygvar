import struct
import datetime
import math
from pygvar.utils import f2int
import numpy as np


"""
This module is an updated mixed optimised port of ancillary ELUG Fortran code that performs GOES GVAR navigation.
It does not handle sounder  instrument but with some changes (constants and conditions for imc) it could
as the instruments are similar

Note there are 3 types of coordinates when it comes to GVAR files

   1. instrument coordinates, points on earth thta corresponds to pixels in camera

    their range is instrument given 1-30680 in x direction and 1-15780 in y direction

  2 relative or file coordinates. As the instruments scan always a part of the whole pixel space every image
  start at a certain line/element in instrument coordinate system. The file coordinates start at  0,0 to nl/nc-1
  Their reslution is identical, the fiel coordinates are just offsetted instrument coordinates
  3. latlon coordinates. These are defined only inside the earth disk. this is why some line/pixels do not have lat/lon
  
The navigation functions use the relative coordinates as input, converts these coordinates to instrument, perform
transformations and then return latlon or convert the instrument coordinates back to file coordinates


Reference
https://goes.gsfc.nasa.gov/text/ELUG0398.pdf

"""



def datetime2mcidas(input_datetime=None):


    tt = input_datetime.timetuple()
    yyyddd =  (input_datetime.year - 1900) * 1000 + tt.tm_yday

    HHMMSSmmm = int(input_datetime.hour*1e4 + input_datetime.minute * 1e2 + input_datetime.second)*1000 + input_datetime.microsecond
    return yyyddd,HHMMSSmmm
'''
def satangle2mcidas(angle=None):
    """
    Converts an attitude or misalignment angle to a form suitable to be injected into GVAR navigation ans required by McIDAS

    An satellite navigation or attitude angle is structure of 55 words. Here is an example of an attitude angle as modeled by C++ class or C struct


    SelFloat is a specific floating point format. Some of GVAR data come in this format. My take is that it is used by the processor abpard thye Imager instrument or
    it is used ofr unknown historical reasons.
    Nevertheless i hvae developed functions to decode it by mainly copying C software like gvartool and McIDAS java code


    AttitudeAngle:
        SelFloat         Exponential_magnitude; /* 523  526*/
        SelFloat         Exponential_time_constant; /* 527  530*/
        SelFloat         Constant_mean_attitude_angle; /* 531  534*/
        uint32           Number_of_sinusoidals_per_angle; /* 535  538*/
        Polar            Sinusoid[15];		/* 539  658*/
        uint32           Number_of_monomial_sinusoids; /* 659  662*/
        MonomialSinusoid Monomial[4]; /* 663  742*/

    Polar:
        SelFloat Mg;			/* 539  542*/
        SelFloat Ph ;			/* 543  546*/
    MonomialSinusoid::
        uint32   Order_of_applicable;	/* 663  666*/
        uint32   Order_of;		/* 667  670 */
        Polar    Sinusoid;		/* 671  678 */
        SelFloat Angle_of_epoch_where_monomial_is_zero; /* 679  682*/




    :param angle:
    :return: an iterable fo same length where each element is altered (usually multiplied ) by a specific number

    """
    assert angle is not None, '"angle" can not be None'
    assert isiter(angle), '"angle" is not iterable'
    assert len(angle) == 55, '"angle" has a wrong length %d . It should be 55' % len(angle)
    #first four elements

    out_angle = [None] * 55
    out_angle[0] = f2int(angle[0]* 1e7 )
    out_angle[1] = f2int(angle[1] * 1e2)
    out_angle[2] = f2int(angle[2] * 1e7)
    out_angle[3] = f2int(angle[3])

    # then fifteen sinusoidals
    for i in range(15):
        j = 4+2*i
        k = j+1
        out_angle[j] = f2int(angle[j]*1e7)
        out_angle[k] = f2int(angle[k]*1e7)

    out_angle[33] = f2int(angle[33])

    #four Monomials, firts two elements are int
    for i_ in range(34, 55, 5):
            monomial_chunk = angle[i_:i_+5]
            out_angle[i_:i_+2] = map(f2int,monomial_chunk[:2])
            out_angle[i_+2: i_+5] = map(f2int, [k_*1e7 for k_ in monomial_chunk[2:]])
    return  out_angle
'''
def prepare_mcidas_nav(b0):
    """
    This function builds a MCIDAS navigation array from a pygvar block0 class
    :param b0:
    :return: list of ints
    """
    #nav = np.zeros((640,), dtype='i4')
    nav = []
    #W1 navigation type; 4 bytes ASCII, GVAR value
    nav.append(struct.unpack('=i', 'GVAR')[0])
    #W2 ASCII string, usually a letter followed by three integers

    nav.append(struct.unpack('=i', 'E000')[0])
    #nav.append(struct.unpack('=i', '000R')[0])
    #nav.append(struct.unpack('=i', '230J')[0])



    # all imc is questionable
    # W3 IMC imager scan status; bits 0-15 are right justified, with 15 the least significant; IMC active flag is bit 8, counting from the least significant bit; 1=active; see OGE Table 3-6, bytes 3-6
    #ss0, ss1, ss2, ss3 = b0.words_value(start_word=3, end_word=6, byte_order=b0.BYTE_ORDER, word_size='4B')
    ss0, ss1, ss2, ss3 = struct.unpack('<4B', b0.iscan.buffer)
    #is0, is1 = struct.unpack('<2H', b0.iscan.buffer)

    #print int(va, 2)
    #print ss0, ss1, ss2, ss3

    #nav.append(int(b0_data.imc))
    #imager scan status; bits 0-15 are right justified, with 15 the least significant; IMC active flag is bit 8, counting from the least significant bit; 1=active; see OGE Table 3-6, bytes 3-6
    nav.append(int(ss1))


    # W4     imager scan status; bits 16-31 are right justified, with 31 the least significant; yaw-flip processing enabled flag is bit 16, counting from the least significant bit; 1=enabled; see OGE Table 3-6, bytes 3-6
    # nav.append(int(b0_data.yaw_flip))
    nav.append(int(ss3))



    #W5 reserved
    nav.append(0)
    # W6 reference longitude, rad*10000000
    nav.append(f2int(b0.reflo.python_value() * 1e7))
    # W7 Ref distance from nominal (km *10^7)
    nav.append(f2int(b0.refra.python_value() * 1e7))
    #W8 Ref Latitude  (rad * 10^7 )
    nav.append(f2int(b0.refla.python_value() * 1e7))
    #W9 Ref Yaw       (rad * 10^7 )
    nav.append(f2int(b0.reforyaw.python_value() * 1e7))
    #W10 Ref  attitude roll (rad * 10^7 )
    nav.append(f2int(b0.refatroll.python_value() * 1e7))
    # W11 Ref attitude pitch (rad * 10^7 )
    nav.append(f2int(b0.refatpitch.python_value() * 1e7))
    # # W12  Ref attitude yaw (rad * 10^7 )
    nav.append(f2int(b0.refatyaw.python_value() * 1e7))
    # #W13-14 Epoch time data BCD format
    #bcd_bytes =  b0.words_value(start_word=323, end_word=330, byte_order=b0.BYTE_ORDER, word_size='2i')
    bcd_bytes =  struct.unpack('>2i',b0.epochdate.buffer)

    nav.append(bcd_bytes[0])
    nav.append(bcd_bytes[1])
    #W15 delta from epoch time, minutes*100
    nav.append(f2int(b0.imc_enable_from_epoch.python_value() * 1e2))
    #W16 image motion compensation roll, rad*10000000
    nav.append(f2int(b0.compens_roll.python_value() * 1e7))
    #W17 image motion compensation pitch, rad*10000000
    nav.append(f2int(b0.compens_pitch.python_value() * 1e7))
    #W18 image motion compensation yaw, rad*10000000
    nav.append(f2int(b0.compens_yaw.python_value() * 1e7))
    #W19-31 longitude delta from reference values, rad*10000000
    #[nav.append(f2int(b0.change_longitude[i] * 1e7)) for i in range(13)]
    nav+=[f2int(e.python_value() * 1e7) for e  in b0.change_lon]

    #W32-42 radial distance delta from reference values, km*10000000
    #[nav.append(f2int(b0_data.change_radial_dist[i])  * 1e7) for i in range(11)]
    nav+=[f2int(e.python_value() * 1e7) for e in b0.change_radist]
    # W43-51 sine of the geocentric latitude delta values, units*10000000
    nav+=[f2int(e.python_value() * 1e7) for e in b0.sin_geoc_lat]
    # W52-60 sine of the orbit yaw delta values, units*10000000
    nav+=[f2int(e.python_value() * 1e7) for e in b0.sin_orbit_yaw]
    #w61 daily solar rate, rad/min*10000000
    nav.append(f2int(b0.daily_solar_rate.python_value() * 1e7))
    #W62 exponential start time from epoch, minutes*100
    nav.append(f2int(b0.exp_start_time_from_epoch.python_value() * 1e2))
    #attitude and missalignment angles follow here
    #W63-117 roll attitude angle (OGE Table 3-6, bytes 523-742)
    #nav+=satangle2mcidas(b0.roll_angle)
    nav += b0.roll_angle.to_mcidas()

    #W118-127 reserved
    [nav.append(0) for i in range(10)]

    #W128 MORE
    nav.append(struct.unpack('=i', 'MORE')[0])
    #W129 GVAR
    nav.append(struct.unpack('=i', 'GVAR')[0])
    #W 130 -184 pitch attitude angle
    nav += b0.pitch_angle.to_mcidas()

    #W 185-239 yaw attitude angle
    nav += b0.yaw_angle.to_mcidas()
    # #W 240-255 reserved
    [nav.append(0) for i in range(16)]

    # W256 MORE
    nav.append(struct.unpack('=i', 'MORE')[0])
    # W257 GVAR
    nav.append(struct.unpack('=i', 'GVAR')[0])
    #
    #W 258-312 roll misalignment angle
    nav+= b0.roll_missalignment_angle.to_mcidas()
    #W 313-367 pitch misalignment angle
    nav += b0.pitch_missalignment_angle.to_mcidas()
    yyyddd, HHMMSSmmm = datetime2mcidas(b0.tched.python_value()) # current time
    # W368 year and Julian day, yyyddd
    nav.append(yyyddd)
    #W369 HHMMSSmmm
    nav.append(HHMMSSmmm)

    #W 370 imager/sounder instrument flag; 1=imager, 2=sounder
    nav.append(1)

    #W 371-379 reserved
    [nav.append(0) for i in range(9)]

    # W380 instrument nadir, north/south cycles; see OGE Table 3-6, byte 6305
    nav.append(b0.ns_cycles)
    #W 381 instrument nadir, east/west cycles; see OGE Table 3-6, byte 6306

    nav.append(b0.ew_cycles)

    #W 382 instrument nadir, north/south increments; see OGE Table 3-6, byte 6307-6308

    nav.append(b0.ns_incr)

    #W 383 instrument nadir, east/west increments; see OGE Table 3-6, byte 6309-6310
    nav.append(b0.ew_incr)

    # W384 MORE
    nav.append(struct.unpack('=i', 'MORE')[0])
    # W385 GVAR
    nav.append(struct.unpack('=i', 'GVAR')[0])
    #W 386 511 reserved
    [nav.append(0) for i in range(126)]
    # W512 MORE
    nav.append(struct.unpack('=i', 'MORE')[0])
    # 513 GVAR
    nav.append(struct.unpack('=i', 'GVAR')[0])

    #W 514-640 reserved
    [nav.append(0) for i in range(127)]

    return nav

def compute_attitude_angle(att_angle=None,solar_orbit_angle=None, exp_time_delay=None):
    """
    
    attitude fields
     _fields_ = (
        ('exp_mag', GouldFloat),
        ('exp_time_const', GouldFloat),
        ('const_mean_attit_angle', GouldFloat),
        ('nsinusoidals', C.c_uint32),
        ('sinusoids', Polar*15),
        ('n_monomial_sinusoids', C.c_uint32),
        ('monomials', MonomialSinusoid*4),
    )
    polar fields 
    _fields_ = (
        ('mg', GouldFloat),
        ('ph', GouldFloat),
    )
    monomial sinusoid fields
    _fields_ = (
        ('order_of_applicable', C.c_uint32),
        ('order_of', C.c_uint32),
        ('sinusoid', Polar),
        ('angle_of_epoch_where_monomial_is_zero', GouldFloat),
    )
    
    Attitude = const + exponential + solart orbit compensated sinusoids + monomials
    
    :param att_angle: 
    :param solar_orbit_angle: 
    :param exp_time_delay: 
    :return: 
    """
    #constant
    att = att_angle.const_mean_attit_angle.python_value()
    #exponential
    if exp_time_delay>=0:
        att += att_angle.exp_mag.python_value() * math.exp(-exp_time_delay/att_angle.exp_time_const.python_value())

    # calculation of sinusoids

    for j in range(att_angle.nsinusoidals):

        att += att_angle.sinusoids[j].mg.python_value() * math.cos(solar_orbit_angle*j+att_angle.sinusoids[j].ph.python_value())

    for l in range(att_angle.n_monomial_sinusoids):
        monomial = att_angle.monomials[l]
        sinusoid = monomial.sinusoid
        # order of sinusoid
        order_of_sinusoid = float(monomial.order_of_applicable)
        #order of monomial
        order_of_monomial = float(monomial.order_of)
        att += sinusoid.mg.python_value() * ((solar_orbit_angle - monomial.angle_of_epoch_where_monomial_is_zero)**order_of_monomial) * math.cos(order_of_sinusoid*solar_orbit_angle+sinusoid.ph.python_value())

    return att

def ll2le(lats=None, lons=None, b0=None, lineres=1., colres=1.):
    """
    Converts  geographic|latlon coordinates  to file or relative coordinates line/element|pixel
    :param lats: numpy array of lats, float, 
    :param lons: numpy array of lons, float
    :param b0: a valid block0 from a GVAR file
    :param lineres: 1
    :param colres: 1
    :return: numpy array (lines, columns) corresponding to input latlon coordinates
    
    
    
    """


    # print '#'*50,
    # print 'GVAR',
    # print '#' * 50

    #init
    #TODO FIGURE OUT the resolution deal inside an area file
    resLine = lineres #setup by user, default value if 1 (vis chn res is used)
    resElement = colres #setup by user, default value if 1 (vis chn res is used)
    magLine = 1.
    magElement = 1.
    startLine = 0
    startElement = 0
    startImageLine = b0.infln #init from b0
    startImageElement = b0.iwfpx #init from b0
    isLineFlipped=False
    lineOffset = 0
    #print resLine, resElement, magLine, magElement, startLine,startElement, startImageLine, startImageElement, isLineFlipped, lineOffset
    PI = math.pi
    #DEG = 180. / PI
    RAD = PI / 180.  # degrees to radians conversion pi/180

    #instr = 1  # 1 imager 2 sounder...we do not handle sounder...goodbye
    isEastPositive = True  # this refers to longitude being positive. MCIDAS can mess with logitudes specially when reprojecting
    # earth

    nomorb = 42164.365  # nominal radial distance of instrument (km)
    aec = 6378.137  # earth equatorial radius (km)
    ferc = 1. - (6356.7533 / aec)  # earth flattening coeff = 1-(be/ae)
    aebe2c = 1. / (1 - ferc) ** 2
    aebe3c = aebe2c - 1.
    aebe4c = ((1. - ferc) ** 4.) - 1.



    # instrument(IMAGER) specific data these data is only for imager

    incmax = 6136  # number of increments per cycle

    elvmax = 0.220896  # bounds in elevation (radians)
    scnmax = 0.24544  # bounds in scan angle (radians)
    elvinc = 8.e-6  # change in elevation angle per increment (rad)
    scninc = 16.e-6  # change in scan angle per increment (radians)
    elvln = 28.e-6  # elevation angle per detector line (radians)
    scnpx = 16.e-6  # scan angle per pixel (radians)
    nsnom = .220896  # north-south center of instrument (=4.5 x incmax x elvinc)

    ewnom = .24544 # east-west center of instrument (=2.5 x incmax x sncinc)

    nadnsc = b0.ns_cycles
    nadnsi = b0.ns_incr
    nadewc = b0.ew_cycles
    nadewi = b0.ew_incr

    if nadnsc != 0 and nadnsi != 0 and nadewc != 0 and nadewi != 0:
        # only imager is handled , sounder is ignored
        elvmax = (nadnsc * incmax + nadnsi) * elvinc
        # update the scan angle in x (first?))
        scnmax = (nadewc * incmax + nadewi) * scninc

    #print nadnsc, nadnsi,nadewc,nadewi, elvmax, scnmax
    # time computation
    jan_1_1950 = datetime.datetime(1950, 1, 1)
    bcdtime = b0.epochdate.python_value()


    # this is time in milisec from epoch jan 1 1950
    diff = bcdtime - jan_1_1950
    epoch = diff.total_seconds() / 60.
    # currect header block
    #firts convert tched into yyyydoy and hhmmssmmm
    ydoy, hhmmssmmm = datetime2mcidas(b0.tched.python_value())
    #convert hhmmssmmm into a float, byte order is hardcoded here for little endian
    hhmmssmmm_f = struct.unpack('{byte_order:s}{size:d}{format_str:s}'.format(byte_order='<', size=1, format_str='f'),struct.pack('{byte_order:s}{size:d}{format_str:s}'.format(byte_order='<', size=1,format_str='i'), hhmmssmmm))[0]
    year = 1900 + ydoy / 1000
    day = ydoy - ydoy / 1000 * 1000
    hour = int(hhmmssmmm_f / 10000)
    minuts = int(hhmmssmmm_f / 100 - hour * 100)
    secs = hhmmssmmm_f - float(100 * minuts) - float(10000 * hour)

    j = day + 1461 * (year + 4799) / 4 - 3 * ((year + 4899) / 100) / 4 - 2465022

    imgtim = float(j * 1440. + hour * 60. + minuts + (secs / 60.))

    #print b0.iscan.imc_active, b0.iscan.yaw_flip
    #if imc is 1 it is not active
    imc = 0 if  b0.iscan.imc_active != 0  else  1
    #imc = b0.iscan.imc_active  # should suffice because it is either 0 or 1
    iflip = -1 if b0.iscan.yaw_flip else 1

    # imc = b0.iscan.imc_active
    # iflip = b0.iscan.yaw_flip
    #print imc, iflip



    # assign reference values to the subsatellite longitude and
    # latitude, the radial distance and the orbit yaw.
    lam = b0.reflo.python_value()
    dr = b0.refra.python_value()

    phi = b0.refla.python_value()

    psi = b0.reforyaw.python_value()
    subpoint = []
    subpoint.append( phi / RAD)
    subpoint.append(lam / RAD)

    # assign reference values to the attitudes and misalignments
    roll = b0.refatroll.python_value()
    pitch = b0.refatpitch.python_value()
    yaw = b0.refatyaw.python_value()


    rma = 0.
    pma = 0.
    if imc != 0:
        # set reference radial distance, latitude and orbit yaw to zero, does it make sens?
        dr = 0.
        phi = 0.
        psi = 0.
        # compute time since epoch (in minutes)
        ts = imgtim - epoch
        # computes orbit angle and the related trigonometric functions.
        # earth rotational rate=.729115e-4 (RAD/s)
        w = 0.729115e-4 * 60.0 * ts
        sw = math.sin(w)
        cw = math.cos(w)
        sw1 = math.sin(0.927 * w)
        cw1 = math.cos(0.927 * w)
        s2w = math.sin(2. * w)
        c2w = math.cos(2. * w)
        sw3 = math.sin(1.9268 * w)
        cw3 = math.cos(1.9268 * w)
        # computes change in the imc longitude from the reference
        change_lon = [e.python_value() for e in b0.change_lon]
        lam = lam + change_lon[0] + (change_lon[1] + change_lon[2] * w) * w \
              + (change_lon[9] * sw1 + change_lon[10] * cw1 + change_lon[3] * sw \
                 + change_lon[4] * cw + change_lon[5] * s2w + change_lon[6] * c2w \
                 + change_lon[7] * sw3 + change_lon[0] * cw3 + w * (change_lon[11] * sw + change_lon[12] * cw)) * 2.

        change_radist = [e.python_value() for e in b0.change_radist]
        # computes change in radial distance from the reference (km)
        dr = dr + change_radist[0] + change_radist[1] * cw + change_radist[2] * sw \
             + change_radist[3] * c2w + change_radist[4] * s2w + change_radist[5] * cw3 \
             + change_radist[6] * sw3 + change_radist[7] * cw1 \
             + change_radist[8] * sw1 + w * (change_radist[9] * cw + change_radist[10] * sw)

        # computes the sine of the change in the geocentric latitude
        change_sin_geoc_lat = [e.python_value() for e in b0.sin_geoc_lat]
        dlat = change_sin_geoc_lat[0] + change_sin_geoc_lat[1] * cw + change_sin_geoc_lat[2] * sw \
               + change_sin_geoc_lat[3] * c2w + change_sin_geoc_lat[4] * s2w + \
               w * (change_sin_geoc_lat[5] * cw + change_sin_geoc_lat[6] * sw) + change_sin_geoc_lat[7] * cw1 + \
               change_sin_geoc_lat[8] * sw1

        # computes geocentric latitude by using an expansion for arcsine
        phi = phi + dlat * (1. + dlat * dlat / 6.)

        change_sin_orbit_yaw = [e.python_value() for e in b0.sin_orbit_yaw]
        # computes sine of the change in the orbit yaw
        dyaw = change_sin_orbit_yaw[0] + change_sin_orbit_yaw[1] * sw + change_sin_orbit_yaw[2] * cw \
               + change_sin_orbit_yaw[3] * s2w + change_sin_orbit_yaw[4] * c2w \
               + w * (change_sin_orbit_yaw[5] * sw + change_sin_orbit_yaw[6] * cw) \
               + change_sin_orbit_yaw[7] * sw1 + change_sin_orbit_yaw[8] * cw1

        # computes the orbit yaw by using an expansion for arcsine.
        psi = psi + dyaw * (1. + dyaw * dyaw / 6.)

        # calculation of changes in the instrument orbit ends here

    # conversion of the imc longitude and orbit yaw to the subinstrument
    # longitude and the orbit inclination (ref: goes-pcc-tm-2473, inputs
    # required for earth location and gridding by sps, june 6, 1988)
    slat = math.sin(phi)
    syaw = math.sin(psi)
    sinoi = slat * slat + syaw * syaw
    cosoi = math.sqrt(1. - sinoi)
    sinoi = math.sqrt(sinoi)

    if slat == 0.0 and syaw == 0.0:
        u = 0.0
    else:
        u = math.atan2(slat, syaw)
    sinu = math.sin(u)
    cosu = math.cos(u)

    # computes longitude of the ascending node
    asc = lam - u
    sinasc = math.sin(asc)
    cosasc = math.cos(asc)

    # computes the subinstrument geographic latitude
    sublat = math.atan(aebe2c * math.tan(phi))

    # computes the subinstrument longitude
    sublon = asc + math.atan2(cosoi * sinu, cosu)

    # computes the spacecraft to earth fixed coordinates transformation
    # matrix:
    #     (vector in ecef coordinates) = b * (vector in s/c coordinates)
    #b = [[np.nan for i in range(3)], [np.nan for i in range(3)], [np.nan for i in range(3)]]
    b = np.empty((3,3))
    b[:] = np.nan

    b[0][1] = -sinasc * sinoi
    b[1][1] = cosasc * sinoi
    b[2][1] = -cosoi
    b[0][2] = -cosasc * cosu + sinasc * sinu * cosoi
    b[1][2] = -sinasc * cosu - cosasc * sinu * cosoi
    b[2][2] = -slat
    b[0][0] = -cosasc * sinu - sinasc * cosu * cosoi
    b[1][0] = -sinasc * sinu + cosasc * cosu * cosoi
    b[2][0] = cosu * sinoi
    # computes the normalized spacecraft position vector in earth fixed
    # coordinates - xs.
    r = (nomorb + dr) / aec
    xs = [np.nan] * 3  # normalized s/c position in ecef coordinates
    # #bt = [[np.nan for i in range(3)], [np.nan for i in range(3)], [np.nan for i in range(3)]]
    # bt = np.empty((3,3))
    # bt[:] = np.nan
    #print bt
    xs[0] = -b[0][2] * r
    xs[1] = -b[1][2] * r
    xs[2] = -b[2][2] * r

    # precomputes q3 (used in lpoint function (now in navToLatLon() )
    #q3 = xs[0] * xs[0] + xs[1] * xs[1] + aebe2c * xs[2] * xs[2] - 1.0

    # computes the attitudes and misalignments if imc is off
    if imc != 0:
        # computes the solar orbit angle
        wa = b0.daily_solar_rate.python_value() * ts
        # computes the difference between current time, ts, and the  # exponential time, iparms(62). note that both times are since epoch.
        te = ts - b0.exp_start_time_from_epoch.python_value()
        # computes roll + roll misalignment

        roll = roll + compute_attitude_angle(att_angle=b0.roll_angle, solar_orbit_angle=wa, exp_time_delay=te)
        # computes pitch + pitch misalignment

        pitch = pitch + compute_attitude_angle(att_angle=b0.pitch_angle, solar_orbit_angle=wa, exp_time_delay=te)
        # computes yaw

        yaw = yaw + compute_attitude_angle(att_angle=b0.yaw_angle, solar_orbit_angle=wa, exp_time_delay=te)
        # computes roll misalignment

        rma = float(compute_attitude_angle(att_angle=b0.roll_missalignment_angle, solar_orbit_angle=wa, exp_time_delay=te))
        # computes pitch misalignment

        pma = float(compute_attitude_angle(att_angle=b0.pitch_missalignment_angle, solar_orbit_angle=wa, exp_time_delay=te))
        # apply the earth sensor compensation if needed
        roll = roll + b0.compens_roll.python_value()
        pitch = pitch + b0.compens_pitch.python_value()
        yaw = yaw + b0.compens_yaw.python_value()
        # #end if (imc...)


    #inst2e
    rpy = np.empty((3,3), dtype=np.float32)
    rpy[:] = np.nan
    # we compute instrument to body coordinates transformation
    # matrix by using a small angle approximation of trigonometric
    # funktions of the roll, pitch and yaw.

    rollsq = roll*roll
    pitchsq = pitch*pitch
    yawsq = yaw*yaw

    rpy[0][0] = 1. - 0.5 * (pitchsq + yawsq)
    rpy[0][1] = -yaw
    rpy[0][2] = pitch
    rpy[1][0] = yaw + pitch * roll
    rpy[1][1] = 1. - 0.5 * (yawsq + rollsq)
    rpy[1][2] = -roll
    rpy[2][0] = -pitch + roll * yaw
    rpy[2][1] = roll + pitch * yaw
    rpy[2][2] = 1. - 0.5 * (pitchsq + rollsq)

    # # multiplication of matrices b and rpy
    # for i in range(3):
    #     for j in range(3):
    #         bt[i][j] = b[i][0] * rpy[0][j] + b[i][1] * rpy[1][j] + b[i][2] * rpy[2][j]
    #previous code is replaced using an elegant numpy solution
    bt = np.dot(b, rpy)



    #init ends together with inst2e
    lats = np.asarray(lats)
    lons = np.asarray(lons)

    latsm = lats != np.nan
    lonsm = lons != np.nan

    hw = lats[latsm].size
    # f = np.empty((3,h_,w_))
    f = np.empty((3, hw))
    # ft = np.empty((3,h_,w_))
    ft = np.empty((3, hw))
    # u = np.empty((3,h_,w_))
    u = np.empty((3, hw))
    f[:] = np.nan
    ft[:] = np.nan
    u[:] = np.nan
    # lats and lons can contain nans
    ff =  float(iflip)
    doff = scnmax - ewnom
    rlat = lats[latsm] * RAD
    rlon = lons[lonsm] * RAD
    if not isEastPositive:
        rlon = -rlon

    # transform lat/lon to elevation and scan angles
    # (used to be the gpoint routine...)

    # computes sinus of geographic (geodetic) latitude
    sing = np.sin(rlat)
    w1 = aebe4c * sing * sing
    # sinus of the geocentric latitude
    slat = ((0.375 * w1 - 0.5) * w1 + 1.) * sing / aebe2c

    # computes local earth radius at specified point
    w2 = slat * slat
    w1 = aebe3c * w2
    w1 = (0.375 * w1 - 0.5) * w1 + 1.

    # computes cartesian coordinates of the point
    u[2] = slat * w1
    w2 = w1 * np.sqrt(1. - w2)
    u[0] = w2 * np.cos(rlon)
    u[1] = w2 * np.sin(rlon)
    # pointing vector from instrument to the earth point
    f[0] = u[0] - xs[0]
    f[1] = u[1] - xs[1]
    f[2] = u[2] - xs[2]
    w2 = u[0] * f[0] + u[1] * f[1] + u[2] * f[2] * aebe2c
    w2m = ~np.isnan(w2)
    wmask = np.zeros_like(w2).astype(np.bool)
    wmask[w2m] = w2[w2m] <= 0.

    # converts pointing vector to instrument coordinates
    ft[0] = bt[0][0] * f[0] + bt[1][0] * f[1] + bt[2][0] * f[2]
    ft[1] = bt[0][1] * f[0] + bt[1][1] * f[1] + bt[2][1] * f[2]
    ft[2] = bt[0][2] * f[0] + bt[1][2] * f[1] + bt[2][2] * f[2]



    # converts pointing vector to scan and elevation angles and
    # corrects for the roll and pitch misalignments
    gam = np.arctan(ft[0] / np.sqrt(ft[1] * ft[1] + ft[2] * ft[2]))
    alf = -np.arctan(ft[1] / ft[2])
    w1 = np.sin(alf)
    w2 = np.cos(gam)
    alpha1 = alf + rma * (1. - np.cos(alf) / w2) + pma * w1 * (doff / w2 + np.arctan(gam))
    gam = gam - ff * rma * w1
    alf = alpha1 + alpha1 * gam * doff
    gam = gam - 0.5 * alpha1 * alpha1 * doff

    # convert elevation and scan angles to line/pixel coordinates

    # compute fractional line number

    tmplin = (elvmax - alf) / elvln
    tmplin = tmplin + 4.5

    # compute fractional pixel number
    tmpele = (scnmax + gam) / scnpx + 1.



    lin = np.where(wmask, tmplin, 0)
    ele = np.where(wmask, tmpele, 0)
    # print startLine, startElement
    # print magLine, magElement, resLine
    # print startImageLine, startImageElement

    #alin, aele = imageCoordToAreaCoord(linele=[lin, ele])




    alin = .5 + startLine + (magLine * (lin - startImageLine)) / resLine
    if isLineFlipped: # the line flip should be detected at the begining but in the original code i did not find the chunk that calculate the line offset. Probably the line offste is the difference
        alin -= lineOffset
    aele = .5 + startElement + (magElement * (ele - startImageElement)) / resElement

    #print alin, aele


    olin = np.zeros_like(lats)
    oele = np.zeros_like(lons)
    olin[latsm] = alin
    oele[lonsm] = aele
    return olin.squeeze(), oele.squeeze()

def le2ll(lines=None, cols=None, b0=None, lineres=1., colres=1.):
    """
    Converts between relative indices of a GVAR file or so called  file coordinates (similar ot numpy, for 0 to nl/nc-1)
    into latitude longitude given the orbit and attitude info provided in a valid block0 
    By default operates in the coordinates of the VISIBLE channel that is 1 KM
    Additionally a different resolution for each axis can be requested and shoud be honoured
    
    :param lines: numpy array of lines
    :param cols: numpy array of cols
    :param b0: block0 class emulating a block0 C struct from pygvar package 
    :param lineres:  1
    :param colres: 1
    :return: the corresponding lat/lon coordinates
    
    NB , McIDAS uses positive longitude instead of negative in western hemi. Just to confuse people from outside US
    
    
    """

    # init
    # TODO FIGURE OUT the resolution deal inside an area file
    resLine = lineres  # setup by user, default value if 1 (vis chn res is used)
    resElement = colres  # setup by user, default value if 1 (vis chn res is used)
    magLine = 1.
    magElement = 1.
    startLine = 0
    startElement = 0.
    startImageLine = b0.infln  # init from b0
    startImageElement = b0.iwfpx  # init from b0
    isLineFlipped = False
    lineOffset = 0


    #startImageLine, startImageElement = 2580, 10098


    PI = math.pi
    DEG = 180. / PI
    RAD = PI / 180.  # degrees to radians conversion pi/180

    instr = 1  # 1 imager 2 sounder...we do not handle sounder...goodbye
    isEastPositive = True  # this refers to longitude being positive. MCIDAS can mess with logitudes specially when reprojecting
    # earth

    nomorb = 42164.365  # nominal radial distance of instrument (km)
    aec = 6378.137  # earth equatorial radius (km)
    ferc = 1. - (6356.7533 / aec)  # earth flattening coeff = 1-(be/ae)
    aebe2c = 1. / (1 - ferc) ** 2
    #aebe3c = aebe2c - 1.
    #aebe4c = ((1. - ferc) ** 4.) - 1.

    # instrument(IMAGER) specific data these data is only for imager

    incmax = 6136  # number of increments per cycle

    elvmax = 0.220896  # bounds in elevation (radians)
    scnmax = 0.24544  # bounds in scan angle (radians)
    elvinc = 8.e-6  # change in elevation angle per increment (rad)
    scninc = 16.e-6  # change in scan angle per increment (radians)
    elvln = 28.e-6  # elevation angle per detector line (radians)
    scnpx = 16.e-6  # scan angle per pixel (radians)
    nsnom = .220896  # north-south center of instrument (=4.5 x incmax x elvinc)

    ewnom = .24544  # east-west center of instrument (=2.5 x incmax x sncinc)

    nadnsc = b0.ns_cycles
    nadnsi = b0.ns_incr
    nadewc = b0.ew_cycles
    nadewi = b0.ew_incr

    if nadnsc != 0 and nadnsi != 0 and nadewc != 0 and nadewi != 0:
        # only imager is handled , sounder is ignored
        elvmax = (nadnsc * incmax + nadnsi) * elvinc
        # update the scan angle in x (first?))
        scnmax = (nadewc * incmax + nadewi) * scninc

    # print nadnsc, nadnsi,nadewc,nadewi, elvmax, scnmax
    # time computation
    jan_1_1950 = datetime.datetime(1950, 1, 1)
    bcdtime = b0.epochdate.python_value()

    # this is time in milisec from epoch jan 1 1950
    diff = bcdtime - jan_1_1950
    epoch = diff.total_seconds() / 60.
    # currect header block
    # firts convert tched into yyyydoy and hhmmssmmm
    ydoy, hhmmssmmm = datetime2mcidas(b0.tched.python_value())
    # convert an int hhmmssmmm into a float, byte order is hardcoded here for little endian
    hhmmssmmm_f = struct.unpack('{byte_order:s}{size:d}{format_str:s}'.format(byte_order='<', size=1, format_str='f'),struct.pack('{byte_order:s}{size:d}{format_str:s}'.format(byte_order='<', size=1,format_str='i'), hhmmssmmm))[0]
    #compute image tima in seconds from epoch
    year = 1900 + ydoy / 1000
    day = ydoy - ydoy / 1000 * 1000
    hour = int(hhmmssmmm_f / 10000)
    minuts = int(hhmmssmmm_f / 100 - hour * 100)
    secs = hhmmssmmm_f - float(100 * minuts) - float(10000 * hour)

    j = day + 1461 * (year + 4799) / 4 - 3 * ((year + 4899) / 100) / 4 - 2465022

    imgtim = float(j * 1440. + hour * 60. + minuts + (secs / 60.))

    # this perhaps need to be specifically checked for
    imc = 0 if (b0.iscan.imc_active != 0) else  1

    iflip = -1 if b0.iscan.yaw_flip else 1
    # print imc, iflip, imgtim

    # assign reference values to the subsatellite longitude and
    # latitude, the radial distance and the orbit yaw.
    lam = b0.reflo.python_value()

    dr = b0.refra.python_value()

    phi = b0.refla.python_value()

    psi = b0.reforyaw.python_value()
    subpoint = []
    subpoint.append(phi / RAD)
    subpoint.append(lam / RAD)

    # assign reference values to the attitudes and misalignments
    roll = b0.refatroll.python_value()
    pitch = b0.refatpitch.python_value()
    yaw = b0.refatyaw.python_value()

    rma = 0.
    pma = 0.
    if imc != 0:
        # set reference radial distance, latitude and orbit yaw to zero, does it make sens?
        dr = 0.
        phi = 0.
        psi = 0.
        # compute time since epoch (in minutes)
        ts = imgtim - epoch
        # computes orbit angle and the related trigonometric functions.
        # earth rotational rate=.729115e-4 (RAD/s)
        w = 0.729115e-4 * 60.0 * ts
        sw = math.sin(w)
        cw = math.cos(w)
        sw1 = math.sin(0.927 * w)
        cw1 = math.cos(0.927 * w)
        s2w = math.sin(2. * w)
        c2w = math.cos(2. * w)
        sw3 = math.sin(1.9268 * w)
        cw3 = math.cos(1.9268 * w)
        # computes change in the imc longitude from the reference
        change_lon = [e.python_value() for e in b0.change_lon]
        lam = lam + change_lon[0] + (change_lon[1] + change_lon[2] * w) * w \
              + (change_lon[9] * sw1 + change_lon[10] * cw1 + change_lon[3] * sw \
                 + change_lon[4] * cw + change_lon[5] * s2w + change_lon[6] * c2w \
                 + change_lon[7] * sw3 + change_lon[0] * cw3 + w * (change_lon[11] * sw + change_lon[12] * cw)) * 2.

        change_radist = [e.python_value() for e in b0.change_radist]
        # computes change in radial distance from the reference (km)
        dr = dr + change_radist[0] + change_radist[1] * cw + change_radist[2] * sw \
             + change_radist[3] * c2w + change_radist[4] * s2w + change_radist[5] * cw3 \
             + change_radist[6] * sw3 + change_radist[7] * cw1 \
             + change_radist[8] * sw1 + w * (change_radist[9] * cw + change_radist[10] * sw)

        # computes the sine of the change in the geocentric latitude
        change_sin_geoc_lat = [e.python_value() for e in b0.sin_geoc_lat]
        dlat = change_sin_geoc_lat[0] + change_sin_geoc_lat[1] * cw + change_sin_geoc_lat[2] * sw \
               + change_sin_geoc_lat[3] * c2w + change_sin_geoc_lat[4] * s2w + \
               w * (change_sin_geoc_lat[5] * cw + change_sin_geoc_lat[6] * sw) + change_sin_geoc_lat[7] * cw1 + \
               change_sin_geoc_lat[8] * sw1

        # computes geocentric latitude by using an expansion for arcsine
        phi = phi + dlat * (1. + dlat * dlat / 6.)

        change_sin_orbit_yaw = [e.python_value() for e in b0.sin_orbit_yaw]
        # computes sine of the change in the orbit yaw
        dyaw = change_sin_orbit_yaw[0] + change_sin_orbit_yaw[1] * sw + change_sin_orbit_yaw[2] * cw \
               + change_sin_orbit_yaw[3] * s2w + change_sin_orbit_yaw[4] * c2w \
               + w * (change_sin_orbit_yaw[5] * sw + change_sin_orbit_yaw[6] * cw) \
               + change_sin_orbit_yaw[7] * sw1 + change_sin_orbit_yaw[8] * cw1

        # computes the orbit yaw by using an expansion for arcsine.
        psi = psi + dyaw * (1. + dyaw * dyaw / 6.)

        # calculation of changes in the instrument orbit ends here

    # conversion of the imc longitude and orbit yaw to the subinstrument
    # longitude and the orbit inclination (ref: goes-pcc-tm-2473, inputs
    # required for earth location and gridding by sps, june 6, 1988)
    slat = math.sin(phi)
    syaw = math.sin(psi)
    sinoi = slat * slat + syaw * syaw
    cosoi = math.sqrt(1. - sinoi)
    sinoi = math.sqrt(sinoi)

    if slat == 0.0 and syaw == 0.0:
        u = 0.0
    else:
        u = math.atan2(slat, syaw)
    sinu = math.sin(u)
    cosu = math.cos(u)

    # computes longitude of the ascending node
    asc = lam - u
    sinasc = math.sin(asc)
    cosasc = math.cos(asc)

    # computes the subinstrument geographic latitude
    sublat = math.atan(aebe2c * math.tan(phi))

    # computes the subinstrument longitude
    sublon = asc + math.atan2(cosoi * sinu, cosu)

    # computes the spacecraft to earth fixed coordinates transformation
    # matrix:
    #     (vector in ecef coordinates) = b * (vector in s/c coordinates)
    # b = [[np.nan for i in range(3)], [np.nan for i in range(3)], [np.nan for i in range(3)]]
    b = np.empty((3, 3))
    b[:] = np.nan

    b[0][1] = -sinasc * sinoi
    b[1][1] = cosasc * sinoi
    b[2][1] = -cosoi
    b[0][2] = -cosasc * cosu + sinasc * sinu * cosoi
    b[1][2] = -sinasc * cosu - cosasc * sinu * cosoi
    b[2][2] = -slat
    b[0][0] = -cosasc * sinu - sinasc * cosu * cosoi
    b[1][0] = -sinasc * sinu + cosasc * cosu * cosoi
    b[2][0] = cosu * sinoi
    # computes the normalized spacecraft position vector in earth fixed
    # coordinates - xs.
    r = (nomorb + dr) / aec
    xs = [np.nan] * 3  # normalized s/c position in ecef coordinates

    xs[0] = -b[0][2] * r
    xs[1] = -b[1][2] * r
    xs[2] = -b[2][2] * r

    # precomputes q3 (used in lpoint function (now in navToLatLon() )
    q3 = xs[0] * xs[0] + xs[1] * xs[1] + aebe2c * xs[2] * xs[2] - 1.0

    # computes the attitudes and misalignments if imc is off
    if imc != 0:
        # computes the solar orbit angle
        wa = b0.daily_solar_rate.python_value() * ts
        # computes the difference between current time, ts, and the  # exponential time, iparms(62). note that both times are since epoch.
        te = ts - b0.exp_start_time_from_epoch.python_value()
        # computes roll + roll misalignment

        roll = roll + compute_attitude_angle(att_angle=b0.roll_angle, solar_orbit_angle=wa, exp_time_delay=te)
        # computes pitch + pitch misalignment

        pitch = pitch + compute_attitude_angle(att_angle=b0.pitch_angle, solar_orbit_angle=wa, exp_time_delay=te)
        # computes yaw

        yaw = yaw + compute_attitude_angle(att_angle=b0.yaw_angle, solar_orbit_angle=wa, exp_time_delay=te)
        # computes roll misalignment

        rma = float(
            compute_attitude_angle(att_angle=b0.roll_missalignment_angle, solar_orbit_angle=wa, exp_time_delay=te))
        # computes pitch misalignment

        pma = float(
            compute_attitude_angle(att_angle=b0.pitch_missalignment_angle, solar_orbit_angle=wa, exp_time_delay=te))
        # apply the earth sensor compensation if needed
        roll = roll + b0.compens_roll.python_value()
        pitch = pitch + b0.compens_pitch.python_value()
        yaw = yaw + b0.compens_yaw.python_value()
        # #end if (imc...)

    # inst2e


    rpy = np.empty((3, 3), dtype=np.float32)
    rpy[:] = np.nan
    # we compute instrument to body coordinates transformation
    # matrix by using a small angle approximation of trigonometric
    # funktions of the roll, pitch and yaw.

    rollsq = roll * roll
    pitchsq = pitch * pitch
    yawsq = yaw * yaw

    rpy[0][0] = 1. - 0.5 * (pitchsq + yawsq)
    rpy[0][1] = -yaw
    rpy[0][2] = pitch
    rpy[1][0] = yaw + pitch * roll
    rpy[1][1] = 1. - 0.5 * (yawsq + rollsq)
    rpy[1][2] = -roll
    rpy[2][0] = -pitch + roll * yaw
    rpy[2][1] = roll + pitch * yaw
    rpy[2][2] = 1. - 0.5 * (pitchsq + rollsq)

    # # multiplication of matrices b and rpy
    # for i in range(3):
    #     for j in range(3):
    #         bt[i][j] = b[i][0] * rpy[0][j] + b[i][1] * rpy[1][j] + b[i][2] * rpy[2][j]
    # previous code is replaced
    bt = np.dot(b, rpy)



    # init ends together with inst2e

    #rl, rp = self.areaCoordtoImageCoord_np(linele=[cols, lines])
    # rl, rp = lines, cols

    rl = lineOffset - lines if isLineFlipped else lines
    rl = startImageLine + (resLine * (rl - startLine)) / magLine
    rp = startImageElement + (resElement * (cols - startElement)) / magElement

    rl = np.asarray(rl, dtype=np.float32)
    rp = np.asarray(rp, dtype=np.float32)

    if rl.ndim == 2 and rp.ndim == 2:

        cshape = 3, rl.shape[0], rp.shape[1]
    elif rl.ndim == 1 and rp.ndim == 1:
        cshape = 3,rl.size,rp.size
    elif rl.ndim == 0 and rp.ndim == 0:
        cshape = 3, rl.size, rp.size
    else:
        raise ValueError('input lines and columns can havve maximum 2 dims. lines have %s and cols have %s' % (rl.ndim, rp.ndim))
    g1 = np.empty(cshape)
    g1[:] = np.nan
    g = np.empty(cshape)
    g[:] = np.nan

    u = np.empty(cshape)
    u[:] = np.nan
    #  compute elevation and scan angles (e,s) related to input
    #  line and pixel numbers
    alpha0 = elvmax - (rl - 4.5) * elvln
    zeta0 = (rp - 1.0) * scnpx - scnmax

    # compute sign of misalignment corrections and origin offset
    ff =  float(iflip)
    doff = scnmax - ewnom

    # add new second order origin offset correction

    alpha = alpha0 - alpha0 * zeta0 * doff
    zeta = zeta0 + 0.5 * alpha0 * alpha0 * doff
    #  transform elevation and scan angles to geographic coordinates
    #  (this is the old 'lpoint' routine...

    # computes trigonometric funktions of the scan and elevation
    # angles corrected for the roll and pitch misalignments
    ca = np.cos(alpha)
    sa = np.sin(alpha)
    cz = np.cos(zeta)
    da = alpha - pma * sa * (ff / cz + np.tan(zeta)) - rma * (1.0 - ca / cz)
    dz = zeta + ff * rma * sa

    # corrected scan angle
    cz = np.cos(dz)
    # computes pointing vector in instrument coordinates
    g[0] = np.sin(dz)
    g[1] = -cz * np.sin(da)
    g[2] = cz * np.cos(da)

    # transforms the pointing vector to earth fixed coordinates
    g1[0] = bt[0][0] * g[0] + bt[0][1] * g[1] + bt[0][2] * g[2]
    g1[1] = bt[1][0] * g[0] + bt[1][1] * g[1] + bt[1][2] * g[2]
    g1[2] = bt[2][0] * g[0] + bt[2][1] * g[1] + bt[2][2] * g[2]

    # computes coefficients and solves a quadratic equation to
    # find the intersect of the pointing vector with the earth
    # surface
    q1 = g1[0] * g1[0] + g1[1] * g1[1] + aebe2c * g1[2] * g1[2]
    q2 = xs[0] * g1[0] + xs[1] * g1[1] + aebe2c * xs[2] * g1[2]
    d = q2 * q2 - q1 * q3

    # if the discriminant of the equation, d, is negative, the
    # instrument points off the earth
    dm = np.abs(d) < 1. - 9
    d = np.where(dm, 0, d)
    dm1 = d >= 0

    d = np.sqrt(d)

    # slant distance from the instrument to the earth point
    h = -(q2 + d) / q1

    # cartesian coordinates of the earth point
    u[0] = xs[0] + h * g1[0]
    u[1] = xs[1] + h * g1[1]
    u[2] = xs[2] + h * g1[2]

    # sinus of geocentric latitude
    d1 = u[2] / np.sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2])

    rlat = np.where(dm1, np.arctan(aebe2c * d1 / np.sqrt(1. - d1 * d1)), np.nan)
    rlon = np.where(dm1, np.arctan2(u[1], u[0]), np.nan)




    rlat = rlat * DEG
    rlon = rlon * DEG

    #  put longitude into mcidas form
    if not isEastPositive:
        rlon = -rlon
    return rlat.squeeze(), rlon.squeeze()



