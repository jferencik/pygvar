import ctypes as C
import datetime
from pprint import pformat
from pygvar.utils import f2int
import math
PI = math.pi
DEG = 180. / PI
RAD = PI / 180.  # degrees to radians conversion pi/180


class GvarStruct(C.BigEndianStructure):

    _pack_ = 1
    _fields_ = ()

    def python_value(self):
        pass

    @property
    def buffer(self):
        return C.string_at(C.byref(self), C.sizeof(self))

    @property
    def fields(self):
        d = {}

        for fdef in self._fields_:
            v = getattr(self, fdef[0])

            try:
                d[fdef[0]] = v.fields
            except AttributeError as ae:
                # print fdef, ae
                cls_name = fdef[1].__class__.__name__
                if 'Array' in cls_name:
                   #d[fdef[0]] = dict([('%s[%d]' % (fdef[0],i), repr(e)) for i, e in enumerate(v)])
                   d[fdef[0]] = [e for e in v]
                else:
                    d[fdef[0]] = v




        return d

    @classmethod
    def size(cls):
        return C.sizeof(cls)

    def __repr__(self):

        return pformat(dict( (e[0], getattr(self, e[0])) for e in self._fields_ ))
        # d ={}
        # for e in self._fields_:
        #     v = getattr(self, e[0])
        #     if hasattr(v, 'python_value',):
        #         d[e[0]] = v.python_value()
        #     else:
        #         d[e[0]] = v
        # return pformat(d)


    def __eq__(self, other):
        rv = True
        for fn, fv in  self.fields.items():
            lv = fv == other.fields[fn]
            rv = rv&lv
        return rv

class BCDDateTime(GvarStruct):
    _field_names_ = names = 'year_1000', 'year_100', 'year_10', 'year_1', 'day_100', 'day_10', 'day_1', 'hour_10', 'hour_1', 'min_10', 'min_1', 'sec_10', 'sec_1', 'msec_100', 'msec_10', 'msec_1'
    _field_types_ = [C.c_uint8]* len(_field_names_)
    _field_l_ = [4]* len(_field_names_)

    _fields_ = tuple(zip(_field_names_, _field_types_, _field_l_))

    # this does not work in python3
    #_fields_ = tuple([(_field_names_[i], _field_types_[i], 4) for i in range(len(_field_names_))])


    @staticmethod
    def get_month_day(year, day, one_based=False):
        if one_based:  # if Jan 1st is 1 instead of 0
                day -= 1
        dt = datetime.datetime(year, 1, 1) + datetime.timedelta(days=day)
        return dt.month, dt.day

    def python_value(self):
        # compose year
        d = self.fields
        year = d['year_1000'] * 1000 + d['year_100'] * 100 + d['year_10'] * 10 + d['year_1']
        yday = d['day_100'] * 100 + d['day_10'] * 10 + d['day_1']
        month, day = self.get_month_day(year, yday, True)
        hour = d['hour_10'] * 10 + d['hour_1']
        minute = d['min_10'] * 10 + d['min_1']
        second = d['sec_10'] * 10 + d['sec_1']
        milisecond = d['msec_100'] * 100 + d['msec_10'] * 10 + d['msec_1']
        return datetime.datetime(year, month, day, hour, minute, second, milisecond*1000)

class ISCAN(GvarStruct):
    """
        ' 0 Frame Start',
        ' 1 Frame End',
        ' 2 Frame break ( lines lost )',
        ' 3 Line break ( Pixels lost )',
        ' 4 Priority 1 frame data ',
        ' 5 Priority 2 frame data ',
        ' 6 East-to-west scan',
        ' 7 South-to-north scan',
        ' 8 IMC active',
        ' 9 Lost header block',
        '10 Lost trailer block',
        '11 Lost telemetry data',
        '12 (Star sense) time break',
        '13 Side 2 (secondary) active',
        '14 Visible normalization active',
        '15 IR calibration active',
        '16 yaw-flipped mode (GOES-10)',
        '17 IR detector 1 data not valid',
        '18 IR detector 2 data not valid',
        '19 IR detector 3 data not valid',
        '20 IR detector 4 data not valid',
        '21 IR detector 5 data not valid',
        '22 IR detector 6 data not valid',
        '23 IR detector 7 data not valid',
        '24 Visible detector 1 data not valid',
        '25 Visible detector 2 data not valid',
        '26 Visible detector 3 data not valid',
        '27 Visible detector 4 data not valid',
        '28 Visible detector 5 data not valid',
        '29 Visible detector 6 data not valid',
        '30 Visible detector 7 data not valid',
        '31 Visible detector 8 data not valid'
    """
    _fields_ = (
        ('frame_start', C.c_uint8, 1), #0
        ('frame_end', C.c_uint8, 1), #1
        ('frame_break', C.c_uint8, 1), #2
        ('line_break', C.c_uint8, 1), #3
        ('priority1', C.c_uint8, 1), #4
        ('priority2', C.c_uint8, 1), #5
        ('west_to_east', C.c_uint8, 1), #6
        ('north_to_south', C.c_uint8, 1), #7
        ('imc_active', C.c_uint8, 1), #8
        ('header_lost', C.c_uint8, 1), #9
        ('trailer_lost', C.c_uint8, 1), #10
        ('telemetry_lost', C.c_uint8, 1), #11
        ('time_break', C.c_uint8, 1), #12
        ('side2', C.c_uint8, 1), #13
        ('vis_normalization', C.c_uint8, 1), #14
        ('ir_calibration', C.c_uint8, 1), #15
        ('yaw_flip', C.c_uint8, 1), #16
        ('irdet1_invalid', C.c_uint8, 1),
        ('irdet2_invalid', C.c_uint8, 1),
        ('irdet3_invalid', C.c_uint8, 1),
        ('irdet4_invalid', C.c_uint8, 1),
        ('irdet5_invalid', C.c_uint8, 1),
        ('irdet6_invalid', C.c_uint8, 1),
        ('irdet7_invalid', C.c_uint8, 1),
        ('visdet1_invalid', C.c_uint8, 1),
        ('visdet2_invalid', C.c_uint8, 1),
        ('visdet3_invalid', C.c_uint8, 1),
        ('visdet4_invalid', C.c_uint8, 1),
        ('visdet5_invalid', C.c_uint8, 1),
        ('visdet6_invalid', C.c_uint8, 1),
        ('visdet7_invalid', C.c_uint8, 1),
        ('visdet8_invalid', C.c_uint8, 1),
    )

    def bstring(self):
        return ''.join([str(getattr(self, e[0])) for e in self._fields_])

class GouldFloat(GvarStruct):

    _fields_ = (
        ('int_value', C.c_uint32),
    )


    def python_value(self):
           #the input is an 4 byte int of a hexadecimal float in Gould format
        #we need to prepare byte masks in order to extract the most significant byte (byte 4)
        # as well the first three bytes. I am not sure but my feeling is that OGE doc implies Gould is
        #big endian byte order and thus the section 3.5.4 labels bytes in an opposite way (byte1 is the most significant byte)
        #however the nav header is in Little endian and I use its labeling for bytes
        #prepare masks
        byte012_mask = int('0xffffff',16) # 0000 0000 1111 1111 1111 1111 1111 1111
        byte3_mask = int('0xFFFFFFFF', 16) # 1111 1111 0000 0000 0000 0000 0000 0000
        sign_bit_mask = int('0x80000000',16) # 10000000000000000000000000000000
        exponent_mask = int('0x7F000000', 16) # '01111111000000000000000000000000
        #compute sign
        sign = -1 if self.int_value & sign_bit_mask else 1
        if sign < 0:
            #two's complement for negative numbers
            number = (self.int_value ^ byte3_mask) + 1
        else:
            number = self.int_value
        #compute exponent
        exponent = (number & exponent_mask) >> 24
        if exponent == 0: #account for bias
            exponent = 64
        #compute matissa
        mantissa = number & byte012_mask
        tempVal = 16 ** float(70 - exponent)
        nativeVal = mantissa / tempVal
        nativeVal = nativeVal * sign

        return nativeVal


class Polar(GvarStruct):
    _fields_ = (
        ('mg', GouldFloat),
        ('ph', GouldFloat),
    )



class MonomialSinusoid(GvarStruct):
    _fields_ = (
        ('order_of_applicable', C.c_uint32),
        ('order_of', C.c_uint32),
        ('sinusoid', Polar),
        ('angle_of_epoch_where_monomial_is_zero', GouldFloat),
    )

class AttitudeAngle(GvarStruct):
    """
    SelFloat         Exponential_magnitude; /* 523  526*/
  SelFloat         Exponential_time_constant; /* 527  530*/
  SelFloat         Constant_mean_attitude_angle; /* 531  534*/
  uint32           Number_of_sinusoidals_per_angle; /* 535  538*/
  Polar            Sinusoid[15];		/* 539  658*/
  uint32           Number_of_monomial_sinusoids; /* 659  662*/
  MonomialSinusoid Monomial[4]; /* 663  742*/
    """
    _fields_ = (
        ('exp_mag', GouldFloat),
        ('exp_time_const', GouldFloat),
        ('const_mean_attit_angle', GouldFloat),
        ('nsinusoidals', C.c_uint32),
        ('sinusoids', Polar*15),
        ('n_monomial_sinusoids', C.c_uint32),
        ('monomials', MonomialSinusoid*4),
    )
    def to_mcidas(self):
        out_angle = list()
        out_angle.append(f2int(self.exp_mag.python_value() * 1e7))
        out_angle.append(f2int(self.exp_time_const.python_value() * 1e2))
        out_angle.append(f2int(self.const_mean_attit_angle.python_value()* 1e7))
        out_angle.append(f2int(self.nsinusoidals))

        # sinusoidals
        for _s in self.sinusoids:

            out_angle.append(f2int(_s.mg.python_value() * 1e7))
            out_angle.append(f2int(_s.ph.python_value() * 1e7))

        out_angle.append(f2int(self.n_monomial_sinusoids))

        # Monomials
        for m in self.monomials:
            out_angle.append(f2int(m.order_of_applicable))
            out_angle.append(f2int(m.order_of))
            out_angle.append(f2int(m.sinusoid.mg.python_value() * 1e7))
            out_angle.append(f2int(m.sinusoid.ph.python_value() * 1e7))
            out_angle.append(f2int(m.angle_of_epoch_where_monomial_is_zero.python_value() * 1e7))

        return out_angle




class BiasStatistics(GvarStruct):


    _fields_ = (
        ('tod_of_bias_stats', BCDDateTime),
        ('total_sample_size', C.c_uint16*7),
        ('filtered_sample_size', C.c_uint16*7),
        ('unfiltered_min_val', C.c_uint16*7),
        ('filtered_min_val', C.c_uint16*7),
        ('unfiltered_max_val', C.c_uint16*7),
        ('filtered_max_val', C.c_uint16*7),
        ('unfiltered_mean_val', GouldFloat*7),
        ('filtered_mean_val', GouldFloat*7),
        ('unfiltered_stddev', GouldFloat*7),
        ('filtered_sigma_counts', GouldFloat*7),
        ('filtered_sigma_radiance', GouldFloat*7),
        ('filtered_sigma_temp', GouldFloat*7),
        ('clamp_mode', C.c_int32),

    )

class Calibration(GvarStruct):

    _fields_ = (
        ('L1', GouldFloat*6),
        ('L2', GouldFloat*6),
        ('BP1', GouldFloat*6*6),
        ('BB', GouldFloat*8*6),
        ('OP', GouldFloat*9*6),
        ('PA', GouldFloat*6),
        ('PB', GouldFloat*6),
    )

class block_header(GvarStruct):

    _fields_ =(

        ('block_id', C.c_uint8),
        ('wsize', C.c_uint8),
        ('wcount', C.c_int16),
        ('product_id', C.c_uint16),
        ('repeat_flag', C.c_uint8),
        ('version', C.c_uint8),
        ('data_valid', C.c_uint8),
        ('ascii_binary', C.c_uint8),
        ('spss_id', C.c_int8),
        ('range', C.c_uint8),
        ('block_count', C.c_uint16),
        ('spares', C.c_uint16),
        ('bcd_bytes', BCDDateTime),
        ('spares1', C.c_uint32),
        ('crc', C.c_uint16)



    )

    MAGIC_CRC = 0x1d0f
    def nbytes(self):

        return int(self.wsize * float(self.wcount-2))/8 + 2


class block0(GvarStruct):
    nbytes = 8040

    _fields_ = (
        #('headers', block_header*3),
        ('spcid', C.c_uint8),
        ('spsid', C.c_uint8),
        ('iscan', ISCAN),
        ('idsub', C.c_uint8*16),
        ('tcurr', BCDDateTime),
        ('tched', BCDDateTime),
        ('tctrl', BCDDateTime),
        ('tlhed', BCDDateTime),
        ('tltrl', BCDDateTime),
        ('tipfs', BCDDateTime),
        ('tinfs', BCDDateTime),
        ('tispc', BCDDateTime),
        ('tiecl', BCDDateTime),
        ('tibbc', BCDDateTime),
        ('tistr', BCDDateTime),
        ('tiran', BCDDateTime),
        ('tiirt', BCDDateTime),
        ('tivit', BCDDateTime),
        ('tclmt', BCDDateTime),
        ('tiona', BCDDateTime),
        ('risct', C.c_uint16),
        ('aisct', C.c_uint16),
        ('insln', C.c_uint16), # northernmost visible detector scan line in the current scan.

        ('iwfpx', C.c_uint16),#start_pixel
        ('iefpx', C.c_uint16),#end_pixel
        ('infln', C.c_uint16),#start line
        ('isfln', C.c_uint16),#end_line
        ('imdpx', C.c_uint16),
        ('imdln', C.c_uint16),
        ('imdct', C.c_uint16),

        #The following four terms (words 171-182) are computed using the current O&A set. If IMC is
        #active, the terms reflect the reference subsatellite point position. If IMC is off, the terms reflect
        #the actual subsatellite point.

        ('igvln', C.c_uint16),
        ('igvpx', C.c_uint16),
        ('subla', GouldFloat),
        ('sublo', GouldFloat),
        ('czone', C.c_uint8),
        ('v1phy', C.c_uint8),#The physical detector number (1-8) assigned to GVAR Block 3, could be important to establish the detectors
        ('g1cnt', C.c_uint16),
        ('g2cnt', C.c_uint16),
        ('pbias', C.c_int16),
        ('lbias', C.c_int16),
        ('iscp1', C.c_uint8),
        ('spare_194', C.c_uint8),
        ('idber', GouldFloat),
        ('range', GouldFloat),
        ('gpath', GouldFloat),
        ('xmsne', GouldFloat),
        ('tgpat', BCDDateTime),
        ('txmsn', BCDDateTime),
        ('istim', C.c_int16),
        ('ifram', C.c_uint8),
        ('imode', C.c_uint8),
        #The following four floating point values are in units of degrees. Off-earth coordinates have a value of 999999.
        ('ifnw1', GouldFloat),
        ('ifnw2', GouldFloat),
        ('ifse1', GouldFloat),
        ('ifse2', GouldFloat),
        ('ig2tn', C.c_uint8),
        ('iscp2', C.c_uint8),
        ('isca2', C.c_uint16), # only first 2 words of iscan
        ('t_frame_start', BCDDateTime),
        ('spares_259_277', C.c_uint8*(277-258)),
        ('lon_parity', C.c_uint8),
        #Q&A
        ('imcid', C.c_char*4),
        ('spare_283_294', C.c_uint8*12),
        ('reflo', GouldFloat),
        ('refra', GouldFloat),
        ('refla', GouldFloat),
        ('reforyaw', GouldFloat),
        ('refatroll', GouldFloat),
        ('refatpitch', GouldFloat),
        ('refatyaw', GouldFloat),
        ('epochdate', BCDDateTime),
        ('imc_enable_from_epoch', GouldFloat),
        ('compens_roll', GouldFloat),
        ('compens_pitch', GouldFloat),
        ('compens_yaw', GouldFloat),
        ('change_lon', GouldFloat*13),
        ('change_radist', GouldFloat*11),
        ('sin_geoc_lat', GouldFloat*9),
        ('sin_orbit_yaw', GouldFloat*9),
        ('daily_solar_rate', GouldFloat),
        ('exp_start_time_from_epoch', GouldFloat),
        ('roll_angle', AttitudeAngle),
        ('pitch_angle', AttitudeAngle),
        ('yaw_angle', AttitudeAngle),
        ('roll_missalignment_angle', AttitudeAngle),
        ('pitch_missalignment_angle', AttitudeAngle),
        #OGE is correct(isca3), cgvar says/deoes nothing
        ('isca3', C.c_uint16),
        ('iscp3', C.c_uint8),
        ('parity_279_1265', C.c_uint8),

        ('coregistration_tbl_id', C.c_char*4),
        ('east_west_hourly_corr_terms', C.c_ubyte*48),
        ('index_of_active_corr_terms', C.c_ubyte),
        ('spares_1680_1690', C.c_ubyte*11),

        ('w1691_782', C.c_uint16 *46), #Current scan raw header data block
        ('w1783_874', C.c_uint16 *46),# Current scan raw trailer data block
        ('w1875_966', C.c_uint16 *46), #Latest lagged scan raw header data block
        ('w1967_2058', C.c_uint16 *46), #Latest lagged scan raw trailer data block

        ('w2059_80', C.c_ubyte* 22), #Block 1 command register 1 report, B(20-129)
        ('w2081_102', C.c_ubyte* 22), #Block 2 command register 1 report, B(20-129)
        ('w2103_24', C.c_ubyte* 22), #Block 3 command register 1 report, B(20-129)
        ('w2125_46', C.c_ubyte* 22), #Block 4 command register 1 report, B(20-129)
        ('block_telemetry', C.c_uint16*2*39),
        ('spares_2303_5', C.c_ubyte*3),
        ('parity_1691_2305', C.c_ubyte),

        ('grid1_detector', C.c_ubyte*512),
        ('grid2_detector', C.c_ubyte*512),
        ('grid1_pixel', C.c_uint16*512),
        ('grid2_pixel', C.c_uint16*512),
        ('gridset1_rev_revel', C.c_uint16),
        ('gridset2_rev_revel', C.c_uint16),
        ('w5383_85', C.c_ubyte*3),
        ('parity_2307_5385', C.c_ubyte),

        ('w5387_478', C.c_uint16 * 46),# Oldest lagged scan raw header data block
        ('w5479_570', C.c_uint16 * 46),# Oldest lagged scan raw trailer data block
        ('tophed', BCDDateTime), #time of Oldest lagged scan raw header data block
        ('totrl', BCDDateTime), #time of Oldest lagged scan raw trailer data block

        ('iwbias', GouldFloat*7),# IR calibration bias term
        ('igain1', GouldFloat*7), #IR calibration 1st order gain
        ('igain2', GouldFloat*7), # IR calibration 2nd order gain
        ('ibrate', GouldFloat*7), #IR calibration bias rate
        # #BIASn = IWBIAS + (N - 1) IBRATE
        ('cda_tod_west_pixel', BCDDateTime),
        ('imager_ir_clamped', BiasStatistics),
        ('imager_ir_drift', BiasStatistics),
        ('ns_halfh_correction_terms', C.c_byte*48),
        ('scan_clamp_e_w_clipping_edge_limb_offset', C.c_byte*4),
        ('w6287_303', C.c_byte*17),
        ('parity_5387_6303', C.c_byte),

            #6305-8040 Factory Parameters/8th IR detector Calibration Data
        ('ns_cycles',  C.c_byte),
        ('ew_cycles',  C.c_byte),
        ('ns_incr',  C.c_uint16),
        ('ew_incr',  C.c_uint16),

        ('vis_detector_x_offset',  C.c_uint16*8),
        ('ir_detector_x_offset',  C.c_uint16*14),
        ('vis_detector_y_offset',  C.c_uint16*8),
        ('ir_detector_y_offset',  C.c_uint16*14),

        ('ivcrb',  GouldFloat*8),
        ('ivcr1',  GouldFloat*8),
        ('ivcr2',  GouldFloat*8),
        ('ivral',  GouldFloat),
        ('iicrb',  GouldFloat *2*7),
        ('iicr1',  GouldFloat *2*7),
        ('iicr2',  GouldFloat *2*7),
        ('iisfb',  GouldFloat *2*7),
        ('iisf1',  GouldFloat *2*7),
        ('ig2it',  GouldFloat *2*7*4),
        ('ig2bp',  GouldFloat *4),
        ('ibbtr',  GouldFloat *2*7*4),

        # ('w7234_366',  C.c_ubyte*124),
        # ('imager_cal',  Calibration),
        # ('patch_control_voltage_gain',  GouldFloat),
        # ('patch_control_voltage_bias',  GouldFloat),
        # ('instr_curr_gain',  GouldFloat),
        # ('instr_curr_bias',  GouldFloat),
        # ('w8031_39', C.c_byte*8),
        # ('parity_6305_8039', C.c_byte),
        # ('crc', C.c_uint16)

        #FIX for GOES14 onward as a result of a bug
        ('w7234_366', C.c_ubyte * 124),
        ('imager_cal', Calibration),
        ('patch_control_voltage_gain', GouldFloat),
        ('patch_control_voltage_bias', GouldFloat),
        ('instr_curr_gain', GouldFloat),
        ('instr_curr_bias', GouldFloat),
        ('w8031_32', C.c_byte * 2),
        ('iofnc', C.c_byte ),
        ('iofec', C.c_byte ),
        ('iofni', C.c_short ),
        ('iofei', C.c_short ),
        ('parity_6305_8039', C.c_byte),
        ('crc', C.c_uint16)

        #finally



    )

    @property
    def subpoint(self):
        lon = self.reflo.python_value() / RAD
        lat = self.refla.python_value() / RAD
        return lat, lon
    @property
    def metadata(self):

        s = 'vis bbox: %s \n' % str(self.vis_bbox)
        s += 'ssp: %s \n' % (str(self.subpoint))
        s += 'channels \n'
        for chn_no, sr in self.channels_shape.items():
            s+= '\t%d : %d lines %d columns\n' % (chn_no, sr[0], sr[1])
        return s
    @property
    def vis_bbox(self):
        return self.infln, self.isfln, self.iwfpx, self.iefpx


    @property
    def vis_lines(self):
        return self.isfln - self.infln

    @property
    def vis_cols(self):
        return self.iefpx - self.iwfpx



class datablock(GvarStruct):
    _fields_ = (


    )

class cds_time(C.BigEndianStructure):
    _fields_ = (
        ('days', C.c_uint16),
        ('msecs', C.c_uint32),
    )

    def python_value(self):
        return  datetime.datetime(1958, 1, 1) + datetime.timedelta(days=self.days, milliseconds=self.msecs)







