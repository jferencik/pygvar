'''
Constants for GOES from

 https://www.ospo.noaa.gov/Operations/GOES/calibration/gvar-conversion.html
 https://github.com/zflamig/gvartool/tree/master/src/include



VIS Calibration
_____________________________________________________________________________________________________
  There are three steps to convert a
   10-bit GVAR count value (0-1023) to albedo

Step 1: Convert the GVAR count value to a radiance.
Step 2: Convert radiance to albedo


 Step 1
Although the visible detector-channels are not calibrated in orbit, calibration coefficients measured by
ITT before launch[1] are transmitted to users in the GVAR data stream.
A factor that converts radiance to reflectance factor, or effective albedo, is included as well.
Since detector responsivities can vary unpredictably, the pre-launch calibration may not be valid after launch.
The calibration equation is either Equation (1)

		R = m X + b			(1)

or Equation (2)

 		R = m (X -Xsp),			(2)

where X is the instrument output in counts, the subscript sp refers to the view of space, and m and b are the calibration
coefficients measured before launch[1].
For each visible-channel detector, the radiance R is the average of the spectral (monochromatic) radiance over the
spectral response function for that detector, i.e.,

  R = (integral of R(lambda)PHI(lambda) )/(integral of PHI(lambda) )    	(3)

where lambda is wavelength in �m, PHI the spectral response function, and R(lambda ) the spectral radiance of the target.
The units of R are W/(m2-sr-�m).

The value of b in Equation (1) depends on the electronic zero level.
This level varies with a standard deviation of approximately one count for the imagers and tens of counts for the sounders.
Therefore, when the satellite is in orbit, the value of b determined in the laboratory (or at any other earlier time)
may not be valid. Equation (1) is preferred[2].
Visible-channel data from the GOES-8 and -9 imagers (but not GOES-10 through GOES-13, nor any sounder) are being normalized.
Since normalization makes the responses from all eight imager detectors the same as that of the reference detector,
users of the pre-launch calibration of the GOES-8 and -9 imagers should apply the calibration coefficients for the
reference detector (identified in Table 1) to the data from all detectors.
With relativization enabled, an instrument's output at the view of space is modified slightly, and the value of b in the
calibration equation (1) needs to be modified accordingly. For this reason also, the best approach is to use Equation (2).
If Equation (1) must be used, then the value of b should be determined from the equation

  	         b = - m Xo,			(4)

in which, for the imagers, m is the slope for the reference detector.
For the sounders, it is the slope for an individual detector.
Values of b determined in this way, as well as the values of m and X0, appear in Tables 1 - 6.

Step 2

The reflectance factor (or effective albedo) is obtained from the radiance by

		A = k R,				(5)  where

		k = p / H,			        (6)

and where H is the solar spectral irradiance H(lamda) averaged over the spectral
response function of the visible detector, i.e.,

 h = (integral of H(lambda)PHI(lambda) )/(integral of PHI(lambda) )      (7)

Values of H were computed by ITT[1] from tables of solar irradiance vs wavelength
provided by Rossow et al.[3], whose values are based on measurements by
Neckel and Labs[4].
The values of A lie between 0 and 1. The value of 1 corresponds to the
radiance of a perfectly reflecting diffuse surface illuminated at
normal incidence when the sun is at its annual-average distance from the Earth.
Values of k appear in Tables 1 - 6.


IR Calibration
_____________________________________________________________________________________________________

  There are three steps to convert a
   10-bit GVAR count value (0-1023) to temperature.

Step 1: Convert the GVAR count value to a radiance.
Step 2: Convert radiance to effective temperature using the inverse of the Planck function.
Step 3: Convert effective temperature Teff to actual temperature T (K)


 1.Conversion of Imager GVAR Count to Scene Radiance

A 10-bit GVAR count value (0-1023) can be converted
to a scene radiance according to the following equation:

R = (X - b)/m,
(1)

where R is radiance (mW/[m2-sr-cm-1]) and
      X is the GVAR count value.

The coefficients m and b are the scaling slope and intercept,
respectively.  The values of m and b are listed in Table 1.
 They depend on the channel selected, but for a given channel
they are constant for all time and are the same for all
satellites of the series.


Step 2: Convert radiance to effective temperature using the inverse
of the Planck function as follows:

                (c2 * n )
Teff =   _____________________________
            ln [1 + (c1 * n 3) / R]
where
        c1 = 1.191066 x 10-5 [mW/(m2-sr-cm-4)]
        c2 = 1.438833 (K/cm-1)


where Teff is effective temperature (K),
      ln stands for natural logarithm, and
      R is radiance.

The coefficients n, c1, and c2 are the central wavenumber of the channel
and the two radiation constants, respectively.

The constants c1 and c2 are invariant, but n depends on the
spectral characteristics of the channel and will vary from
instrument to instrument.

Step 3: Convert effective temperature Teff to
actual temperature T (K) using the following equation:

T = a + b * Teff                     (3)
where a and b are two conversion coefficients.

Note in the conversions that:

The values of n (cm-1) in step 2 and constants a and b
in step 3 depend on channel and instrument. Their
values are listed below in Tables 2-1 through 2-5.

The term side 1 or side 2 in the table headings indicates
the operation of one of the two redundant sets of detectors
and electronics on each imager. The coefficients n, a, and b
depend on the choice of side.

The GOES-8, -9, -11, -12 and -13 imagers have always operated on side 1.
The GOES-10 imager is operating on side 2.
The GOES-O imager is expect to operate on side 1 after launch on April 28, 2009.

'''


import numpy as np

X0 = 29
CHANNELS = range(1,7)
VIS_CALIB_COEFFS = {
                'GOES08': {
                    'm':0.5501873,
                    'x0':X0,
                    'k': 1.92979e-3,
                },
                'GOES09': {
                    'm':0.5492361,
                    'x0':X0,
                    'k': 1.94180e-3,
                },
                'GOES10': {
                    'm': (0.5605602, 0.5563529, 0.5566574, 0.5582154, 0.5583361, 0.5571736, 0.5563135, 0.5613536),
                    'x0': X0,
                    'k': 1.98808e-3,
                },
                'GOES11': {
                    'm': (0.5605602, 0.5563529, 0.5566574, 0.5582154, 0.5583361, 0.5571736, 0.5563135, 0.5613536),
                    'x0': X0,
                    'k': 2.01524e-3,
                },
                'GOES12': {
                    'm': (0.5771030, 0.5761764, 0.5775825, 0.5790699, 0.5787051, 0.5755969, 0.5753973, 0.5752099 ),
                    'x0': X0,
                    'k': 1.97658e-3,
                },
                'GOES13': {
                    'm': (0.6120196, 0.6118504, 0.6096360, 0.6087055, 0.6132860, 0.6118208, 0.6122307, 0.6066968),
                    'x0': X0,
                    'k': 1.89544e-3,
                },
                'GOES14': {
                    'm': ( 0.5874693, 0.5865367, 0.5862807, 0.5864086, 0.5857146, 0.5852004, 0.5860814, 0.5841697),
                    'x0': X0,
                    'k': 1.88772e-3,
                },
                'GOES15': {
                    'm': ( 0.5851966, 0.5879772, 0.5856793, 0.5854250, 0.5866992, 0.5836241, 0.5846555, 0.5843753),
                    'x0': X0,
                    'k': 1.88852e-3,
                },

}



#step 1
IR_SLOPEOFFS_COEFFS = {

    2: {'m': 227.3889, 'b': 68.2167},
    3: {'m':  38.8383, 'b': 29.1287},
    4: {'m':  5.2285, 'b': 15.6854 },
    5: {'m':  5.0273, 'b': 15.3332 },
    6: {'m': 5.5297, 'b': 16.5892 }

}
#step 2
c1 = 1.191066e-5  #[mW/(m2-sr-cm-4)]
c2 = 1.438833 # (K/cm-1)

#step 3

#GOES08
#https://www.ospo.noaa.gov/Operations/GOES/calibration/tables/table3_1.htm

IR_TEFF2TEMP_COEFFS = {

    'GOES08': {

            2: [
                (2556.71, -0.618007, 1.001825, -6.021442e-07),
                (2558.62, -0.66864, 1.002221, -1.323758e-06)
            ],

            3: [
                (1481.91, -0.656443, 1.001914 , -9.535221e-07)
            ],
            4: [
                (934.30, -0.519333, 1.002834, -3.005194e-06),
                (935.38, - 0.553431,1.002894, -3.077855e-06)
            ],

            5: [
                (837.06, -0.383077, 1.000856, 6.026892e-07),
                (837.00, -0.351510, 1.000340,1.761416e-06)
            ]
    },

    'GOES09': {

        2: [
            (2555.18,-0.592268,1.00104,-1.882973E-07),
            (2555.18,-0.592268,1.00104,-1.882973E-07)
        ],

        3: [
            (1481.82,-0.559306,1.001602,-1.010812E-06)
        ],
        4: [
            (934.59,-0.525515,1.002411,-2.148433E-06),
            (934.28,-0.532929,1.002616,-2.584012E-06)
        ],

        5: [
            (834.02,-0.317704,1.001058,-2.245684E-07),
            (834.09,-0.346344,1.001261,-6.031501E-07)
        ]
    },

    'GOES10': {

        2: [
            (2552.9845,-0.63343694,1.0013206,-4.2038547e-07),
            (2552.9845,-0.63343694,1.0013206,-4.2038547e-07)
        ],

        3: [
            (1486.2212,-0.66500842,1.0017857,-7.3885254e-07)
        ],
        4: [
            (936.10260,-0.36939128,1.0017466,-1.4981835e-06),
            (935.98981,-0.41013697,1.0020766,-2.1303556e-06)
        ],

        5: [
            (830.88473,-0.32763317,1.0014057,-9.5563444e-07),
            (830.89691,-0.32184480,1.0013828,-9.3581045e-07)
        ]
    },
    'GOES11': {

        2: [
            (2562.07,-0.651377,1.000828,-1.002675e-07),
            (2562.07,-0.651377,1.000828,-1.002675e-07)
        ],

        3: [
            (1481.53,-0.620175,1.002104,-1.171163e-06)
        ],
        4: [
            (931.76,-0.546157,1.003175,-3.656323e-06),
            (931.76,-0.546157,1.003175,-3.656323e-06)
        ],

        5: [
            (833.67,-0.329940,1.000974,5.000439e-08),
            (833.04,-0.307032,1.000903,1.233306e-07)
        ]
    },
    # GOES12 is the only sat to use 2 sides. As a result its value is a list not a dict of channels, where element
    # 0 is side 1 and element 1 is side 2. The block0 can be queried to sheck id side2 is active
    # block0.iscan.side2. If the value is 0 is inactove, a value of 1 menas it is active
    'GOES12': [
                    {
                            2: [
                                (2562.45,-0.727744,1.002131,-1.173898e-06),
                                (2562.45,-0.727744,1.002131,-1.173898e-06)
                            ],

                            3: [
                                (1536.43,-5.278975,1.016476,-7.754348e-06),
                                (1536.95,-5.280110,1.016383,-7.607900e-06)
                            ],
                            4: [
                                (933.21,-0.534982,1.002693,-2.667092e-06),
                                (933.21,-0.534982,1.002693,-2.667092e-06)
                            ],

                            6: [
                                (751.91,-0.177244,1.000138,1.163496e-06)

                            ]
                    },
                    {
                        2: [
                            (2562.45,-0.704652,1.001948,-8.244692e-07),
                            (2562.45,-0.704652,1.001948,-8.244692e-07)
                        ],

                        3: [
                            (1536.43,-5.287978,1.016547,-7.888533e-06),
                            (1536.27,-5.278729,1.016471,-7.810888e-06)
                        ],
                        4: [
                            (933.21,-0.535578,1.002698,-2.677418e-06),
                            (933.21,-0.535578,1.002698,-2.677418e-06)
                        ],

                        6: [
                            (751.77,-0.178664,1.000158,1.121666e-06)

                        ]
                    }
    ],

    'GOES13': {

                2: [
                    (2561.7421,-1.4755462,1.0028656,-5.8203946e-07),
                    (2561.7421,-1.4755462,1.0028656,-5.8203946e-07)
                ],

                3: [
                    (1522.5182,-4.1556932,1.0142082,-8.0255086e-06),
                    (1521.6645,-4.1411143,1.0142255,-8.0755893e-06)
                ],
                4: [
                    (937.23449,-0.52227011,1.0023802,-2.0798856e-06),
                    (937.27498,-0.51783545,1.0023789,-2.1027609e-06)
                ],

                6: [
                    (749.82589,-0.16089410,1.0006896,-3.9853774e-07)

                ]
    },

    #GOES14 has multiple revisions of coeffs, the latest is set here
    'GOES14': {

        2: [
            (2577.3518,-1.5565294,1.0027731,-4.0683469e-07),
            (2577.3518,-1.5565294,1.0027731,-4.0683469e-07)
        ],

        3: [
            (1519.3488,-3.9656363,1.0133248,-7.5834376e-06),
            (1518.5610,-3.9615475,1.0135737,-7.9139638e-06)
        ],
        4: [
            (933.98541,-0.50944128,1.0029288,-3.3213255e-06),
            (934.19579,-0.51352435,1.0027813,-2.9825801e-06)
        ],

        6: [
            (752.88143,-0.16549136,1.0001953,9.0998038e-07),
            (752.82392,-0.16396181,1.0002291,8.1000947e-07),

        ]
    },
    #GOES15 has multiple revisions of coeffs, the latest is set here
    'GOES15': {

        2: [
            (2562.7905,-1.5903939,1.0026700,-3.1926330e-07),
            (2562.7905,-1.5903939,1.0026700,-3.1926330e-07)
        ],

        3: [
            (1521.1988,-3.9614910,1.0132094,-7.4310172e-06),
            (1521.5277,-3.9580619,1.0130975,-7.3039584e-06)
        ],
        4: [
            (935.89417,-0.51798741,1.0025141,-2.3893260e-06),
            (935.78158,-0.51371643,1.0025320,-2.4517009e-06)
        ],

        6: [
            (753.72229,-0.16676232,1.0002673,7.3287547e-07),
            (753.93403,-0.17147513,1.0001237,1.1424349e-06),

        ]
    },

}



def get_ir_teff2tem_coeffs(satellite=None, channel=None,  detector_side=None):
    """
    Retrieves the coefficients needed to convert the effective temp to absolute temp for
    GOES 08-15 IR channels.
    Al sats feature two sets of detectors. However, except GOES12 they use only one set (1)
    So , when the sat is GOES12 you need to suppy also the side


    :param satellite: strm  sat name
    :param channel: str, chn name
    :param detector: active detector number (0 or 1) in doc is called a or b. Some IR channels have only one detector
    and a double resolution in x dimension because of this
    :param detector_side, 1 or 2
    :return:
    """

    assert satellite is not None, f'satellite={satellite} is invalid'
    assert 'GOES' in satellite, f'satellite={satellite} has to contain string GOES'
    assert channel is not None, f'channel={channel} is invalid'
    assert channel in CHANNELS, f'channel={channel} does not exist, Valid values are {list(CHANNELS.keys())}'
    assert channel>1, f'channel={channel} is not of infrared type'
    if satellite == 'GOES12':
        assert detector_side is not None, f'detector_side arg is needed for {satellite}'
        assert detector_side in (1, 2), f'invalid detector_side={detector_side}. Valid values are 1 or 2'

    try:
        sat_no = int(satellite[-2:])
        if not 8<=sat_no<=15:raise ValueError(f'satellite={satellite} is incorrect. Valid number are in range 08-15')

        chn_coeffs = IR_TEFF2TEMP_COEFFS[satellite]
        if isinstance(chn_coeffs, list):
            chn_coeffs = chn_coeffs[detector_side-1]
        names= 'wavenumber', 'b0', 'b1', 'b2'
        vals = chn_coeffs[channel]

        rd = dict()
        for i in range(len(vals)):
            v = vals[i]
            det_name = f'detector_{i}'
            rd[det_name] = dict(zip(names, v))
        return rd

    except Exception:
        raise TypeError(f'satellite={satellite} has to end with a number')


def vis_calibrate(counts=None, channel_detectors=None, detector_lines = None,
                  satellite=None,  calibrate_to='reflectance'):
    assert calibrate_to in ('radiance', 'reflectance'), f'calibrate_to={calibrate_to} is invalid.' \
                                                        f'Valid values are "radiance" or "reflectance"'

    vis_coeffs = VIS_CALIB_COEFFS[satellite]

    slopes, x0, k = vis_coeffs['m'], vis_coeffs['x0'], vis_coeffs['k']


    out_data = np.full_like(counts, dtype='f4', fill_value=np.nan)


    detectors = list(set(detector_lines))

    for detector_no in detectors:
        if detector_no not in channel_detectors:
            detector_lines[detector_lines == detector_no] = 1
            continue
        try:
            m = slopes[detector_no]
        except Exception:
            m = slopes
        #b = -m * x0
        det_lines = np.argwhere(detector_lines == detector_no).squeeze()
        radiance = m*(counts[det_lines, :] - x0)
        #radiance = m*counts[det_lines, :] + b
        if calibrate_to == 'radiance':
            out_data[det_lines, :] = radiance
        else:
            out_data[det_lines, :] = radiance*k
    return out_data

def ir_calibrate(counts=None,  satellite=None,
                 detector_lines=None,
                 channel=None, side=None, calibrate_to='radiance'):

    assert calibrate_to in ('radiance', 'temperature'), f'calibrate_to={calibrate_to} is invalid.'\
    f'Valid values are "radiance" or "temperature"'

    channel_coeffs = get_ir_teff2tem_coeffs(satellite=satellite, channel=channel,
                                       detector_side=side)

    out_data = np.full_like(counts, dtype='f4', fill_value=np.nan)
    m = IR_SLOPEOFFS_COEFFS[channel]['m']
    b = IR_SLOPEOFFS_COEFFS[channel]['b']


    detectors = list(set(detector_lines))

    for detector_no in detectors:
        det_name = f'detector_{detector_no}'
        if not det_name in channel_coeffs:
            detector_lines[detector_lines == detector_no] = 0
            continue

        det_lines = np.argwhere(detector_lines == detector_no).squeeze()
        radiance = (counts[det_lines, :] - b)/m
        detector_coeffs = channel_coeffs[det_name]

        if calibrate_to == 'radiance':
            out_data[det_lines, :] = radiance
        else:
            t_eff = (c2 * detector_coeffs['wavenumber']) / np.log(((c1 * (detector_coeffs['wavenumber'] ** 3)) / radiance) + 1)
            return detector_coeffs['b0'] + (detector_coeffs['b1'] * t_eff) + (detector_coeffs['b2'] * (t_eff ** 2))

    return  out_data

if __name__ == '__main__':

    c = get_ir_teff2tem_coeffs(satellite='GOES13', channel=6)

    print(c)