from pygvar.reader import parse_gvar
from pygvar.calibration import vis_calibrate, ir_calibrate
import os
import numpy as np

def view_gvar(gvar_file, channel=None, output=None):
    assert output in (None, 'radiance', 'reflectance', 'temperature')
    satellite = os.path.basename(gvar_file).split('.')[0].upper()

    b0, data = parse_gvar(gvar_file, channels_to_extract=[channel], fill_missing_lines=True)

    chn_data = data[channel]

    if output is not None:
        if channel ==  1:
            calib_data = vis_calibrate(counts=chn_data,channel_detectors=b0.detectors[channel],
                                       detector_lines=b0.detector_lines[channel],
                                     satellite=satellite, calibrate_to=output)
        else:
            side = 2 if b0.iscan.side2 == 1 else 1
            calib_data = ir_calibrate(counts=chn_data, detector_lines = b0.detector_lines[channel],
                                      satellite=satellite, channel=channel, side=side, calibrate_to=output)
    else:
        calib_data = chn_data
    try:
        import pylab
    except ImportError:
        print('failed to import pylab')
    print(np.nanmin(calib_data), np.nanmax(calib_data))
    pylab.imshow(calib_data, interpolation='nearest')
    pylab.show()


if __name__ == '__main__':
    import logging
    logging.basicConfig()
#    import pylab
    logger = logging.getLogger()
    logger.setLevel('INFO')

    src_file = '/data/sample_gvar/goes12.2003.152.223144'

    view_gvar(src_file, channel=4, output='temperature')






