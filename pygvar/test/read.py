from pygvar.reader import parse_gvar
import glob
import numpy as np
PI = np.pi
DEG = 180. / PI
RAD = PI / 180.  # degrees to radians conversion pi/180

def view_gvar(gvar_file, channel_to_extract=None):
    b0, data = parse_gvar(gvar_file, channels_to_extract=[channel_to_extract], fill_missing_lines=True)

    chn_data = data[channel_to_extract]

    try:
        import pylab
        pylab.imshow(chn_data, interpolation='nearest')
        pylab.show()

    except ImportError:
        print('failed to import pylab')




if __name__ == '__main__':
    import logging
    logging.basicConfig()
#    import pylab
    logger = logging.getLogger()
    logger.setLevel('INFO')

    src_file = '/data/sample_gvar/goes12.2003.152.223144'

    view_gvar(src_file, 1)






