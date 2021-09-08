from pygvar.reader import parse_gvar
import glob
import numpy as np
PI = np.pi
DEG = 180. / PI
RAD = PI / 180.  # degrees to radians conversion pi/180

def view_gvar(gvar_file, channel_to_extract=None):
    b0, data = parse_gvar(gvar_file, channels_to_extract=[channel_to_extract], fill_missing_lines=True)
    nl, nc, linres, eleres = b0.channels_shape[channel_to_extract]
    # print (nl, nc, linres, eleres)
    # print (b0.vis_bbox)
    # print (b0.g1cnt, b0.g2cnt)
    # #print b0.iscan
    #
    # lam = b0.reflo.python_value()
    #
    #
    # phi = b0.refla.python_value()
    #
    #
    # subpoint = []
    # subpoint.append(phi / RAD)
    # subpoint.append(lam / RAD)
    #
    #
    # la, lo = b0.subla.python_value(), b0.sublo.python_value()
    # print ('original LAT %f LON %f' % (la, lo))
    # # la, lo = subpoint[0], subpoint[1]
    # # print 'original LAT %f LON %f' % (la, lo)
    #
    # from pygvar.navigation import ll2le, le2ll
    # l, e = ll2le(la,lo,b0=b0, lineres=linres, colres=eleres)
    # bla, blo = le2ll(l,e,b0, lineres=linres, colres=eleres)
    # print ('transformed back LAT %f LON %f' % (bla, blo))
    #
    # edge_lat, edge_lon = 27.848373, -115.081904
    # # edge_lat, edge_lon = 20.890275, -157.061706
    # # edge_lat, edge_lon = 46.726106, -92.161731
    # edge_lat, edge_lon = 19.388956, -155.105946
    # edge_lin, edge_ele = ll2le(lats=edge_lat, lons=edge_lon, b0=b0, lineres=linres, colres=eleres)
    # print(int(edge_lin), int(edge_ele), float(nc)/int(edge_ele))
    #
    # # this perhaps need to be specifically checked for
    # imc = 0 if (b0.iscan.imc_active != 0) else  1
    # # imc = b0.iscan.imc_active  # should suffice because it is either 0 or 1
    # iflip = -1 if b0.iscan.yaw_flip else 1
    # #print imc, iflip
    chn_data = data[channel_to_extract]
    # data[int(edge_lin), int(edge_ele)] = 1022
    try:
        import pylab
        pylab.imshow(chn_data, interpolation='nearest')
        pylab.show()

    except ImportError:
        print('failed to import pylab')




if __name__ == '__main__':
    import logging
    import re
    import os
    logging.basicConfig()
#    import pylab
    logger = logging.getLogger()
    logger.setLevel('INFO')

    #print src_file
    src_file = '/data/sample_gvar/goes12.2003.152.223144'

    view_gvar(src_file, 2)
    #b0, datad = parse_gvar(gvar_file=src_file, channels_to_extract=[chn])

    # data = datad[chn]
    # #data[602:data.shape[0]-8,:] = data[610:,:]
    # pylab.imshow(data, cmap='jet', interpolation='nearest')
    # pylab.title(name)
    # pylab.show()






