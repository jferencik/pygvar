from pygvar.reader import parse_gvar
import os
import logging
def extract_imgbox(gvar_file):
    assert os.path.exists(gvar_file), 'gvar_file %s does not exist' % gvar_file
    b0 = parse_gvar(gvar_file=gvar_file)
    return b0.vis_bbox


if __name__ == '__main__':
    logger = logging.getLogger()
    logger.setLevel('DEBUG')
    gvar_file = '/work/data/goese_hybrid/OPERATIONAL/CLASS/INCOMING/092/2817777983.goes13.2017.092.034520'
    print extract_imgbox(gvar_file=gvar_file)
