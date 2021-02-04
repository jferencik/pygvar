"""
Some background is necessary to understand the navigation concepts in GOES pre 16 and after 7 (GVAR AAA)

All GOES instruments (after 7) transmit data in GVAR format.
GVAR is a block based format that streams data from the sensor to the ground stations and is then assembled and disseminated.
Because it is a stream based format, block0 that carries metadata (Q&A, navigation, etc) is is send in a repeated manner, every tents block is block 0
as the instumernts progresses in scan swatghs from north to south and west to east

Because GVAR is indeed complex (holds IMAGER and SOUNDER) a simpler format called AREA was developed.

An AREA file is a collection of information that defines the image and its associated ancillary data.
Area files consist of these six blocks:

Directory block
Navigation block (NAV block)
Calibration block (CAL block)
Supplemental block
Data block
Comment block (AUDIT block)


The Navigation block holds info required to perform  navigation of GOEs images, that is, tranformation from line|row/columns to latlon and reverse


These tranformations have been developed historically as Fortran code, incorporated into McIDAS X, then ported to Java by SSEC and finally ported to python
by Solargis. The porting was done in the comntext of Area file librray, therefore, the navigatin is implemented/folows java code as a python superclass (AreaNav)
 and specific subclasses that correspond to various instruments/sensors.
 
 
As pygvar features a modular/procedural design, it is necessary to refactor the navigation/tranformation code from a class based structure to a function based one.
It is important to note that because of historical reasons, the navigation classes are murky, obscure and object oriented. 
This refactoring should simplify them. To be specific, for example,  the navigation info is extracted from the Area file navigation block, converted to integer
 muktiplied/divided with various constrants, packed into bunary and the unpacked and converted back inside the navigation classes.
 All this is redundant, pygvar attempts to just read the info from block0 and then perform the navigation.




The goal of this module is to test the Area Navigation vs the same area navigation implemented as a fucntion inside the pygvar module.

Several steps are tested

1. equivalence of navigation, that is, given a soecuific scan time, teh same image should exist in Area and GVAr format.
The test shall evaluate if the vanigation bytes in the form the Area Navigatin class exects them are equiavelent


"""

import unittest
import os
import sys
import datetime
from general_utils import daytimeconv



class Testpygvar(unittest.TestCase):
    goese_folder = '/work/dev/python/goese/src'
    def test_navigation_bytes(self):
        """
        For this function to work, it is necessary to have the Area file python library.
        Because of Solargis dev standard and style, make sure the AreaFile lib is the one form goese package.
        
        The files shall be picked from CLASS/CYCLONE for the same time.
        
       
        """
        #first update path so goese can be imported

        sys.path.append(self.goese_folder)
        from goese.area.imager.file import AreaFile
        from goese.area.imager.navigation import AreaNav
        from goese.util import decompresFile
        from pylab import imshow, show
        import numpy as np
        from pygvar.reader import parse_gvar
        from pygvar import navigation as gvar_nav


        gvar_file = '/work/data/goese_hybrid/OPERATIONAL/TEST/GVAR/2824983383.goes13.2017.121.001519'
        gvar_file = '/work/data/goese_hybrid/OPERATIONAL/TEST/GVAR/2825044283.goes13.2017.121.061518'
        gvar_file = '/work/data/goese_hybrid/OPERATIONAL/TEST/GVAR/2825197333.goes13.2017.121.201519'
        area_file = '/work/data/goese_hybrid/OPERATIONAL/TEST/AREA/G13-201705010015.gstn.vis.enhr.ara.Z'
        area_file = '/work/data/goese_hybrid/OPERATIONAL/TEST/AREA/G13-201705012015.gstn.vis.enhr.ara.Z'

        #b0, gvar_vis_data = parse_gvar(gvar_file=gvar_file,channels_to_extract=[1])
        b0 = parse_gvar(gvar_file=gvar_file)
        print b0.tched.python_value()
        print b0.subla.python_value()
        print b0.sublo.python_value()
        dec_area_file = decompresFile(area_file)

        af = AreaFile(dec_area_file)
        an = af.navigation

        lat = 6, 7
        lon = -80, -82 #McIDAS uses positive longitude in easter hemisphere
        #use the simple original fucntion to test the transformations
        print lat, lon
        # ele, lin =  an.toLinEle(latlon=[[lat], [lon]])
        # print lin[0], ele[0]
        # clon, clat = an.toLatLon(linele=[ele, lin])
        # print clon[0], clat[0]

        #tehe are numpy based functions
        lats = np.array([lat])
        lons = np.array([lon])
        lin, ele = an.toLinEle2(lats=lats, lons=lons)
        print lin, ele
        clat, clon =  an.toLatLon2(lines=lin, cols=ele)
        print clat, clon

        #data = af.data.squeeze()

        # imshow(gvar_vis_data[1])
        # show()
        # imshow(data)
        # show()

        w3 = af.nav
        #print len(w3)


        #n1 = gvar_nav.prepare_mcidas_nav1(b0=b0)

        # navigation = AreaNav(navblock=n1, byte_order='<', int_format='i')
        # navigation.setImageStart(b0.infln, b0.iwfpx)
        # navigation.setRes(1, 1)
        # navigation.setMag(1, 1)
        # navigation.setStart(0, 0)
        # e1, l1 = navigation.toLinEle(latlon=[[lat], [lon]])
        # print l1, e1
        olin, oele =  gvar_nav.ll2le(lats=lats, lons=lons, b0=b0)
        print olin, oele
        olat, olon = gvar_nav.le2ll(lines=olin, cols=oele, b0=b0)
        print olat, olon

        af.close()
        os.remove(dec_area_file)




if __name__ == '__main__':

    unittest.main()
