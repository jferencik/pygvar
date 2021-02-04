import netCDF4 as nc
from osgeo import ogr, osr
import os



def test(subsatellite_point):
    geos_proj_str = '+proj=geos  +h=35785831.0 +lon_0=%.2f +units=m +ellps=WGS84' % float(subsatellite_point)
    geos_proj = osr.SpatialReference()
    geos_proj.ImportFromProj4(geos_proj_str)

    center = ogr.Geometry(ogr.wkbPoint)
    center.AddPoint(0,0)
    #circle = center.Buffer(5500000, 1000)
    circle = center.Buffer(5500000-82000, 3000)

    driver = ogr.GetDriverByName('ESRI Shapefile')
    ds = driver.CreateDataSource('/tmp/geos_circle.shp')
    if ds.GetLayerByName('circ') is not None:
        ds.DeleteLayer('circ')
    lyr = ds.CreateLayer('circ', geos_proj, ogr.wkbPolygon)
    layerDef = lyr.GetLayerDefn()
    feature = ogr.Feature(layerDef)
    feature.SetGeometry(circle)
    feature.SetFID(0)
    lyr.CreateFeature(feature)

def test1(ncpath):
    geos_proj_str = '+proj=geos  +h=35785831.0 +lon_0=%.2f +units=m +ellps=WGS84' % float(-75)
    _, n = os.path.split(ncpath)

    d = nc.Dataset(ncpath)

    v = d.variables['Rad']
    v.set_auto_maskandscale(False)
    #v.set_auto_maskandscale(False)

    data = v[:].squeeze()
    print(type(data)), data.dtype
    #data = data>>5


    from pygvar.utils import toGDAL
    toGDAL(fpath='/tmp/%s.tif' % n,array=data,gt=[-5500000., 2000., 0., 5500000., 0., -2000.],proj=geos_proj_str)








if __name__ == '__main__':
    from pygvar import reader
    import logging
    import pylab
    from pygvar.navigation import ll2le
    logging.basicConfig()
    logger = logging.getLogger()
    logger.setLevel('DEBUG')
    gvar_path = '/media/d/data/goes08/goes08.1999.019.171515'

    # b0, data = reader.parse_gvar(gvar_file=gvar_path, channels_to_extract=[1])
    # chn_data = data[1]
    b0= reader.parse_gvar(gvar_file=gvar_path, )
    lon0 = b0.sublo.python_value()
    test(lon0)
    #test1('/media/d/data/satellite_goes_r/http/incoming/OR_ABI-L1b-RadF-M3C07_G16_s20181690600431_e20181690611209_c20181690611243.nc')

    from mtsat.hrit import geometry_correction
    print geometry_correction.LAKES_SHAPEFILE
    #geometry_correction.compute_mask()
    # pylab.imshow(chn_data, interpolation='nearest')
    # pylab.show()
