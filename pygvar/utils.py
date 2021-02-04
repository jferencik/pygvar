import logging
import numpy as np
import os
import struct
import time


_ , mname = os.path.split(__file__)
logger = logging.getLogger(mname)



crc16tab = [
    0x0000,0x1021,0x2042,0x3063,0x4084,0x50a5,0x60c6,0x70e7,
    0x8108,0x9129,0xa14a,0xb16b,0xc18c,0xd1ad,0xe1ce,0xf1ef,
    0x1231,0x0210,0x3273,0x2252,0x52b5,0x4294,0x72f7,0x62d6,
    0x9339,0x8318,0xb37b,0xa35a,0xd3bd,0xc39c,0xf3ff,0xe3de,
    0x2462,0x3443,0x0420,0x1401,0x64e6,0x74c7,0x44a4,0x5485,
    0xa56a,0xb54b,0x8528,0x9509,0xe5ee,0xf5cf,0xc5ac,0xd58d,
    0x3653,0x2672,0x1611,0x0630,0x76d7,0x66f6,0x5695,0x46b4,
    0xb75b,0xa77a,0x9719,0x8738,0xf7df,0xe7fe,0xd79d,0xc7bc,
    0x48c4,0x58e5,0x6886,0x78a7,0x0840,0x1861,0x2802,0x3823,
    0xc9cc,0xd9ed,0xe98e,0xf9af,0x8948,0x9969,0xa90a,0xb92b,
    0x5af5,0x4ad4,0x7ab7,0x6a96,0x1a71,0x0a50,0x3a33,0x2a12,
    0xdbfd,0xcbdc,0xfbbf,0xeb9e,0x9b79,0x8b58,0xbb3b,0xab1a,
    0x6ca6,0x7c87,0x4ce4,0x5cc5,0x2c22,0x3c03,0x0c60,0x1c41,
    0xedae,0xfd8f,0xcdec,0xddcd,0xad2a,0xbd0b,0x8d68,0x9d49,
    0x7e97,0x6eb6,0x5ed5,0x4ef4,0x3e13,0x2e32,0x1e51,0x0e70,
    0xff9f,0xefbe,0xdfdd,0xcffc,0xbf1b,0xaf3a,0x9f59,0x8f78,
    0x9188,0x81a9,0xb1ca,0xa1eb,0xd10c,0xc12d,0xf14e,0xe16f,
    0x1080,0x00a1,0x30c2,0x20e3,0x5004,0x4025,0x7046,0x6067,
    0x83b9,0x9398,0xa3fb,0xb3da,0xc33d,0xd31c,0xe37f,0xf35e,
    0x02b1,0x1290,0x22f3,0x32d2,0x4235,0x5214,0x6277,0x7256,
    0xb5ea,0xa5cb,0x95a8,0x8589,0xf56e,0xe54f,0xd52c,0xc50d,
    0x34e2,0x24c3,0x14a0,0x0481,0x7466,0x6447,0x5424,0x4405,
    0xa7db,0xb7fa,0x8799,0x97b8,0xe75f,0xf77e,0xc71d,0xd73c,
    0x26d3,0x36f2,0x0691,0x16b0,0x6657,0x7676,0x4615,0x5634,
    0xd94c,0xc96d,0xf90e,0xe92f,0x99c8,0x89e9,0xb98a,0xa9ab,
    0x5844,0x4865,0x7806,0x6827,0x18c0,0x08e1,0x3882,0x28a3,
    0xcb7d,0xdb5c,0xeb3f,0xfb1e,0x8bf9,0x9bd8,0xabbb,0xbb9a,
    0x4a75,0x5a54,0x6a37,0x7a16,0x0af1,0x1ad0,0x2ab3,0x3a92,
    0xfd2e,0xed0f,0xdd6c,0xcd4d,0xbdaa,0xad8b,0x9de8,0x8dc9,
    0x7c26,0x6c07,0x5c64,0x4c45,0x3ca2,0x2c83,0x1ce0,0x0cc1,
    0xef1f,0xff3e,0xcf5d,0xdf7c,0xaf9b,0xbfba,0x8fd9,0x9ff8,
    0x6e17,0x7e36,0x4e55,0x5e74,0x2e93,0x3eb2,0x0ed1,0x1ef0
]
crct_array = np.array(crc16tab)

BIT_SHIFTS = [i for i in range(6, -2, -2)]
TEN_BIT_MASK = [1023 << i for i in BIT_SHIFTS]

def crc16(data=None, crc=0xffff):
    """Calculate CRC16 using the existing table.
    `data`      - data for calculating CRC, must be a string
    `crc`       - initial value
    `table`     - table for caclulating CRC (list of 256 integers)
    Return calculated value of CRC
    """
    idata = struct.unpack('%dB' % len(data), data)

    for byte in idata:
        crc = ((crc<<8)&0xff00) ^ crct_array.take(((crc>>8)&0xff)^byte)
    return 65535-crc

def f2int(val):
    """
    Traditional statistical rounding
    :param val:
    :return:
    """

    return int(round(val+10**(-len(str(val))-1)))

def isiter(a):
    try:
        (x for x in a)
        return True
    except TypeError:
        return False

def eight_bit2ten_bit(in_array=None):

    """
        Converts a int eight bit byte array into a 10 bit byte array. Numpy optimised  version

        :param, in_array, numpy array of int type
        :return, numpy array of ype u2 where only 10 bits are used


        TEN_BIT_MASK is the mask used to convert an array of u1 or 8 bit ints to an array of 10 bit ints as follows
        given the  binary string  consisting of 5 8 bit words
        convert it to a 4 10 bit words without using C

        The issue is python per se does not "know" about 10bit integers.
        So what I am doing here is to concatenate every 2 consecutive bytes|elements into a short and then  bitwise and & with the correponding mask and then right shift the necessary number of bits to get the final 10bit int

        theoretical example





        1. binary string
        s = '\xc77e\xf79'
        2. convert to 1 byte
        sb = np.frombuffer(s, dtype='u1')
        print [199  55 101 247  57]
        #bit  representation
        ['{:08b}'.format(e) for e in sb]
        ['11000111', '00110111', '01100101', '11110111', '00111001']

        3. shifts
        shifts = [i for i in range(6, -2, -2)]
        print shifts
        [6, 4, 2, 0]

        3. create masks,
         TEN_BIT_MASK = [1023 << i for i shifts]

        #mask representation
        ['{:08b}'.format(e) for e in TEN_BIT_MASK]
        ['1111111111000000', '0011111111110000', '0000111111111100', '0000001111111111']

        #convert 5 8 bit words (40 bits) into 4 10bit words 40 bits

        ten_bit = np.zeros(4, dtype='u2') # 16 bits ints
        ten_bit[::4]=(((sb[::5]*2**8) | sb[1::5])&TEN_BIT_MASK[0])//2**shifts[0]
        ten_bit[1::4]=(((sb[1::5]*2**8) | sb[2::5])&TEN_BIT_MASK[1])//2**shifts[1]
        ten_bit[2::4]=(((sb[2::5]*2**8) | sb[3::5])&TEN_BIT_MASK[2])//2**shifts[2]
        ten_bit[3::4]=(((sb[3::5]*2**8) | sb[4::5])&TEN_BIT_MASK[3])//2**shifts[3]

        or in a loop
        ten_bit = np.zeros(4, dtype='u2') # 16 bits ints
        for i in range(4):
             ten_bit[i::4]=(((sb[i::5]*2**8) | sb[i+1::5])&TEN_BIT_MASK[i])//2**shifts[i]


        using the loop it is more efficient


    """
    #we need a numpy array of
    #handle byte order
    assert hasattr(in_array, 'dtype'), 'A numpy array is required'
    assert in_array.size > 0, 'The input array is empty'
    asize = in_array.size
    #assert in_array.size % 5 == 0, 'The input array has to have multiple 5 number of elements'

    five_rem = asize % 5
    n_out_elems = (asize // 5) * 4
    if five_rem != 0:
        in_array = in_array[:asize-five_rem]
    out_array = np.zeros(n_out_elems, dtype='>u2') # has to stay like this
    for i in range(4):
        out_array[i::4] = (((in_array[i::5] * 2 ** 8) | in_array[i + 1::5]) & TEN_BIT_MASK[i]) // 2 ** BIT_SHIFTS[i]
    return  out_array




def toGDAL(fpath=None,array=None, fmt='GTiff', gt=(0,1,0,0,0,-1), no_data=-9999, proj=None):
    """
        Writes a numpy array to a file format suported by GDAL. By default GTiff is used.
        This function is flexible, it has two modus operandi.
        @args:
            @folder, the forlder where the data will be written. Make sure it is writable and exists.
            @name, the name of the file
            @fmt, GDAL driver, defaults to HFA, nissue here with extension as GDAL can write an Imagine file with a different extension.
            @gt, GDAL geotransform ,In case of north up images, the gt(2) and gt(4) coefficients are zero, and the gt(1) is pixel width, and gt(5) is pixel height. The (gt(0),gt(3)) position is the top left corner of the top left pixel of the raster.
            @no_data, the values to be used as NODA    all_files = sorted(os.listdir(incoming_area_folder))
            @epsg_code, int, the projection code
        @returns:
            None

    """
    from osgeo import gdal, osr, gdalnumeric
    driver = gdal.GetDriverByName(fmt)
    ds = driver.Create( fpath, array.shape[1], array.shape[0], 1, gdalnumeric.NumericTypeCodeToGDALTypeCode(array.dtype.type) )
    band = ds.GetRasterBand(1)
    band.WriteArray(array)
    if no_data is not None:
        band.SetNoDataValue(no_data)
    ds.SetGeoTransform(gt)
    prj = osr.SpatialReference() # custom spatial reference
    try:
        prj.ImportFromEPSG(int(proj))
    except ValueError:
        prj.ImportFromProj4(proj)
    except TypeError:
        return


    ds.SetProjection(prj.ExportToWkt())


def flush_cache(passwd):
    """
    Flush Linux VM caches. Useful for doing meaningful tmei measurements for
    NetCDF or similar libs.
    Needs sudo password
    :return: bool, True if success, False otherwise
    """
    logger.debug('Clearing the OS cache using sudo -S sh -c "sync; echo 3 > /proc/sys/vm/drop_caches')
    #ret = os.system('echo %s | sudo -S sh -c "sync; echo 3 > /proc/sys/vm/drop_caches"' % passwd)
    ret = os.popen('sudo -S sh -c "sync; echo 3 > /proc/sys/vm/drop_caches"', 'w').write(passwd)
    return not bool(ret)

def timeit(func=None,loops=1,verbose=False, clear_cache=False, sudo_passwd=None):
    #print 0, func, loops, verbose, clear_cache, sudo_passwd
    if func != None:
        if clear_cache:
            assert sudo_passwd, 'sudo_password argument is needed to clear the kernel cache'

        def inner(*args,**kwargs):
            sums = 0.0
            mins = 1.7976931348623157e+308
            maxs = 0.0
            logger.debug('====%s Timing====' % func.__name__)
            for i in range(0,loops):
                if clear_cache:
                    flush_cache(sudo_passwd)
                t0 = time.time()
                result = func(*args,**kwargs)
                dt = time.time() - t0
                mins = dt if dt < mins else mins
                maxs = dt if dt > maxs else maxs
                sums += dt
                if verbose == True:
                    logger.debug('\t%r ran in %2.9f sec on run %s' %(func.__name__,dt,i))
            logger.debug('%r min run time was %2.9f sec' % (func.__name__,mins))
            logger.debug('%r max run time was %2.9f sec' % (func.__name__,maxs))
            logger.info('%r avg run time was %2.9f sec in %s runs' % (func.__name__,sums/loops,loops))
            logger.debug('==== end ====')
            return result

        return inner
    else:
        def partial_inner(func):
            return timeit(func,loops,verbose, clear_cache, sudo_passwd)
        return partial_inner

def mode(ndarray, axis=0):
    # Check inputs
    ndarray = np.asarray(ndarray)
    ndim = ndarray.ndim
    if ndarray.size == 1:
        return (ndarray[0], 1)
    elif ndarray.size == 0:
        raise Exception('Cannot compute mode on empty array')
    try:
        axis = range(ndarray.ndim)[axis]
    except:
        raise Exception('Axis "{}" incompatible with the {}-dimension array'.format(axis, ndim))

    # If array is 1-D and numpy version is > 1.9 numpy.unique will suffice
    if all([ndim == 1,
            int(np.__version__.split('.')[0]) >= 1,
            int(np.__version__.split('.')[1]) >= 9]):
        modals, counts = np.unique(ndarray, return_counts=True)
        index = np.argmax(counts)
        return modals[index], counts[index]

    # Sort array
    sort = np.sort(ndarray, axis=axis)
    # Create array to transpose along the axis and get padding shape
    transpose = np.roll(np.arange(ndim)[::-1], axis)
    shape = list(sort.shape)
    shape[axis] = 1
    # Create a boolean array along strides of unique values
    strides = numpy.concatenate([np.zeros(shape=shape, dtype='bool'),
                                 np.diff(sort, axis=axis) == 0,
                                 np.zeros(shape=shape, dtype='bool')],
                                axis=axis).transpose(transpose).ravel()
    # Count the stride lengths
    counts = np.cumsum(strides)
    counts[~strides] = np.concatenate([[0], numpy.diff(counts[~strides])])
    counts[strides] = 0
    # Get shape of padded counts and slice to return to the original shape
    shape = np.array(sort.shape)
    shape[axis] += 1
    shape = shape[transpose]
    slices = [slice(None)] * ndim
    slices[axis] = slice(1, None)
    # Reshape and compute final counts
    counts = counts.reshape(shape).transpose(transpose)[slices] + 1

    # Find maximum counts and return modals/counts
    slices = [slice(None, i) for i in sort.shape]
    del slices[axis]
    index = np.ogrid[slices]
    index.insert(axis, np.argmax(counts, axis=axis))
    return sort[index], counts[index]