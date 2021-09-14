import logging
import os
import datetime
from pygvar.utils import  crc16, eight_bit2ten_bit, toGDAL, mode
import numpy as np
from pygvar.classes import block0
from binascii import crc_hqx
import math


PI = math.pi
DEG = 180. / PI
RAD = PI / 180.  # degrees to radians conversion pi/180
_ , mname = os.path.split(__file__)
logger = logging.getLogger(__name__)




#************** START HARD CODED CONSTANTS AS PER OGE ************************************

MAGIC_CRC =  0x1d0f #
N_HEADER_BYTES = 30 # one header = 30 bytes
N_LDOC_BYTES = 16 # words|bytes in line doc
MAX_MISSING_IRLINES_TO_INTERPOLATE = 4 # interpolate missing lines if their continous sequence is less then or equal ot this value
MAX_MISSING_VISLINES_TO_INTERPOLATE = 8 # interpolate missing lines if their continous sequence is less then or equal ot this value
MIN_USABLE_IMAGE_PERCENTAGE = 5 # if image is corrupted and at leat this percentage of image was read return the image and not raise error

BO = '>'
HEADER = [('BlockID', 'B'),
          ('WordSize', 'B'),
          ('WordCount', 'H' ),
          ('ProductID', 'H' ),
          ('RepeatFlag', 'B' ),
          ('VersionNumber', 'B' ),
          ('DataValid', 'B'),
          ('AsciiBin', 'B' ),
          ('SPS_ID', 'B'),
          ('RangeWord', 'B'),
          ('BlockCount', 'H'),
          ('Spares_1', '2B'),
          ('SPS_Time', '8B'),
          ('Spares_2', '4B'),
          ('CRC', 'H')]

HEADER = [(e[0], BO + e[1]) for e in HEADER]

LDOC = ('SPCID', 'SPSID', 'LSIDE', 'LIDET', 'LICHAA', 'RISCT', 'L1SCAN', 'L2SCAN', 'LPIXLS', 'LWORDS', 'LZCORR', 'LLAG', 'LSPAR')


#************** END HARD CODED CONSTANTS AS PER OGE ************************************

def check_header_crc(buff_str, crc_value=None):
    """
    Given a binary string representing a GVAR imager block header minus  the last two words, computes the crc of the header
    to establish if the header was compromised by noise

    :param buff_str: input binary str
    :param crc_value: int the CRC value storse in the last 2 words of the header interpreted as int16 integer
    :return: bool
    """
    return crc_value == crc16(buff_str)

def wrap_ldoc(ldoc_bytes):
    """
    Creates a dict representing a line doc object. A line doc object consists of first 160 bits
    of a line

    :param ldoc_bytes:
    :return:
    """

    out = np.zeros(len(LDOC), dtype='u2')
    out[:] = -1
    out[0] = ldoc_bytes[0]
    out[1] = ldoc_bytes[1]
    out[2] = ldoc_bytes[2]
    out[3] = ldoc_bytes[3]
    out[4] = ldoc_bytes[4]
    out[5] = ldoc_bytes[5] << 10 | ldoc_bytes[6]
    out[6] = ldoc_bytes[7]
    out[7] = ldoc_bytes[8]
    out[8] = ldoc_bytes[9] << 10 | ldoc_bytes[10]
    out[9] = ldoc_bytes[11] << 10 | ldoc_bytes[12]
    out[10] = ldoc_bytes[13] << 10 | ldoc_bytes[14]
    out[11] = ldoc_bytes[15]

    return dict(zip(LDOC, out))


def parity(v):
    par = 0xff
    for e in v:
        par ^= e
    return par
def b0_parity_check(b0data):
    """
    The deal here is that any block 0 has 6 sections
     we do the parity check for section one, where the most important info exists
     Theoretically, the other block should be checked as well but in reality
     I found out that many times they are not used, so it does not really matter if they are corrupted
    :param b0data:
    :return:
    """
    #these are the offsets intoa block 0 as per OGE
    a = 0,  278, 1626, 2306, 5386, 6304
    b =   277, 1625, 2305, 5385, 6303, 8039
    return parity(b0data[a[0]:b[0]]) == b0data[b[0]]

def header2dict(h):
    return dict(zip(h.dtype.names, *h))

def find_header(fmem=None, fmem_offset=0, nhdrs=5):
    """
    Finds a triplet or the best approximation of 5 consecutive block headers in a byte 1D array
    starting at specific offset

    If this function  does not return a header the data is really bad! ...

    The fucntion works in following way:

    1. tke the first 90 bytes from the offset and tryc to convert them to three identical headers.
    2. there can be 3 good headers (ideal), two goo headers, one good header, or not even one good heasder

    3. if there is at least 1 good header this are still OK
    4. if there is none, the function will advance 1 byte at a time for up to 60 bytes and try to find the headers.
    5 if this failed the function return to the begining (offset) and collapses the three headers into one using majority vote, because the chance that all three
    headers have the same element (positional) corrpted is really small.
    6 I step five does not produce a valid header


    :param fmem: 1D numpy array (byte) usually  resulted by mem mapping the GVAR imager file
    :param fmem_offset: int, the offste where to start searching
    :param nhdrs: int, defaults to 5
    :return:
    """
    #logger.debug('searching for new header')
    # an idea, test if fmemoffset is safe here
    #print fmem_offset
    block_offset = None#
    hindex = None
    tracker = 0
    local_mem_ffset = fmem_offset

    while True:
        #read using numpy
        #some sanity check,  just in case

        if len(fmem) -  fmem_offset < nhdrs*N_HEADER_BYTES:
            raise IOError('Not enough bytes %d in file to find headers ' % (nhdrs*N_HEADER_BYTES))
        headers = np.frombuffer(fmem, dtype=HEADER,count=nhdrs, offset=local_mem_ffset)
        good_headers = [np.frombuffer(h_.data, dtype=HEADER) for h_ in headers]

        #compute crcs, the bytes object is necessary to ensure compatibility between python2 and 3
        #as the data atrr differs very much between python2 and three
        hcrc = np.array([check_header_crc(bytes(h.data)[:-2], h['CRC']).item() for h in headers])
        #number of good headers
        n_good_headers = hcrc[hcrc].size
        # logger.debug(hcrc)
        #handle 3 or more
        if n_good_headers >= 3:


            if hcrc[0] and hcrc[1] and hcrc[2]: #normal case
                block_offset =  0
                hindex = 0
            elif hcrc[1] and hcrc[2] and hcrc[3]: #middle
                block_offset = 1
                hindex = 1
            elif hcrc[2] and hcrc[3] and hcrc[4]:
                block_offset = 2
                hindex=2

        elif n_good_headers == 2: # handle 2

            if (hcrc[0] and hcrc[1]) or (hcrc[0] and hcrc[2]):

                block_offset = 0
                hindex = 0
            elif hcrc[1] and hcrc[3]:
                block_offset = 1
                hindex=1
            elif hcrc[2] and hcrc[4]:
                block_offset = 2
                hindex=2
            elif  hcrc[1] and hcrc[2]:

                #is [0,1,2] or [1,2,3] gh = gheaders[1]
                d1 = good_headers[1].view('u1') == good_headers[0].view('u1')
                d2 = good_headers[1].view('u1') == good_headers[3].view('u1')

                if d1[d1].size > d2[d2].size:
                    block_offset = 0
                    hindex = 1
                else:
                    block_offset = 1
                    hindex = 1
            elif hcrc[2] and hcrc[3]:
                # is [1,2,3] or [2,3,4]
                d1 = good_headers[2].view('u1') == good_headers[1].view('u1')
                d2 = good_headers[2].view('u1') == good_headers[4].view('u1')
                if d1[d1].size > d2[d2].size:
                    block_offset = 1
                    hindex = 2
                else:
                    block_offset = 2
                    hindex = 2

        elif n_good_headers == 1: # only one good header

            if hcrc[0]:
                #[0,1,2] is the only possibility
                d1 = good_headers[0].view('u1') == good_headers[1].view('u1')
                d2 = good_headers[0].view('u1') == good_headers[2].view('u1')
                if d1[d1].size>=20 and d2[d2].size >= 20: # should be valid headers, or is better than and, however it is weak
                    block_offset = 0
                    hindex = 0
                #this is improvement over gvar tool
                elif d1[d1].size>=20 and d2[d2].size <= 20:

                    #print fmem[fmem_offset-N_HEADER_BYTES:fmem_offset]
                    block_offset = -1
                    hindex = 0

                elif d2[d2].size >= 20 and d1[d1].size <= 20:
                    block_offset = 0
                    hindex = 0
            elif hcrc[1]:#[False  True False False False]
                #possibilities are [0,1,2] or [1,2,3]
                d1 = good_headers[1].view('u1') == good_headers[0].view('u1')
                d2 = good_headers[1].view('u1') == good_headers[3].view('u1')
                if d1[d1].size >= d2[d2].size: #  first is better then third
                    block_offset = 0
                    hindex = 1

                else:
                    block_offset = 1
                    hindex = 1


            elif hcrc[2]: #the most complicated case
                #[0,1,2], [1,2,3], [2,3,4]
                gh = good_headers[2].view('u1')
                gh_0 = gh == good_headers[0].view('u1')
                gh_1 = gh == good_headers[1].view('u1')
                gh_3 = gh == good_headers[3].view('u1')
                gh_4 = gh == good_headers[4].view('u1')

                a = gh_0[gh_0].size + gh_1[gh_1].size
                b = gh_1[gh_1].size + gh_3[gh_3].size
                c = gh_3[gh_3].size + gh_4[gh_4].size

                if  a >= b and a >=c: #[0,1,2]
                    block_offset = 0
                    hindex=2
                elif  b >= a and b >=c:#[1,2,3]
                    block_offset = 1
                    hindex = 2
                else: #[2,3,4]
                    block_offset = 2
                    hindex = 2



        if block_offset is None:
            #advance one byte forward and do the whole search again

            tracker+=1

            #logger.debug('tracker %d' % tracker)

            if tracker<=nhdrs*N_HEADER_BYTES-1:
                #logger.debug(headers['BlockCount'])
                #logger.debug('Moving forward one byte to offset %d' % tracker)
                local_mem_ffset+=tracker
                continue
            else:
                #all options have been exhausted, try majority vote on headers, there is noting we loose. If it works we continue if not.... maktub
                block_offset = 0
                tracker = 0
                #create only three headers this time, they are not really neded , we just do it to be consistent
                headers = np.frombuffer(fmem, dtype=HEADER, count=nhdrs - 2, offset=fmem_offset)

                hmchunk = fmem[fmem_offset:fmem_offset + 3 * N_HEADER_BYTES]
                hmchunk_cp = np.zeros(N_HEADER_BYTES, dtype='u1')


                for i in range(N_HEADER_BYTES):
                    #print i, hmchunk[i:N_HEADER_BYTES], mode(hmchunk[i::N_HEADER_BYTES])[0]
                    hmchunk_cp[i] = mode(hmchunk[i::N_HEADER_BYTES])[0]

                header = np.frombuffer(hmchunk_cp, dtype=HEADER, count=1)
                # I am not 100% sure this is a good idea....it could still be the case the headers CRC does not check but the header is mostly valid and all these files would be unreadable

                header_is_good = check_header_crc(buff_str=header.data[:-2],crc_value=header['CRC'].item())
                if header_is_good:
                    return headers, header, 150-(2-block_offset)*N_HEADER_BYTES + tracker
                else:
                    #raise Exception('Failed to find a valid header at offset %s' % fmem_offset)
                    return None, None, None


        else:
            if block_offset == -1:
                logger.debug('*'*100)
                logger.debug(150-(2-block_offset)*N_HEADER_BYTES + tracker)
            return headers[hcrc], good_headers[hindex], 150-(2-block_offset)*N_HEADER_BYTES + tracker



def inspect_gvar(fmem=None):
    """
    Inspects a 1D u1 (byte) array representing an in memory GOES GVAR file
    and collects information about the channels present in the GVAR as well as their dimensions.
    This function could be easily modified to return detector configuratin of the instrument.
    should that information be required


    :param fmem: numpy 1D byte array, the mem mampped GVAR data
    :return: a dict with infor on channels dimension
    """
    nbits_in_byte = 8
    fmem_offset = 0
    ir_channels_det = dict()
    detectors = dict()
    chn_dims = dict()
    ir_channels_npix = []
    logger.debug('Inspecting...')
    break_bcount = 0
    nb0 = 0  # this is the number of time b0 was encountered. WE nend thi to return corect vislines because every subsequent b0 has smaller vislines (- 8 ) - - 1 scan
    while True:

        # logger.debug(fmem_offset)
        headers, header, block_offset = find_header(fmem=fmem, fmem_offset=fmem_offset)  # nhdrs=5, default
        # logger.debug(headers)
        if header is None:
            # one can assume in case a header was not found the end of fiel was reached
            logger.debug('Finished inspecting. Bad header')
            break

        # update fmemoffset
        fmem_offset += block_offset

        data_valid = header['DataValid']
        bid = 0 if header['BlockID'] == 240 else header['BlockID']
        bcount = header['BlockCount'].item()
        w_size = header['WordSize'].item()
        w_count = header['WordCount'].item()
        version = header['VersionNumber'].item()
        #logger.debug('%s %s' % (data_valid, bid))

        bytes_to_read = int((w_count - 2) * (w_size / float(nbits_in_byte)) + 2)
        if len(fmem) - fmem_offset < bytes_to_read:
            raise IOError('GVAR %s has a bad size' )
        block_data = fmem[fmem_offset:fmem_offset + bytes_to_read]

        if not data_valid:
            fmem_offset += bytes_to_read
            logger.debug('Skipping block %d count no %d. Invalid data' % (bid, bcount))
            if bcount < break_bcount:
                break_bcount += 1
                logger.debug('Setting break bcount because of invalid data to %d in %d' % (break_bcount, bcount))
            continue  # skip, makes no sense to go on

        good_block_data = crc_hqx(block_data.data, 0xffff) == MAGIC_CRC
        block_words = eight_bit2ten_bit(in_array=block_data)



        # update fmemoffset
        fmem_offset += bytes_to_read


        if bid == 0:
            nb0 += 1
            # we have to scan the file here, so no block90 replacement can occur
            # create block 0
            b0 = block0.from_buffer_copy(block_data.data)
            if good_block_data:
                break_bcount = bcount + 4
            else:
                if not b0_parity_check(block_data):  # section 1 parity check passed, lets be forginving
                    break_bcount == bcount + 14
                else:
                    break_bcount = bcount + 4
            logger.debug('Setting break bcount to %d in %d' % (break_bcount, bcount))

            # he have a good b0 block
            # let's extract vis lines and cols
            vis_nlines = b0.vis_lines
            vis_ncols = b0.vis_cols
            logger.debug('vislines %d viscols %d' % (vis_nlines, vis_ncols))

            # lines are multiple of 8 and cols are of 4, we assume

            h = vis_nlines - vis_nlines % 8 + 8
            w = vis_ncols - vis_ncols % 4 + 4



        elif bid in (1, 2):
            if not data_valid:
                logger.debug('block %d count %d has invalid data' % (bid, bcount))
                break_bcount += 1  # increment so next block is scanned
                continue
            # IR channels
            line_in_block = 0
            words_to_go = block_words.size
            chn_offset = 0

            while words_to_go > 2:  # as per gvar tool but why? really


                ldoc_16 = block_words[chn_offset:chn_offset + N_LDOC_BYTES]
                if ldoc_16.size != N_LDOC_BYTES:# perhaps some extra bytes
                    logger.debug('Skiping a line at bcount %d' % bcount)
                    break

                ldoc = wrap_ldoc(ldoc_16)
                n_words = ldoc['LWORDS']
                chn = ldoc['LICHAA']
                lidet = ldoc['LIDET']
                npix = ldoc['LPIXLS']
                #print chn, lidet, npix

                #TODO find a way to add validation cond for channel. from predefined constants hard wired...
                if n_words == 0 or npix == 0  or lidet not in range(1, 9):
                    logger.debug('Invalid line documentation encountered at block number %s of type %s' % (bcount, bid))
                    break
                if not chn in detectors:
                    detectors[chn] = [lidet]
                else:
                    if not lidet in detectors[chn]:
                        detectors[chn].append(lidet)
                if not chn in ir_channels_det:
                    ir_channels_det[chn] = [lidet]
                else:
                    ir_channels_det[chn].append(lidet)
                ir_channels_npix.append(npix)
                line_in_block += 1
                words_to_go -= n_words
                chn_offset += n_words

        elif bid > 2 and bid < 11:
            if not data_valid:
                break_bcount += 1  # increment so next block is scanned
                logger.debug('block %d count %d has invalid data' % (bid, bcount))
                break_bcount+=1
                continue
            ldoc_16 = block_words[:N_LDOC_BYTES]
            if ldoc_16.size != N_LDOC_BYTES:  # perhaps some extra bytes
                logger.debug('Skiping vis block at bcount %d' % bcount)
                break_bcount+=1
                continue

            vis_ldoc = wrap_ldoc(block_words[:N_LDOC_BYTES])
            chn = vis_ldoc['LICHAA']
            vis_npix = vis_ldoc['LPIXLS']
            lidet = int(vis_ldoc['LIDET'])

            if vis_npix == 0  or lidet not in range(1, 9):
                logger.debug('Invalid line documentation encountered at block number %s of type %s' % (bcount, bid))
                break_bcount+=1
                continue

            if not chn in detectors:
                detectors[chn] = [lidet]
            else:
                if not lidet in detectors[chn]:
                    detectors[chn].append(lidet)

            logger.debug('vis npix %d' % vis_npix)
        logger.debug('bid, %d bcount %d, good_data %s, break at %d' % (bid, bcount, good_block_data, break_bcount))
        #print 'inspect', len(ir_channels_npix)
        #if bcount == break_bcount and len(ir_channels_det) == 4:

        if bcount == break_bcount and len(detectors[chn]) >= 8:
            logger.debug('Finished inspecting at %d %d ' % (bid, bcount))
            break



    # now we have info to create a dict of channels and their resolutions
    # we got vis lines and cols from b0
    # match the ncols with number of pixels in first vis block (3, 10)
    # compute the n of detectors for each ir channel. AS GOEs has 5 channels one channel will have 1 detector and have
    # thus half of the line or vertical resolution of the other channels. Because we know the IR channels have 4 times smaller res than vis
    # we know the line resolution of ir channels and can calculate the nlines for each IR channel
    # finally we collect the npixels from ldoc for IR channels
    if nb0 == 1 and bid < 3:
        raise Exception('Corruped GVAR file at byte offset %s ' % fmem_offset )
    vislines = h if nb0 <= 1 else nb0 * 8 + h
    #vislines = vis_nlines
    viscols = w if w == vis_npix else vis_npix
    #viscols =vis_ncols
    #
    # n_ir_detectors = [len(e) for e in ir_channels_det.values()]
    # max_ndet = max(n_ir_detectors)
    # line_ir_ratios = [e / float(max_ndet) / 4 for e in n_ir_detectors]
    # chn_dims[1] = vislines, viscols, 1,1
    # for i, chn in enumerate(ir_channels_det):
    #     chn_dims[chn] = int(vislines * line_ir_ratios[i]), ir_channels_npix[i], int(round(vis_nlines/(vislines * line_ir_ratios[i]))), 4
    #
    # b0.channels_shape = chn_dims



    no_detectors = [len(e) for e in detectors.values()]

    ir_npix = list(set(ir_channels_npix))
    if len(ir_npix) == 1:
        ir_npix = ir_npix[0]
    else:
        npd = dict()
        for e in ir_npix:
            npd[e] = ir_npix.count(e)
        max_v = max(npd.keys())
        ir_npix = npd[max_v]


    ratios = [ e/8 for e in no_detectors]
    yres = [ 8/e for e in no_detectors]
    for i, chn in enumerate(detectors):
        chn_nlines = int(vislines * ratios[i])
        chn_yreskm = int(yres[i])
        chn_xreskm = 1 if chn == 1 else int(round(viscols / ir_npix))
        chn_ncols = int(round(viscols/chn_xreskm))

        chn_dims[chn] = chn_nlines,chn_ncols, chn_yreskm, chn_xreskm

    b0.channels_shape = chn_dims
    b0.detectors = detectors



    #compute theoretical scan duration
    #it looks like the ISTIM can  be used to compute the scan duration in seconds
    # in folloing way (vis_nlines // 8. * (b0.istim+200)) / 1000.

    # the 200 milisecs is the observed delay for every swath  and every image disregarding to scene type so it is added to the ISTIM
    #however for scenes different than FD there seems to be something else going on every roughly 340 visible lines
    #this event takes in general around four seconds. I have noticed some scenes have a value of 4.5 secs . I have no idea what this delay really is
    #it could be many things.

    unknown_delay = (vis_nlines // 340) * 4 # seconds

    if vis_ncols > 20000: # pretty safe to assume FD scene
        sc = (vis_nlines / 8. * (b0.istim+200)) / 1000.
    else:
        sc = (vis_nlines / 8. * (b0.istim + 200)) / 1000. + unknown_delay

    b0.theoretical_scan_duration = datetime.timedelta(0,sc)


    '''
    31 Jan 2018. GOES14 %15 fix.
    It looks like for some reason the  instrument data was moved from words 6304-6308
    to 8033-8036

    Also some fix for GOES14 specifically (GVAR version 3)

    3.  Imager Factory Coefficients in GVAR:   Beginning with GOES-O, the
        imagers will have an additional infrared detector. There will be eight,
        instead of seven, IR detectors.   In block 0, the amount of space for
        the drift correction coefficients (currently words 5587-6304) will need
        to be increased by approximately 15% to accommodate the data for the
        eighth detector.  To make room in block 0, the factory coefficients
        (currently words 6305-8040) will be removed and will be sent instead in
        a new type of block 11.  [The reason we retain the drift correction data
        but not the factory data in block 0 is that the drift correction data
        are updated with each scan line, whereas the factory coefficients remain
        the same for a given satellite for all time.  (And, of course, blocks 0
        are sent with each scan line, whereas blocks 11 are sent less often.)]

        ***The following was added in March 2001: The new "Imager Factory
        Coefficients" block 11 will actually be introduced with the GOES-M
        satellite.  This block 11 will be identified by the value "20" in words
        5 and 6 ("Product ID") of the GVAR block header.  However, the changes
        to block 0, including the removal of the Imager factory data and the
        addition of drift correction data for the eighth IR detector, will not
        become effective until GOES-O.***

        4. Versions of GVAR: The changes to GVAR outlined above will not be
        retroactive.  Data from GOES-8 through GOES-11 will always be
        transmitted to users with the current version of GVAR.  The revisions to
        GVAR described above will become effective with the GOES-M satellite.

        ***The following was added in March 2001: Users can identify the version
        of GVAR they receive from word 8 ("Version Number") of the GVAR block
        header.   When this word has the value "0," it indicates the original
        (GOES-8 through GOES-11) version of GVAR.  The value "1" will also
        indicate the GOES-8 through -11 version of GVAR, but with certain spare
        words now utilized to transmit the status of a new algorithm to mitigate
        certain calibration anomalies that occur near satellite midnight.  NOAA
        expects to make version 1 operational in mid- to late-2001.  Its receipt
        will not require any modifications or actions by current users.   For
        the GVAR versions for GOES-M/N and GOES-O/P, word 8 of the GVAR block
        header will take the values "2" and "3,"  respectively. ***


    '''

    #fix nadir location
    if b0.spcid >=14:
        b0.ns_cycles = b0.iofnc
        b0.ew_cycles = b0.iofec
        b0.ns_incr = b0.iofni
        b0.ew_incr = b0.iofei

    #sepcific fix for GOES 14
    if b0.spcid == 14 and (b0.iofnc == 0 or b0.iofec == 0 or b0.iofni==0 or b0.iofei == 0):
        b0.ns_cycles = 4
        b0.ew_cycles = 2
        if b0.iscan.imc_active:
            b0.ns_incr = 3013
            b0.ew_incr = 3192
        else:
            b0.ns_incr = 3123
            b0.ew_incr = 2944

    return b0

#@timeit()
def parse_gvar(gvar_file=None, channels_to_extract=None, fill_missing_lines=True, create_own_index=False,read_all=False):
    """
    A functional interface to a GOES GVAR imager file. Pure python (numpy) and decently fast

    GVAR is a real time stream based format and thus needs a specific approach, typical for handling streams.
    Conceptually,  information that comes through the stream at a given point in time provides clues into the future or what is about to come through the pipe.
    For example the header of a block contains information related to how much data/bytes the block contains.
    If this information is corrupted  or lost it is almost impossible to continue reading the stream. This is why parsing GVAR data is difficult.

    Simply put the stream consists of an array of blocks. Some of these blocks  carry only metadata (block 0),
    while other  blocks carry data (1-10 for Imager).

    A block consists of a header and payload. The block header comes in triple redundancy precisely to facilitate the parsing.

    GOES Imager works by scanning the space vertically (North south cycles) and horizontally (east west cycles)
    Each NS cycle consists of a series of scan swaths. A swath can be described as a one time snapshot(in reality the scan swath duration varies with the scene width) where one block0 and
    10 blocks carrying data are generated by the onboard cameras.
    Thus a swath contains 8 visible lines of data in 8 blocks (3-10) and one line of data for  four infrared channels,
    that is a Block 1 or 2 contains data for two channels or more than one  imfrared image line.


    This function approaches the parsing in the following way:

    1. first the file is mapped into memory through memmap to boost the speed of parsing
    2. the stream is inspected  at its very begining (until all medatatda that is required is collected) and metadata is extracted through inspect_gvar function.
    This medatada is used  further to parse the file
    3. Using a while loop the file is parsed block by block until there are no more bytes to parse. Blocks are extracted at every iteration
    There a three types of situations/blocks the code handles:
        a) parses a block 0. This basically gives out information if the stream is not corrupted al all other metadata (scan time , etc)
        b) blocks that carry visible channel data (3-10)
        c) block that carry IR channels (2,3,4,5,6) data
    The data is then extracted, converted from 8 bit words(python reads 8 bit words)  to 10 bit words (GVAR uses 10 bit words) and stored
    into numpy arrays.
    Compared to other software that can read GVAR streams pygvar uses a pythonic approach (it expects things to be as they should)
    and does a lot of guessing and trying).
    McIDAS X and gvar tool use more a  brut force approach and rely on speed of  C and Fortran
    to achieve the same.



    :param gvar_file: str,
    :param channels_to_extract:
    :param fill_missing_lines: bool
    :param create_own_index: bool, False
    :return:
    """
    if create_own_index:
        index_file = f'{gvar_file}.mindx'
    assert os.path.abspath(gvar_file), 'The "gvar_file"  %s  must be an absolute path' % gvar_file
    assert os.path.exists(gvar_file), 'The "gvar_file" %s a does not exist' % gvar_file
    nbits_in_byte = 8
    logger.debug('Parsing %s' % gvar_file)
    try:
        f = np.memmap(gvar_file,mode='r')
        if create_own_index:
            mif = open(index_file, 'w')
    except Exception as e:
        fsize  = os.path.getsize(gvar_file)
        raise Exception('Looks like something went wrong while mapping %s to memory %s. The size of the GVAR is %d' % (gvar_file, e, fsize))
    fmem_offset = 0
    try:
        #fetch channels info as well as block0
        fb0 = inspect_gvar(fmem=f)
    except Exception as inspe:
        logger.error('GVAR %s is invalid. %s' % (gvar_file, inspe))
        raise
    b0 = fb0
    #return b0 of user did not ask for data
    channels_dims = fb0.channels_shape
    channels_data = {}
    detector_lines = {}
    missing_lines = {}
    prev_lines = {}
    for chn_n, dims in channels_dims.items():
        logger.debug(' channel %s dimensions %s ' % (chn_n, str(dims)))

    if not channels_to_extract:
        if not read_all:
            fb0.scan_duration = fb0.theoretical_scan_duration
            return fb0
        else:
            channels_to_extract = list(channels_dims.keys())
    for channel_to_extract in channels_to_extract:
        assert channel_to_extract in channels_dims, 'Invalid channel %s. Valid channels are %s' % (channel_to_extract, channels_dims.keys())
        channel_lines, channel_cols, chn_lineres, chn_colres = channels_dims[channel_to_extract]
        channels_data[channel_to_extract] = np.zeros((channel_lines, channel_cols), dtype='int16')
        missing_lines[channel_to_extract] = list()
        detector_lines[channel_to_extract] = np.empty((channel_lines,), dtype='u1')


    #main loop
    nb = 0
    while True:

        perc = 100. * fmem_offset / f.size

        if int(perc) == 100:
            break
        #logger.debug(fmem_offset)
        headers, header, block_offset  = find_header(fmem=f, fmem_offset=fmem_offset) #nhdrs=5, default
        if header is None:
            if perc > MIN_USABLE_IMAGE_PERCENTAGE:
                logger.debug('Failed to find a valid header at %.2f %% in %s' % (perc, gvar_file))
                break
            logger.error('Failed to find a valid header at %.2f %% in %s' % (perc, gvar_file))
            raise(Exception('Failed to find a valid header at %.2f %% in %s' % (perc, gvar_file)))
            break
        #print(header2dict(header))
        #logger.debug(headers)
        if header is None:
            # one can assume in case a header was not found the end of file was reached
            logger.debug('Finished scanning')
            break

        data_valid = header['DataValid'].item()
        bid = 0 if header['BlockID'] == 240 else header['BlockID'].item()
        bcount = header['BlockCount'].item()
        w_size = header['WordSize'].item()
        w_count = header['WordCount'].item()


        nb+=1

        if not bcount%1000:
            logger.debug('read %.2f %%' % perc)


        # u[date fmemoffset
        fmem_offset += block_offset

        bytes_to_read = int((w_count - 2) * (w_size / float(nbits_in_byte)) + 2)



        if len(f) - fmem_offset < bytes_to_read:
            raise IOError('GVAR %s has a bad size' % gvar_file)
        block_data = f[fmem_offset:fmem_offset+bytes_to_read]

        #logger.debug('bid %d  bcount  %d boffset %d' % (bid, bcount, block_offset))

        #convert 8 bit data to 10 bit
        block_words = eight_bit2ten_bit(in_array=block_data)

        #update fmemoffset
        fmem_offset+=bytes_to_read

        if bid == 0:

            good_block0_data = crc_hqx(block_data.data, 0xffff) == MAGIC_CRC
            #print bid, bcount, bool(data_valid), good_block0_data
            # create block 0
            if good_block0_data:
                last_good_b0 = b0

            pb0 = b0
            b0 = block0.from_buffer_copy(block_data.data)
            if create_own_index:
                mif.write(f'{os.path.split(gvar_file)[1]}  {(bytes_to_read + block_offset) * nbits_in_byte} {bid} {b0.insln} \n')

            '''
            KISS. As long as the code is simple things ought to work
            It's not a good idea to mess around with block 0. That is why the next code is commented out
            '''

            # parc =  b0_parity_check(block_data)
            # if not parc:
            #     logger.info('bid %s  good %s, parc %s' % (bid, good_block0_data, parc))
            # if good_block0_data:
            #
            #
            # else:
            #
            #     # if something goes wrong in this scan, in case parity check failes then
            #     # here is the place to set blindly an invalid b0 block
            #     '''
            #         1-278 Instrument and Scan Status
            #         279-1626 Instrument Orbit and Attitude (O&A) Data
            #         1627-2306 Scan Reference Data
            #         2307-5386 Grid Data
            #         5387-6304 Scan Reference and Calibration Data
            #         6305-8040 Factory Parameters
            #
            #
            #         '''
            #     #looks like this parity check is  the way to goo.
            #     if b0_parity_check(block_data):
            #         # section 1 parity check passed, lets be forgiving. The section 2 wich contains nav can not be in my opinion
            #         # identical because it HAS different scan time. But is is still better to continue with original block0
            #         continue
            #     else:
            #         '''
            #         This is a tricky situation.....
            #         crc and parity check failed, normally the main loop breaks.
            #
            #         However, this resulted in too many gvar files being read only partially. So instead I choose to
            #         assign the invalid b0 to last valid b0 and increment by 8 it's northermost visible line as relative line /scan count by 1
            #         For many gvar this did the trick. But I have encountered gvars where the result of this jump was the main loop never returned.
            #         for example goes13.2017.080.084518 is one such gvar file
            #
            #
            #
            #         '''
            #
            #         #break
            #         b0 = pb0
            #                 #print 'bid %d bcount %f relvline %d  lidet %d vis_line %d bl %d' % (bid, bcount, rel_vis_line, lidet, vis_line, bl)
                #logger.debug('bid %d bcount %f relvline %d vis_line %d bl %d' % (bid, bcount, rel_vis_line, vis_line, bl))
            #         b0.insln += 8 #adjust scan northemost line, this is important because otherwise the line index si not computed correctly
            #         b0.risct +=1
            #         b0.aisct +=1
            #     # and the result is a black line


        #logger.info('bl %d bid %s' % (bl, bid))
        # lon0 = (fb0.reflo.python_value() / RAD)
        # print(lon0)
        if bid in (1, 2): # IR data
            if not data_valid:
                logger.debug('block %d count %d has invalid data' % (bid, bcount))
                continue
            if create_own_index:
                mif.write(
                    f'{os.path.split(gvar_file)[1]}  {(bytes_to_read + block_offset) * nbits_in_byte} {bid}  \n')
            #IR channels clocks contains data from mutiple channels
            line_in_block = 0
            words_to_go = block_words.size
            chn_offset = 0
            while words_to_go-2 > 0: # the main condition to conetinue parsing
                #raw line doc
                ldoc_16 = block_words[chn_offset:chn_offset + N_LDOC_BYTES]
                #check its lenght
                if ldoc_16.size != N_LDOC_BYTES:# perhaps some extra bytes
                    logger.debug('Skiping a line at bcount %d' % bcount)
                    break
                #looks like it is
                ldoc = wrap_ldoc(ldoc_16)
                n_words = int(ldoc['LWORDS'])
                chn = int(ldoc['LICHAA'])
                lidet = int(ldoc['LIDET'])
                npix = int(ldoc['LPIXLS'])
                #validate line, pixels, chn
                if n_words == 0 or npix == 0 or chn not in fb0.channels_shape or lidet not in range(1,9):
                    logger.debug('Invalid line documentation encountered at block number %s of type %s' % (bcount, bid))
                    break

                if ldoc['RISCT'] != b0.risct:  # necessary
                    logger.debug('risct out of sync %s %s for channel %s' % (ldoc['RISCT'], b0.risct, chn))
                    break

                #logger.info('bcount %s chn %s  nwords %s, words_to_go %s' % (bcount,chn,n_words, words_to_go))
                #looks like we have a good line doc structure. This is not always the case and that is why we need some care inside next loop but if we skip this out too many lines will be  not read
                if chn in channels_to_extract:
                    chn_res = channels_dims[chn][2]
                    chn_nl, chn_nc = channels_dims[chn][:2]
                    #code from gvar tool
                    det = (lidet+3)&7 if chn == 1 else (lidet-1)&1
                    if b0.iscan.yaw_flip: #my idea about flipping
                        det = ~det+2 #negate and append 2 because there are 2 channels per ir line . It seems to be working. TODO handle ldoc lagged scan
                    #let's compute the visible lie number
                    rel_vis_line = chn_res * det + chn_res/2.-.5
                    vis_line = (int(b0.insln + rel_vis_line)) - b0.infln
                    #get the infrared line number by dividing the vis using resolution
                    ir_line = vis_line//chn_res


                    # this is the pythonic alternative to brute force. By keeping the track of current and previous line number
                    #missing lines can be detected and interpolated which greatly improves the quality of images

                    try:
                        pir_line = prev_lines[chn]
                        if ir_line < pir_line:
                            tdl = ir_line - pir_line
                            if tdl == 1:
                                ir_line = pir_line+1
                            else:
                                break
                        dl = ir_line-pir_line
                        if dl > 1 and dl <= MAX_MISSING_IRLINES_TO_INTERPOLATE:
                            tml = ir_line-dl, ir_line-1
                            missing_lines[chn].append(tml)
                    except KeyError as ke:
                        pass

                    prev_lines[chn] = ir_line
                    if ir_line >= chn_nl: #something went bad
                        logger.debug('line is out of sync')
                        break # ??? this should be safe but will it read all data?? I guess yes


                    start = chn_offset + N_LDOC_BYTES
                    end = start + npix
                    ir_line_data = block_words[start:end]
                    #if ir_line >= 20 and ir_line <= 29:
                    #print({'channel':chn, 'line_no':ir_line, 'detector':lidet, 'det_no':det}, b0.iscan.side2)


                    #make sure again things are under control
                    if npix >= chn_nc: #there are more columsn than header says so we subset the pixels with difference
                        npix +=chn_nc-npix
                    if ir_line_data.size != npix:
                        pass

                    #finally read the data
                    channels_data[chn][ir_line, :npix] = ir_line_data[:npix]
                    detector_lines[chn][ir_line] = det
                    #logger.debug('ir_line %s chn_nl %s, vis_line %s ' % (ir_line, chn_nl, vis_line))
                line_in_block += 1
                words_to_go -= n_words
                chn_offset += n_words

        if bid > 2 and bid < 11: #vis channel data blocks
            if not data_valid:
                logger.debug('block %d count %d has invalid data' % (bid, bcount))
                continue
            if create_own_index:
                mif.write(
                    f'{os.path.split(gvar_file)[1]}  {(bytes_to_read + block_offset) * nbits_in_byte} {bid}  \n')
            #make sure we have a valid line doc structure
            ldoc = wrap_ldoc(block_words[:N_LDOC_BYTES])
            if ldoc_16.size != N_LDOC_BYTES:  # perhaps some extra bytes
                logger.debug('Skiping a line at bcount %d' % bcount)
                continue

            npix = int(ldoc['LPIXLS'])
            chn = int(ldoc['LICHAA'])
            lidet = int(ldoc['LIDET'])

            if npix == 0 or chn not in fb0.channels_shape or lidet not in range(1,9):
                logger.debug('Invalid line documentation encountered at block number %s of type %s' % (bcount, bid))
                continue
            #everything looks good so let's read
            if chn in channels_to_extract:
                chn_res = channels_dims[chn][2]
                chn_nl, chn_nc = channels_dims[chn][:2]

                det = (lidet+3)&7 if chn == 1 else (lidet-1)&1
                #in yaw flip mode for vis channel it suffices to reverse the detector TODO handle ldoc lagged scan
                if b0.iscan.yaw_flip:
                    det = 7-det
                #get the visible line number
                rel_vis_line = chn_res * det + chn_res/2.-.5
                vis_line = (int(b0.insln + rel_vis_line)) - b0.infln

                # this is the pythonic alternative to brute force. By keeping the track of current and previous line number
                # missing lines can be detected and interpolated which greatly improves the quality of images
                try:

                    if vis_line < pvis_line:
                        tvdl = vis_line - pvis_line
                        if tvdl == 1:
                            vis_line = pvis_line + 1
                        else:
                            break
                    vdl = vis_line - pvis_line
                    tml = vis_line - vdl + 1, vis_line-1
                    if vdl > 1 and vdl <= MAX_MISSING_VISLINES_TO_INTERPOLATE:
                        missing_lines[chn] .append(tml)
                except UnboundLocalError:
                    pass
                pvis_line = vis_line
                #it's bettre to check than be surprised
                if vis_line >= chn_nl:
                    break  # ??? this should be safe but will it read all data i guess yes
                line_data = block_words[N_LDOC_BYTES:N_LDOC_BYTES+npix]
                if pb0.iscan.frame_end:
                    logger.debug('frame end %s ' % bcount)
                    break
                # print ir_line_data.shape, npix
                if npix >= chn_nc:  # there are more columns than the header says so we subset the pixels using their delta
                    npix += chn_nc - npix
                if line_data.size != npix:# skip thi line, something went bad
                    pass # just do nothing

                #print({'channel':chn, 'line_no':vis_line, 'detector':lidet, 'det_no':det}, b0.iscan.side2)

                #finally read data
                channels_data[chn][vis_line, :npix] = line_data[:npix]
                detector_lines[chn][vis_line] = det






    #scan duration

    try:
        if abs(vis_line - fb0.channels_shape[1][0]) > 40:  # we accept +- 5 swaths = 5*8 lines per swath
            fb0.scan_duration = fb0.theoretical_scan_duration

        else:
            end_scan = b0.tched.python_value()
            fb0.scan_duration = end_scan - fb0.tched.python_value()
    except Exception as e1:
        fb0.scan_duration = fb0.theoretical_scan_duration


    # if missing lines were encountered  interpolate them
    if fill_missing_lines: # is user requested and missing lines exist interpolate linearly the missing lines
        for ccc, miss_lines in missing_lines.items():
            if miss_lines:
                for min_ml, max_ml in miss_lines:
                    start_l = min_ml - 1
                    end_l = max_ml + 1
                    start_l_data = channels_data[ccc][start_l, :]
                    end_l_data = channels_data[ccc][end_l, :]
                    for mml in range(min_ml, max_ml+1):
                        channels_data[ccc][mml,:] = start_l_data+ np.round((end_l_data-start_l_data)/(float(end_l)-mml)).astype(np.int32)
    del perc, fmem_offset, b0, headers, header, missing_lines # leaks because python 2 leaves the vars in the namespoace after any loop has finished
    if create_own_index:
        mif.close()
    fb0.detector_lines = detector_lines
    return fb0, channels_data

def parse_gvars_multichannel(folder=None, channels_to_extract=[], view=False, save_to_tiff_folder=None):
    """
    Util function that allow visualisation of GOES   channels from a GOES GVAR file (only Imager) locate din a folder

    :param folder:  str, folder where gvar resides
    :param channels_to_extract: iterable of channels
    :param view: bool, if true the channels will be showed using matplotlib
    :param save_to_tiff_folder:, str path to the fodler where each channel will be writtine to tiff format . requires gdal python  bindings
    :return:
    """
    if view:
        from pylab import imshow, show, title
    for gvar_name in sorted(os.listdir(folder)):
        # if not 'goes13' in gvar_name:
        #     continue
        sample_gvar = os.path.join(folder, gvar_name)

        b0, channels_data = parse_gvar(sample_gvar, channels_numbers=channels_to_extract)

        if channels_data is None:
            continue


        if view:
            for channel_to_extract, data in channels_data.items():
                title('Channel %d from %s ' % (channel_to_extract, gvar_name))
                imshow(data[::1, ::1], interpolation='nearest')
                show()

        if save_to_tiff_folder:
            assert os.path.abspath(save_to_tiff_folder), 'save_to_tiff_folder has to be an absolute path to a folder'
            assert os.path.exists(save_to_tiff_folder), 'save_to_tiff_folder does not exist'
            for channel_to_extract, data in channels_data.items():
                tif_path = os.path.join(save_to_tiff_folder, gvar_name+'B%s.tif' % channel_to_extract)
                logger.debug('Writing ' % tif_path)
                toGDAL(tif_path, array=data )
