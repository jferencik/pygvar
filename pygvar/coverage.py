COVERAGES_EAST = dict(
                        FD=(2472, 13288, 5852, 24829),
                        US=(2928, 6368, 9052, 17752),
                        USE=(2696, 7342, 9052, 19052),
                        NH=(3096, 8345, 9052, 22500),
                        NHE=(2696, 10000, 9052, 22900),
                        SH=(9928, 12192, 9052, 22900),
                        SAC=(8345,10104,14050,20700),
                        SAS = (10104,12184,15200,20700)
                    )

COVERAGES_WEST = dict(
                        NH=(2485, 7905, 7829, 21077),
                        US=(2605, 6537, 10037, 21461),
                        SH=(7837, 12089, 7829, 18885),
                        SUS=(3069, 5457, 15941, 20917),
                        FD=(2485, 13305, 5789, 26637)
                    )
REGIONS = ('east', 'west')

def get_region(ssp=None):
    if ssp <=-105: #105 is (135-75)/2
        return 'west'
    else:
        return 'east'


# ########################IMPORTS ########################################
def get_theoretical_coverage(h=None, m=None, region=None):
    """
    GOES <=15 theoretical coverage
    The imaging of a specific group of standard sectors executed in a particular sequence constitutes
    a GOES imaging operational scenario. The NOAA NESDIS/NWS Study Group has defined a
    set of sectors of the imager full-earth field of view and three imager operational scenarios that
    satisfy NWS requirements for the collection of Imager data. The operational modes are
    designated Routine, Rapid Scan and Full Disk. Tables 3-10 and 3-11 depict GOES in Routine
    operational mode.
    The three imaging modes correspond, respectively, to operation of the GOES system under
    normal or typical meteorological conditions, operation under conditions of one or more evolving
    severe weather or tropical storm conditions, and operation during periods of degraded GOES
    system performance, such as during system component failure or maintenance, or during the
    spring and fall eclipse season when GOES scanning is halted to avert solar intrusion on the
    radiometer.
    In the GOES routine schedule scan mode, two views at approximately 15 minute intervals of the
    CONUS (GOES-East) or PACUS (GOES-West) are provided in a half hour period. A northern
    hemisphere scan is also included in the 30 minute cycle.

    During GOES Rapid Scan Operations (RSO), four views of the CONUS (GOES-8) or
    Sub-CONUS (GOES-9) are provided at approximately 7.5 minute intervals in a half hour period.
     A northern hemisphere scan for both GOES is also included in the 30 minute cycle. This yields
    eight views of the continental U.S. per hour.

    During GOES Super Rapid Scan Operations (SRSO), approximately 10 one-minute interval
    scans are provided every half hour using prescribed 1000 km x 1000 km sectors. The remaining
    time in the half hour cycle is devoted to scans of the northern hemisphere and CONUS
    (GOES-East) or Sub-CONUS (GOES-West).

    When GOES RSO or SRSO is utilized, most of the southern hemisphere is not
    scanned.

    --------------------------------------------
    | Full Earth | Earth Edge 26:10 0000, 0300, etc
    ----------------------------------------------
    Northern Hemisphere 0-66N/90W-170E 9:00 xx00, xx30
    Southern Hemisphere 0-45S/115W-170E 7:00 xx22, xx52
    PACUS 12-60N/90-175W 5:00 xx15, xx45

    :param h: hour
    :param m: minute
    :return: str, the coverage of the frame
    """
    assert region.lower() in REGIONS, f'Invalid region {region}. Valid values are {"".join(REGIONS)}'
    reg = region.lower()
    h, m = map(int, (h, m))
    if reg == 'west': #GOES WEST
        if m == 0:
            if h%3 == 0:
                return 'FD'
            else:
                return 'NH'
        elif m in list(range(10-2, 10+3)) + list(range(40-2,40+3)):
            return 'SUS'
        elif m == 30:
            return 'NH'
        elif m in (15, 45):
            return 'US'
        elif m in list(range(22-3, 22+4)) + list(range(52-3, 52+4)):
            return 'SH'

    else: #GOES EAST
        if m == 45:
            if h%3-2 == 0:
                return 'FD'
            else:
                return 'NH'
        elif m in range(15-2,15+3):
            return 'NH'
        elif m in [0,1,2]+ list(range(30-2,30+3)):
            return 'US'
        elif m in list(range(7-3, 7+3)) + (range(37-3, 37+3)):
            return 'SH'


def get_empirical_coverage(sil=None, eil=None, sie=None, eie = None, region=None):
    """
    Computes the empirical coverage of a GOES<=15 GVAR file using the bounding box
    of the VIS channel in image coordinates(line, ele)
    :param sil:
    :param eil:
    :param sie:
    :param eie:
    :param region: str, region, east or west
    :return:
    """

    coverages = COVERAGES_EAST if region.lower() == 'east' else COVERAGES_WEST
    d = {}
    for cname, cextent in coverages.items():
        covsil, coveil, covsie, coveie = cextent
        covw  = float(coveie-covsie)
        covh  = float(coveil-covsil)
        covhh = covsil + covh//2
        covp = 2*covw+2*covh
        w  = eie-sie
        h  = eil-sil
        hh = sil+h//2
        p = 2*w+2*h
        r1 = abs(1-covp/p)
        r2 = abs(hh-covhh)
        d[cname] = int(r1*r2)

    if len(d)>0:
        min_d = min(d.values())
        return [k for k, v in d.items() if v == min_d].pop()

