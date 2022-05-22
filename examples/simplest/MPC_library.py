# -*- coding: utf-8 -*-
import math
import numpy as np
import scipy

class Constants:
    GMsun = 2.9591220828559115e-04 
    Rearth_km = 6378.1363
    au_km = 149597870.700 # This is now a definition 
    Rearth_AU = Rearth_km/au_km
    ecl = (84381.4118*(1./3600)*np.pi/180.) # Obliquity of ecliptic at J2000
    #ecl = (84381.448*(1./3600)*np.pi/180.) # Obliquity of ecliptic at J2000
    speed_of_light = 2.99792458e5 * 86400./au_km

def rotate_matrix(ecl):
    ce = np.cos(ecl)
    se = np.sin(-ecl)
    rotmat = np.array([[1.0, 0.0, 0.0],
                  [0.0,  ce,  se],
                  [0.0, -se,  ce]])
    return rotmat

from novas import compat as novas
from novas.compat import eph_manager
from novas.compat import solsys
# This opens DE405 by default.
jd_start, jd_end, number = eph_manager.ephem_open()

class Observatory:

    # Parses a line from the MPC's ObsCode.txt file
    def parseObsCode(self, line):
        code, longitude, rhocos, rhosin, ObsName = line[0:3], line[4:13], line[13:21], line[21:30], line[30:].rstrip('\n')
        if longitude.isspace():
            longitude = None
        if rhocos.isspace():
            rhocos = None
        if rhosin.isspace():
            rhosin = None
        return code, longitude, rhocos, rhosin, ObsName

    def __init__(self):

        self.observatoryPositionCache = {} # previously calculated positions to speed up the process

        # Convert ObsCodes.txt lines to geocentric x,y,z positions and
        # store them in a dictionary.  The keys are the observatory
        # code strings, and the values are (x,y,z) tuples.
        # Spacecraft and other moving observatories have (None,None,None)
        # as position.
        ObservatoryXYZ = {}
        with open('ObsCodes.txt', 'r') as f:
            next(f)
            for line in f:
                code, longitude, rhocos, rhosin, Obsname = self.parseObsCode(line)
                if longitude and rhocos and rhosin:
                    rhocos, rhosin, longitude = float(rhocos), float(rhosin), float(longitude)
                    longitude *= np.pi/180.
                    x = rhocos*np.cos(longitude)
                    y = rhocos*np.sin(longitude)
                    z = rhosin
                    ObservatoryXYZ[code]=(x,y,z)
                else:
                    ObservatoryXYZ[code]=(None,None,None)
        self.ObservatoryXYZ = ObservatoryXYZ

    # The routine below calculates the geocentric position of the observatory
    # in equatorial cartesian coordinates.
    def getObservatoryGeoPosition(self, obsCode, jd_utc, xyz=None):

        # If obsCode==None, set the observatory to be the geocenter
        if not obsCode:
            obsCode = '500'

        if (obsCode, jd_utc) in self.observatoryPositionCache:
            return self.observatoryPositionCache[(obsCode, jd_utc)]
    
        obsVec = self.ObservatoryXYZ[obsCode]

        if obsVec[0]==None and xyz==None:
            #print("problem with obscode: ", obsCode)
            return None, None, None
        else:

            jd_tdb  = EOP.jdTDB(jd_utc)
            pos = getEarthPosition(jd_tdb)

            if obsCode=='500':
                geocentric_vec = np.zeros(3)
            elif xyz != None:
                geocentric_vec = xyz
            else:
                delta_t = EOP.delta_t(jd_utc)
                jd_ut1  = EOP.jdUT1(jd_utc)
                xp, yp  = EOP.pmx(jd_utc), EOP.pmy(jd_utc)
                geocentric_vec = Constants.Rearth_AU*np.array(novas.ter2cel(jd_ut1, 0.0, delta_t, xp, yp, obsVec))

            return geocentric_vec

    # The routine below calculates the heliocentric position of the observatory
    # in equatorial cartesian coordinates.
    def getObservatoryPosition(self, obsCode, jd_utc, xyz=None):

        # If obsCode==None, set the observatory to be the geocenter
        if not obsCode:
            obsCode = '500'

        if (obsCode, jd_utc) in self.observatoryPositionCache:
            return self.observatoryPositionCache[(obsCode, jd_utc)]
    
        obsVec = self.ObservatoryXYZ[obsCode]

        if obsVec[0]==None and xyz==None:
            #print("problem with obscode: ", obsCode)
            return None, None, None
        else:

            jd_tdb  = EOP.jdTDB(jd_utc)
            pos = getEarthPosition(jd_tdb)

            if obsCode=='500':
                geocentric_vec = np.zeros(3)
            elif xyz != None:
                geocentric_vec = xyz
            else:
                delta_t = EOP.delta_t(jd_utc)
                jd_ut1  = EOP.jdUT1(jd_utc)
                xp, yp  = EOP.pmx(jd_utc), EOP.pmy(jd_utc)
                geocentric_vec = Constants.Rearth_AU*np.array(novas.ter2cel(jd_ut1, 0.0, delta_t, xp, yp, obsVec))

            heliocentric_vec = pos+geocentric_vec
            self.observatoryPositionCache[(obsCode, jd_utc)] = heliocentric_vec
            return heliocentric_vec


def getSatellitePosition(xyz, jd_utc):
    x,y,z = parseXYZ(xyz)
    x /= Constants.au_km
    y /= Constants.au_km
    z /= Constants.au_km
    jd_tdb  = EOP.jdTDB(jd_utc)
    pos = getEarthPosition(jd_tdb)
    heliocentric_vec = pos[0]+x, pos[1]+y, pos[2]+z
    return heliocentric_vec        
    
class EarthAndTime:
    # Dealing with leap seconds and polar motion
    # ### Relating this to MPC data
    # 
    # I believe that the MPC observations have dates in UTC, which is the conventional thing to do.   
    # According to Gareth Williams, prior to 1972 Jan 1 the times are probably UT1.
    # 
    # So we take a time from an MPC observation.  If it's prior to 1972 Jan 1 we assume it's UT1 and 
    # we get delta_t (TT-UT1) from the historical table.  If it's on or after 1972 Jan 1 we determine 
    # the number of leap seconds from function below and then calculate delta_t.  
    # 
    def __init__(self, filename1='finals2000A.all', filename2='tai-utc.dat'):
        _xydeltat = {}
        with open(filename1) as f:
            _mjd    = slice(7, 15)
            _UT1_UTC= slice(58, 68)
            _pmx    = slice(18, 27)
            _pmy    = slice(37, 46)
    
            for line in f:
                if not line[_UT1_UTC].strip() == '':
                    _jd     = float(line[_mjd]) + 2400000.5
                    UT1_UTC= float(line[_UT1_UTC])
                    pmx    = float(line[_pmx])
                    pmy    = float(line[_pmy])
                    _xydeltat[_jd] = (UT1_UTC, pmx, pmy)

        jds = sorted(_xydeltat.keys())
        ut1_utcs = [_xydeltat[jd][0] for jd in jds]
        pmxs     = [_xydeltat[jd][1] for jd in jds]
        pmys     = [_xydeltat[jd][2] for jd in jds]

        nitems = len(jds)
        from scipy import interpolate        
        self.ut1_utc_func = interpolate.interp1d(jds, ut1_utcs, fill_value=(ut1_utcs[0], ut1_utcs[nitems-1]), bounds_error=False )
        self.pmx_func     = interpolate.interp1d(jds, pmxs, fill_value=(pmxs[0], pmxs[nitems-1] ), bounds_error=False )
        self.pmy_func     = interpolate.interp1d(jds, pmys, fill_value=(pmys[0], pmys[nitems-1] ), bounds_error=False )

        # Get the TAI-UTC data from:
        # http://maia.usno.navy.mil/ser7/tai-utc.dat
        self.tai_minus_utc_dict = {}
        with open(filename2) as f:
            _jd      = slice(16, 28)
            _tai_minus_utc = slice(36, 49)
            _tref    = slice(59, 66)
            _coeff   = slice(69, 79)

            for line in f:
                tai_minus_utc = float(line[_tai_minus_utc])
                jd      = float(line[_jd])
                tref    = float(line[_tref])
                coeff   = float(line[_coeff])
                self.tai_minus_utc_dict[jd] = tai_minus_utc, tref, coeff

    def pmx(self, jd_utc):
        #if jd_utc<244168:
        #    return self.pmx_func(2441684.5)
        return self.pmx_func(jd_utc)

    def pmy(self, jd_utc):
        #if jd_utc<2441684.5:
        #    return self.pmy_func(2441684.5)
        return self.pmy_func(jd_utc)

    # TDT = UTC + (TAI-UTC) + 32.184 sec
    # UT1 = UTC + delta_t
    # TT - TDB is less than 2 milliseconds.

    def ut1_utc(self, jd_utc):
        #if jd_utc<2441684.5:
        #    return self.ut1_utc_func(2441684.5)
        return self.ut1_utc_func(jd_utc)

    def tai_utc(self, jd_utc):
        if jd_utc<2437300.5:
            return 0.0
        ks = sorted(self.tai_minus_utc_dict.keys())
        m = max(i for i in ks if (i-jd_utc)<0.0)
        base, tref, coeff = self.tai_minus_utc_dict[m]
        tai_m_utc = base + coeff*(jd_utc-240000.5-tref)
        return tai_m_utc
    
    def jdTT(self, jd_utc):
        leaps = self.tai_utc(jd_utc)
        jd_tt = jd_utc + (32.184 + leaps)/(24.0*60*60)
        return jd_tt

    def jdUT1(self, jd_utc):
        DUT1  = self.ut1_utc(jd_utc)
        jd_ut1 = jd_utc + DUT1/(24.0*60*60)
        return jd_ut1

    def jdTDB(self, jd_utc):
        jd_tt = self.jdTT(jd_utc)
        _, tdb_tt = novas.tdb2tt(jd_tt)
        jd_tdb = jd_tt + tdb_tt/(24.0*60*60)
        return jd_tdb

    def delta_t(self, jd_utc):
        leaps = self.tai_utc(jd_utc)
        DUT1  = self.ut1_utc(jd_utc)
        delta_t = 32.184 + leaps - DUT1
        return delta_t
        
def getEarthPosition(jd_tdb):
    #pos, _ = solsys.solarsystem(jd_tdb, 3, 1)
    pos, _ = solsys.solarsystem(jd_tdb, 3, 0)    
    return pos

# This routine parses the section of the second line that encodes the geocentric satellite position.
def parseXYZ(xyz):
    xs = xyz[0]
    x = float(xyz[1:11])
    if xs=='-':
        x = -x
    ys = xyz[12]
    y = float(xyz[13:23])
    if ys=='-':
        y = -y
    zs = xyz[24]
    z = float(xyz[25:])
    if zs=='-':
        z = -z
    return x, y, z    


# Parses the date string from the 80-character record
def parseDate(dateObs):
    yr = dateObs[0:4]
    mn = dateObs[5:7]
    dy = dateObs[8:]
    return yr, mn, dy

# Converts the date string to a JD floating point number, using the NOVAS routine
# The time scale of the returned value will be same as that of the input date.
def date2JD(dateObs):
    yr, mn, dy = parseDate(dateObs)
    if ' ' in dy:
        dy=dy.split(' ')[0]
    hr, dy = math.modf(float(dy))
    jd = novas.julian_date(int(yr), int(mn), int(dy), 24.*hr)
    return jd

# These routines convert the RA and Dec strings to floats.
def RA2degRA(RA):
    if RA[0:2].strip() != '':
        hr = float(RA[0:2])
    else:
        hr = 0.0
    if RA[3:5].strip() != '':        
        mn = float(RA[3:5])
    else:
        mn = 0.0
    if RA[6:].strip() != '':        
        sc = float(RA[6:])
    else:
        sc = 0.0
    degRA = 15.0*(hr + 1./60. * (mn + 1./60. * sc))
    return degRA

def Dec2degDec(Dec):
    s = Dec[0]
    if Dec[1:3].strip() != '':
        dg = float(Dec[1:3])
    else:
        dg = 0.0
    if Dec[4:6].strip() != '':
        mn = float(Dec[4:6])
    else:
        mn = 0.0
    if Dec[7:].strip() != '':
        sc = float(Dec[7:])
    else:
        sc = 0.0
    degDec = dg + 1./60. * (mn + 1./60. * sc)
    if s == '-':
        degDec = -degDec
    return degDec

def convertEpoch(Epoch):
    yr0 = Epoch[0]
    yr1 = Epoch[1:3]
    mn  = Epoch[3]
    dy  = Epoch[4]
    if yr0=='I':
        yr0 = 1800
    elif yr0=='J':
        yr0 = 1900
    elif yr0=='K':
        yr0 = 2000
    else:
        print("Year error in convertEpoch")
        yr0 = 2000
    yr = yr0+int(yr1)
    if mn.isdigit():
        mn = int(mn)
    elif mn=='A':
        mn = 10
    elif mn=='B':
        mn = 11
    elif mn=='C':
        mn = 12
    else:
        print("Month error in convertEpoch")
        mn = 0
    if not dy.isdigit():
        dy = 10 + ord(dy) - ord('A')
    return yr, mn, int(dy)
    
def yrmndy2JD(yrmndy):
    yr, mn, dy = yrmndy
    hr, dy = math.modf(float(dy))
    jd = novas.julian_date(int(yr), int(mn), int(dy), 24.*hr)
    return jd

def deg2dms(v):
    minus_flag = (v<0.0)
    v = abs(v)
    v_deg = int(math.floor(v))
    v_min = int(math.floor((v-v_deg)*60.0))
    v_sec = (v-v_deg-v_min/60.)*3600
    if minus_flag:
        v_sgn = '-'
    else:
        v_sgn = "+"
    return v_sgn, v_deg, v_min, v_sec

def convert2MPC1992(trackID, line, mpNum="     ", disc=False):
    comps = line.split()
    detID, mJD, raDeg, decDeg, mag, filt, obsCode, vac, rej = \
        comps[0], comps[1], comps[3], comps[4], comps[5], comps[6], comps[7], comps[14], comps[15]
    JD = float(mJD)+2400000.5
    mag = float(mag)
    yr, mn, dy, hr = novas.cal_date(JD)
    raDeg = float(raDeg)
    decDeg = float(decDeg)
    day = dy + hr/24.0
    ra_sgn, ra_hr, ra_min, ra_sec = self.deg2dms(raDeg/15.0)
    dec_sgn, dec_deg, dec_min, dec_sec = self.deg2dms(decDeg)
    filt_trunc = filt.split('.')[0]
    if disc:
        disc_ast = "*"
    else:
        disc_ast = " "
    note1 = " "
    note2 = "C"
    
    mpc_format1 = "%5s%7s%1s%1s%1s%4d %02d %08.5lf "
    result_string = mpc_format1 % (mpNum, trackID, disc_ast, note1, note2, yr, mn, day)
    
    mpc_format2 = "%02d %02d %06.3lf%1s%02d %02d %05.2lf"
    result_string += mpc_format2 % (ra_hr, ra_min, ra_sec, dec_sgn, dec_deg, dec_min, dec_sec)
    
    #filt_trunc = ' '
    obsCode = '566'
    mpc_format3 = "       %6.1lf %1s      %3s"
    result_string += mpc_format3 % (mag, filt_trunc, obsCode)
    
    #result_string += "  %1s %1s" % (vac, rej)
    
    #print trackID, yr, mn, day, ra_hr, ra_min, ra_sec, dec_deg, dec_min, dec_sec, mag, filt_trunc, obsCode
    
    return result_string, int(vac), int(rej)

def H_alpha(H, G, alpha):
    # H is the absolute magnitude 
    # alpha is the solar phase angle in degrees
    A_1 = 3.332
    B_1 = 0.631
    C_1 = 0.986
    A_2 = 1.862
    B_2 = 1.218
    C_2 = 0.238

    alp = alpha*np.pi/180.0
    ta = np.tan(0.5*alp)
    sa = np.sin(alp)

    W = np.exp(-90.56*ta*ta)

    Phi_1_s = 1.0 - C_1 * np.sin(alp)/(0.119 + 1.341*sa - 0.754*sa*sa)
    Phi_1_l = np.exp(-A_1*np.power(ta, B_1))

    Phi_2_s = 1.0 - C_2 * np.sin(alp)/(0.119 + 1.341*sa - 0.754*sa*sa)
    Phi_2_l = np.exp(-A_2*np.power(ta, B_2))

    Phi_1 = W*Phi_1_s + (1.0-W)*Phi_1_l
    Phi_2 = W*Phi_2_s + (1.0-W)*Phi_2_l

    Ha = H - 2.5*np.log10((1.0-G)*Phi_1 + G*Phi_2)
  
    return Ha

Observatories = Observatory()
EOP = EarthAndTime()
