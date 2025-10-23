# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 13:25:41 2025

@author: AndrewMiller
"""

import config # https://docs.python.org/3/faq/programming.html#how-do-i-share-global-variables-across-modules

from sgp4.api import Satrec
from sgp4.api import SGP4_ERRORS

from datetime import datetime, date, timezone, timedelta

from spacetrack import SpaceTrackClient #this pulls the TLE from a satellite database
st = SpaceTrackClient('andymiller@gmail.com', 'spacetrackpassword')
print('Spacetrack: ')
print(st.gp(norad_cat_id=[25544, 41335], format='tle'))
#        TLE data for most spacecraft can be downloaded from https://www.space-track.org/
#       Tracking map for TLE lookup: https://in-the-sky.org/satmap_radar.php?town=5809844
# set time into the future, look for in view target, pull its info (implement lookup for single number input to grab TLE?), arm telescope which should sit until sat emerges above horizon
s = '1 25544U 98067A   19343.69339541  .00001764  00000-0  38792-4 0  9991' #first line of two line ephemeris (TLE)
t = '2 25544  51.6439 211.2001 0007417  17.6667  85.6398 15.50103472202482' #this is the ISS
s = '1 52335U 22045E   25168.14995947  .00001137  00000-0  89539-4 0  9991'
t = '2 52335  53.2181 321.1759 0001263  88.2138 271.8998 15.08840359175550'
s='1 52643U 22052AX  25167.82616276  .00002229  00000-0  15815-3 0  9998'
t='2 52643  53.2186  12.5999 0001438  95.1413 264.9743 15.08832272170990'

s='1 53927U 22119AT  25168.31106461  .00085274  00000-0  14213-2 0  9994'
t='2 53927  53.1931  17.1406 0002427 115.7151 244.4110 15.51292049152018'
s = '1 36395U 10005A   25168.45492097 -.00000067  00000-0  00000+0 0  9997'
t = '2 36395  33.8852  91.9242 0000558  95.8239 140.0671  1.00270500 56373'

satellite = Satrec.twoline2rv(s, t) #this converts the TLE into a track that can get fed into astropy to produce the position at a time and then get angles for that

#see if astroplan works better for this than DIY
#get first point above horizon
    #advance time until it is above, then binary search between last time below and then?
#report time to that point
#move scope to that point
#wait for in view and then start tracking (include countdown)

    
#get if in view currently
#if yes, only find out of view end
#if no, start advancing with 1 min, 2 min, 4 min etc




from astropy.time import Time, TimeDelta, TimeDeltaSec

    
def get_az_el(t = Time.now()):

    t.format = 'jd' #SGP4 input is julian days and fractional days
    error_code, teme_p, teme_v = satellite.sgp4(t.jd1, t.jd2)  # in km and km/s
    if error_code != 0:
        raise RuntimeError(SGP4_ERRORS[error_code])
        
    from astropy.coordinates import TEME, CartesianDifferential, CartesianRepresentation
    from astropy import units as u
    teme_p = CartesianRepresentation(teme_p*u.km)
    teme_v = CartesianDifferential(teme_v*u.km/u.s)
    teme = TEME(teme_p.with_differentials(teme_v), obstime=t)
            
    from astropy.coordinates import ITRS
    itrs_geo = teme.transform_to(ITRS(obstime=t))
    location = itrs_geo.earth_location
    
    from astropy.coordinates import EarthLocation, AltAz
    siding_spring = EarthLocation.of_site('aao')
    siding_spring = my_locs['EarthLocation']

    topo_itrs_repr = itrs_geo.cartesian.without_differentials() - siding_spring.get_itrs(t).cartesian
    itrs_topo = ITRS(topo_itrs_repr, obstime = t, location=siding_spring)
    aa = itrs_topo.transform_to(AltAz(obstime=t, location=siding_spring))
    return aa



t = Time.now()
for secs in [1000, 100, 10, 1]:
    print('steps of ' + str(secs) +' seconds')
    aa = get_az_el(t)
    while (aa.alt >= 0): #get set time
        t = t + TimeDelta(secs, format='sec')
        aa = get_az_el(t)
        print(aa.alt)
    t = t - TimeDelta(secs, format='sec') #go back in time to the last position it was still up

t.format = 'isot'
print('Set time: ' + t.value)
    
t = Time.now()
for secs in [1000, 100, 10, 1]:
    print('steps of ' + str(secs) +' seconds')
    aa = get_az_el(t)
    while (aa.alt <= 0): #get set time
        t = t - TimeDelta(secs, format='sec')
        aa = get_az_el(t)
        print(aa.alt)
    t = t - TimeDelta(secs, format='sec') #go back in time to the last position it was still up
t.format = 'isot'
print('Rise time: ' + t.value)  
    

if aa.alt > 0:
#            print('ITRS: ' + str(Time.now()) + ' ' + str(aa.alt) + ' ' + str(aa.az))
    #print('ITRS: ' + str(Time.now()) + ' ' + str(aa.alt.deg) + ' ' + str(aa.az.deg))

    print('deg az: ' + str(aa.az.deg) + ' deg el ' + str(aa.alt.deg))



print('Terminated satellite tracking thread')
