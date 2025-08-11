# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 10:52:25 2025

@author: AndrewMiller
"""
import config # https://docs.python.org/3/faq/programming.html#how-do-i-share-global-variables-across-modules

# import astropy.units as u
from datetime import datetime, date, timezone, timedelta
from pygeomag import GeoMag

import time    

from sgp4.api import Satrec
from sgp4.api import SGP4_ERRORS

from astropy.time import Time, TimeDelta, TimeDeltaSec
from astropy import units as u
from astropy.coordinates import EarthLocation, AltAz, ITRS, TEME, CartesianDifferential, CartesianRepresentation,  SkyCoord, get_body, get_sun
from spacetrack import SpaceTrackClient #this pulls the TLE from a satellite database

 
        
#point and shoot mode:
#aim at [center or top left] and raster scan as described
def point_and_shoot(my_loc, queues, FOV = 1, az_steps = 10, el_steps = 5): #units degrees

    config.tracking_ready = True # signal that the tracker is now producing tracking coordinates

    
    print('Starting point and shoot thread')

    #compute steps as degrees to send

    for el in range (0,el_steps): #ADD: need to relate steps to FOV for decent coverage with overlaps for stitching
        for az in range(0,az_steps):
            send_commands(my_loc, az, el, queues['telescope_q'])
            
            while not config.az_in_position and not config.el_in_position:
                pass #wait until both axes in position
                
            time.sleep(1) # Sleep to let jitter settle
            #add photo request to camera queue
            queues['camera_q'].put('point_and_shoot_'+str(az)+'_'+str(el))
            #wait until queue size is zero (photo taken) before moving on
            while not config.camera_ready:
                pass #wait for camera to finish
            
    print('Completed point and shoot thread')
    
#include way to scan by start/end lat/long coordinates

#astronomy mode:
def astronomy(my_locs, queues, target_name = 'moon', camera_period = 10):
    # import config # https://docs.python.org/3/faq/programming.html#how-do-i-share-global-variables-across-modules

    # import astropy.units as u
    # from astropy.time import Time
    # from datetime import datetime, date, timezone, timedelta
    # from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body, get_sun
    # from pygeomag import GeoMag

    print('Starting astronomy thread')

    
    compass_dir = my_locs['compass_dir']   
    mag_declination = my_locs['mag_declination']
    my_loc = my_locs['EarthLocation']

    x = 1
    time_old = datetime.now(timezone.utc)
    
    #locate sky target
    try:
        target = SkyCoord.from_name(target_name) #sky objects (outside of solar system)
    except Exception as e:
        print('could not find target, trying body search')
        try:
            target = get_body(target_name, Time(time_old), location = my_loc) #in solar system bodies
        except Exception as e:
            print('could not find target or body')
            print(e)
            config.end_program = True
    #moon = get_body("moon", Time(telescope_time), location = my_loc)
    #sun = get_body("sun", Time(telescope_time), location = my_loc)

    
    
    while not config.end_program:
            config.tracking_ready = True # signal that the tracker is now producing tracking coordinates

   #     if     datetime.now(timezone.utc) > time_old + timedelta(0,5): #rate limiter, remove after testing

            #https://docs.astropy.org/en/latest/coordinates/example_gallery_plot_obs_planning.html
            telescope_time = Time.now() #datetime.now(timezone.utc) #have to make sure datetime is in utc for all the astro tools unless you specify it in them indivudally
            #end input

            #need to convert moon to SkyCoord before proceeding to get alt/az
            
            #altitude:  angle between the object and the observer's local horizon (0-90 deg)
            #aziumuth: measured from true north and increasing eastward (N = 0, E = 90, S= 180, W = 270)
    
            targetaltaz = target.transform_to(AltAz(obstime=telescope_time, location=my_loc))
            #print(target_name + f"'s Altitude = {targetaltaz.alt:.2}")
            #print(target_name + f"'s Azimuth = {targetaltaz.az:.2}")
            
            #adjust azimuth for position of telescope base
            az = targetaltaz.az.deg + compass_dir + mag_declination.d
            el = targetaltaz.alt.deg
            
            print('az: ' + str(az) + ' ' + 'el: ' + str(el))
            x = x + 1
            send_commands(my_locs, az, el, queues['telescope_q'])
            
            if datetime.now(timezone.utc) > time_old + timedelta(seconds=camera_period): #check if time to take a photo
                queues['camera_q'].put('az ' + str(az) + ' ' + 'el ' + str(el))#take photo
                time_old = datetime.now(timezone.utc)
        
    #print('el: ' + str(el))
    print('Terminated astronomy thread')

def satellite_tracking(my_locs, queues, target_name = 'ISS', camera_period = 1):
    
    #if this is too slow, make a table to timetagged angles to run through instead of trying to calculate in realtime?
    
    # import config # https://docs.python.org/3/faq/programming.html#how-do-i-share-global-variables-across-modules

    # from sgp4.api import Satrec
    # from sgp4.api import SGP4_ERRORS
    
    # from datetime import datetime, date, timezone, timedelta
    # from astropy.time import Time, TimeDelta, TimeDeltaSec
    timing = []

        
    def get_az_el(t):
    
        t.format = 'jd' #SGP4 input is julian days and fractional days
        error_code, teme_p, teme_v = satellite.sgp4(t.jd1, t.jd2)  # in km and km/s
        if error_code != 0:
            raise RuntimeError(SGP4_ERRORS[error_code])
            
        # from astropy.coordinates import TEME, CartesianDifferential, CartesianRepresentation
        # from astropy import units as u
        teme_p = CartesianRepresentation(teme_p*u.km)
        teme_v = CartesianDifferential(teme_v*u.km/u.s)
        teme = TEME(teme_p.with_differentials(teme_v), obstime=t)
                
        # from astropy.coordinates import ITRS
        itrs_geo = teme.transform_to(ITRS(obstime=t))
        location = itrs_geo.earth_location
        
        # from astropy.coordinates import EarthLocation, AltAz
        siding_spring = EarthLocation.of_site('aao')
        siding_spring = my_locs['EarthLocation']
    
        topo_itrs_repr = itrs_geo.cartesian.without_differentials() - siding_spring.get_itrs(t).cartesian
        itrs_topo = ITRS(topo_itrs_repr, obstime = t, location=siding_spring)
        aa = itrs_topo.transform_to(AltAz(obstime=t, location=siding_spring))
        return aa


    # from spacetrack import SpaceTrackClient #this pulls the TLE from a satellite database
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
    
    #get the in view times
    t = Time.now()
    for secs in [1000, 100, 10, 1]:
        print('steps of ' + str(secs) +' seconds')
        aa = get_az_el(t)
        while (aa.alt >= 0): #get set time
            t = t + TimeDelta(secs, format='sec')
            aa = get_az_el(t)
            print(aa.alt)
        t = t - TimeDelta(secs, format='sec') #go back in time to the last position it was still up
    set_time = t
    t.format = 'isot'
    print('Set time: ' + t.value)
    set_time = set_time.datetime #astropy Time class slow, avoid using
        
    t = Time.now()
    for secs in [1000, 100, 10, 1]:
        print('steps of ' + str(secs) +' seconds')
        aa = get_az_el(t)
        while (aa.alt <= 0): #get set time
            t = t - TimeDelta(secs, format='sec')
            aa = get_az_el(t)
            print(aa.alt)
        t = t - TimeDelta(secs, format='sec') #go back in time to the last position it was still up
    rise_time = t
    t.format = 'isot'
    print('Rise time: ' + t.value)  
    rise_time = rise_time.datetime #astropy Time class slow, avoid using

    
        
    time_old = datetime.now()
    time_print = datetime.now()
    
    while not config.end_program:
        config.tracking_ready = True # signal that the tracker is now producing tracking coordinates
        
        timing.append(datetime.now())
        aa = get_az_el(Time.now())

        if aa.alt > 0:
            time1 = datetime.now()

            countdown = set_time - datetime.now()
            if datetime.now() > time_print + timedelta(seconds=10):
                print(str(countdown.seconds) + '.' + str(countdown.microseconds) + 's until out of view')
                print('current tracker generated angle az: ' + str(aa.az.deg) + ' el: ' + str(aa.alt.deg))

                time_print = datetime.now()

            
#            print('ITRS: ' + str(Time.now()) + ' ' + str(aa.alt) + ' ' + str(aa.az))
            #print('ITRS: ' + str(Time.now()) + ' ' + str(aa.alt.deg) + ' ' + str(aa.az.deg))
            
            time1 = datetime.now()
            send_commands(my_locs, aa.az.deg, aa.alt.deg, queues['telescope_q'])
            #don't use astropy Time print((datetime.now() - time1).quantity_str + ' send commands')
#            print(str((datetime.now() - time1).seconds) + '.' + str((datetime.now() - time1).microseconds) + 's send commands')
            #print('deg az: ' + str(aa.az.deg) + ' deg el ' + str(aa.alt.deg))

            if datetime.now() > time_old + timedelta(seconds=camera_period): #check if time to take a photo
                queues['camera_q'].put('az ' + str(aa.az.deg) + ' ' + 'el ' + str(aa.alt.deg),block=False)#take photo
                time_old = datetime.now()
            #don't use astropy Time print((datetime.now() - time1).quantity_str + ' time to send')
#            print(str((datetime.now() - time1).seconds) + '.' + str((datetime.now() - time1).microseconds) + 's time to send')

        else:
            countdown = datetime.now() - rise_time
            #don't use astropy Time print(countdown.quantity_str + ' until in view')
            print(str(countdown.seconds) + '.' + str(countdown.microseconds) + 's until in view')

            pass
    #print(timing)
    config.x = timing
    print('Terminated satellite tracking thread')


#Keep out zones (after getting desired commanded position)
def send_commands(my_locs, az, el, telescope_q):
   #  import config # https://docs.python.org/3/faq/programming.html#how-do-i-share-global-variables-across-modules

   # # import astropy.units as u
   #  from astropy.time import Time
   #  from datetime import datetime, date, timezone, timedelta
   #  from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body, get_sun
   #  from pygeomag import GeoMag

    my_loc = my_locs['EarthLocation']
    compass_dir = my_locs['compass_dir']   
    mag_declination = my_locs['mag_declination']
        
    #Check if will point into ground
    send_command = True
    if el < 0:
        send_command = False
        print("Keep out: pointing into ground")
        
    #get sun location
    telescope_time = datetime.now(timezone.utc) #have to make sure datetime is in utc for all the astro tools unless you specify it in them indivudally
    sun = get_body("sun", Time(telescope_time), location = my_loc)
    sunaltaz = sun.transform_to(AltAz(obstime=telescope_time, location=my_loc))
    sun_az = sunaltaz.az.deg + compass_dir + mag_declination.d
    sun_el = sunaltaz.alt.deg
    keep_out_radius = 5 #5 degree keep out radius to start
    
    config.az_angle = az
    config.el_angle = el
    
    #Check if will point into sun
    if (az-sun_az)**2 + (el - sun_el)**2 < keep_out_radius**2:
        send_command = False
        print("Keep out: pointing at sun")
    
    
    if send_command:
        telescope_q.put({'az':az, 'el':el, 'cmd':'move_az_el'}, block=False)
        config.coordinates.append({'time':telescope_time, 'az':az, 'el':el, 'cmd':'move_az_el'})
        #print('generated new angle az: ' + str(az) + ' el: ' + str(el))
    
    # config.x = config.x + 1
    # if config.x > 100:
    #     config.end_program = True
        

#tracking mode:
#def tracking:
    
    #realtime tracking of object in view at current moment
    #identify with selection or just pick thing in the center