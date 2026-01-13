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
from astropy import coordinates as coord

from spacetrack import SpaceTrackClient #this pulls the TLE from a satellite database

import pandas as pd
 
        
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
def astronomy(my_locs, target_name = 'moon', timespacing = timedelta(seconds=1), timestart = datetime.now(timezone.utc), timeend = datetime.now(timezone.utc) + timedelta(seconds=10)):
    # import config # https://docs.python.org/3/faq/programming.html#how-do-i-share-global-variables-across-modules

    # import astropy.units as u
    # from astropy.time import Time
    # from datetime import datetime, date, timezone, timedelta
    # from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body, get_sun
    # from pygeomag import GeoMag
    
    print('Starting astronomy thread, target name is ' + target_name)
    print('Start time: ' + str(timestart))
    print('End time:   ' + str(timeend))
    
    target_time_az_el_list = []
    
    compass_dir = my_locs['compass_dir']   
    mag_declination = my_locs['mag_declination']
    my_loc = my_locs['EarthLocation']
    
    x = 1
    #time_old = datetime.now(timezone.utc)
    
    #locate sky target
    try:
        target = SkyCoord.from_name(target_name) #sky objects (outside of solar system)
    except Exception as e:
        print('could not find target, trying body search')
        try:
            target = get_body(target_name, Time(timestart), location = my_loc) #in solar system bodies
        except Exception as e:
            print('could not find target or body')
            print(e)
            return ['Astronomy object not found']
#            return ['invalid', 0, 0]
    #moon = get_body("moon", Time(telescope_time), location = my_loc)
    #sun = get_body("sun", Time(telescope_time), location = my_loc)
    
    #     if     datetime.now(timezone.utc) > time_old + timedelta(0,5): #rate limiter, remove after testing
    
    #https://docs.astropy.org/en/latest/coordinates/example_gallery_plot_obs_planning.html
    telescope_time = timestart
    target_time_az_el_list = []
    while telescope_time < timeend:
        
        #end input
        
        #need to convert moon to SkyCoord before proceeding to get alt/az
        
        #altitude:  angle between the object and the observer's local horizon (0-90 deg)
        #aziumuth: measured from true north and increasing eastward (N = 0, E = 90, S= 180, W = 270)
        
        targetaltaz = target.transform_to(AltAz(obstime=telescope_time, location=my_loc))
        #print(target_name + f"'s Altitude = {targetaltaz.alt:.2}")
        #print(target_name + f"'s Azimuth = {targetaltaz.az:.2}")
        
        #adjust azimuth for position of telescope base
        az = (targetaltaz.az.deg + compass_dir + mag_declination.d ) % 360
        el = targetaltaz.alt.deg
        
        print('az: ' + str(az) + ' ' + 'el: ' + str(el))
        x = x + 1
        target_time_az_el_list.append([target_name, telescope_time, az, el, 0, 0]) #trailing zeros are the placeholder velocities, computed outside of this function
        
        telescope_time = telescope_time + timespacing#Time.now() #datetime.now(timezone.utc) #have to make sure datetime is in utc for all the astro tools unless you specify it in them indivudally

    return target_time_az_el_list
    print('Terminated astronomy thread')

def satellite_tracking(my_locs, target_name = 'ISS', timespacing = timedelta(seconds=1), timestart = datetime.now(timezone.utc), timeend = datetime.now(timezone.utc) + timedelta(seconds=10)):

    print('Starting satellite_tracking, target is ' +target_name)
    #if this is too slow, make a table to timetagged angles to run through instead of trying to calculate in realtime?
    
    # import config # https://docs.python.org/3/faq/programming.html#how-do-i-share-global-variables-across-modules

    # from sgp4.api import Satrec
    # from sgp4.api import SGP4_ERRORS
    
    # from datetime import datetime, date, timezone, timedelta
    # from astropy.time import Time, TimeDelta, TimeDeltaSec

    #load local copy of NORAD satellite list and find NORAD catid by name
    norad = pd.read_csv(r'/Users/andrewmiller/telescope/satcat.csv')
    #if target_name is only digits assume this is the NORAD catalog ID
    target_name_original = target_name
    target_name = target_name.upper() #NORAD catalog is allcaps
    print('Desired target name is ' + target_name)
    if target_name.isdigit():
        norad_cat_id = int(target_name)
        print('Target name from NORAD catalog ID is ')
        print(norad[['OBJECT_NAME', 'OBJECT_ID', 'NORAD_CAT_ID']].loc[norad['NORAD_CAT_ID']==(int(target_name))].to_string(index=False))

    else:
        #check for exact match
        name_matches_full = norad[norad['OBJECT_NAME'].str.fullmatch(target_name)]
        print('Found ' + str(len(name_matches_full)) + ' name_matches_full for ' + target_name)
        print(name_matches_full)
        #check for partial match        
        name_matches = norad[norad['OBJECT_NAME'].str.contains(target_name)]
        print('Found ' + str(len(name_matches)) + ' name_matches for ' + target_name)
        print(name_matches)
        #if an exact match, use NORAD (search for 'moonlighter')
        #if a partial match but not an exact match use NORAD (search for 'moo', get moonlighter)
        #if neither use astronomy (search for 'moon')
        if len(name_matches_full) == 1:
            print(name_matches_full['NORAD_CAT_ID'])
            norad_cat_id = int(name_matches_full['NORAD_CAT_ID'].iloc[0])
            #only a single match, so return that
        elif len(name_matches_full) < 1: #see if a full name match exists for the astronomy targets
            #use NORAD
            print('Could not find object. May be decayed or does not exist. If a celestial body, use astronomy mode instead of satellite tracking.')
            target_time_az_el_list = astronomy(my_locs, target_name_original, timespacing, timestart, timeend) #presume it's a celestial or solar system body
            if target_time_az_el_list[0] == 'Astronomy object not found':
                print('Could not find body in astronomy, putting in partial search of NORAD database')
            else:
                return target_time_az_el_list #succesfully got the values out of astronomy
        
        if (len(name_matches) >= 1): #multiple partial matches, don't knwo which one so display everything that isn't decayed
            nondecay_matcheds = name_matches[~name_matches['OPS_STATUS_CODE'].str.contains('D')]
            if len(nondecay_matcheds.index)>1:
                print('Found multiple matches for target \'' + target_name + '\', please specify which and re-run.')
                nondecay_matcheds = nondecay_matcheds[['OBJECT_NAME', 'OBJECT_ID', 'NORAD_CAT_ID']].reset_index(drop = True)
                print(nondecay_matcheds.to_string(index=True))
                val = input('Please enter desired index: ')
                print('You have selected ')
                print(nondecay_matcheds.iloc[int(val)].to_string())
                norad_cat_id = int(nondecay_matcheds['NORAD_CAT_ID'].iloc[int(val)]) #int(ser.iloc[0])
            else:
                print(nondecay_matcheds['NORAD_CAT_ID'])
                norad_cat_id = int(nondecay_matcheds['NORAD_CAT_ID'].iloc[0])
        else:
            print('unknown object, did not find in astronomy or satellite tracking')
    print('rolling with norad_cat_id = ' + str(norad_cat_id))
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
                
        # itrs_geo = teme.transform_to(ITRS(obstime=t))
        # location = itrs_geo.earth_location
        # location.geodetic
        
        
        # from astropy.coordinates import ITRS
        itrs_geo = teme.transform_to(ITRS(obstime=t)) #location of object over earth in lat/long
        location = itrs_geo.earth_location
        #print(location.geodetic)
        
        # from astropy.coordinates import EarthLocation, AltAz
        #siding_spring = EarthLocation.of_site('aao')
        obs_location = my_locs['EarthLocation']
    
        topo_itrs_repr = itrs_geo.cartesian.without_differentials() - obs_location.get_itrs(t).cartesian
        itrs_topo = ITRS(topo_itrs_repr, obstime = t, location=obs_location)
        aa = itrs_topo.transform_to(AltAz(obstime=t, location=obs_location))
        return aa

    def get_lat_long(t): #report satellite location as lat/long
        t.format = 'jd' #SGP4 input is julian days and fractional days
        error_code, teme_p, teme_v = satellite.sgp4(t.jd1, t.jd2)  # in km and km/s
        if error_code != 0:
            raise RuntimeError(SGP4_ERRORS[error_code])
            
        teme_p = CartesianRepresentation(teme_p*u.km)
        teme_v = CartesianDifferential(teme_v*u.km/u.s)
        teme = TEME(teme_p.with_differentials(teme_v), obstime=t)
                
        itrs_geo = teme.transform_to(ITRS(obstime=t))
        location = itrs_geo.earth_location
        return location.geodetic


    # from spacetrack import SpaceTrackClient #this pulls the TLE from a satellite database
    try:
        print('Connecting to Spacetrack')
        st = SpaceTrackClient('andymiller@gmail.com', 'spacetrackpassword')
        print('Spacetrack TLE: ')
        #print(st.gp(norad_cat_id=[25544, 41335], format='tle')) #can get more than one TLE but not doing that here
        x = st.gp(norad_cat_id=[norad_cat_id], format='tle')
        tle1 = x.splitlines()[0] #first line of two line ephemeris (TLE)
        tle2 = x.splitlines()[1] #second line of two line ephemeris (TLE)
        print(x)
    except Exception as e:
        print(e)
        print('Error loading object in SpaceTrack!')
        #return
#        TLE data for most spacecraft can be downloaded from https://www.space-track.org/
#       Tracking map for TLE lookup: https://in-the-sky.org/satmap_radar.php?town=5809844
# set time into the future, look for in view target, pull its info (implement lookup for single number input to grab TLE?), arm telescope which should sit until sat emerges above horizon



    
    satellite = Satrec.twoline2rv(tle1, tle2) #this converts the TLE into a track that can get fed into astropy to produce the position at a time and then get angles for that
    print('Satellite track computed, comparing in view times')
    #see if astroplan works better for this than DIY
    #get first point above horizon
        #advance time until it is above, then binary search between last time below and then?
    #report time to that point
    #move scope to that point
    #wait for in view and then start tracking (include countdown)
    
    
    ########## IN WORK ##########
    def get_inview_times_TLL():
        #direct calculation approach: https://digitalcommons.usu.edu/cgi/viewcontent.cgi?referer=&httpsredir=1&article=2075&context=smallsat
        #using 2.1 to find approximate close locations, then use SGP4 directly to get details
        lat = config.lat
        long = config.long
        n = float(tle2.split()[7]) #orbital mean motion
        orbital_period = 1/(n/(60*60*24) )#orbits per day to 1/orbits per second = seconds per orbit
        sat_latlong = get_lat_long(Time.now())
        sat_lat = sat_latlong.lat.deg #arbitrary starting point in the past
        #find time where satellite crosses ground station latitude
        t = Time.now() - TimeDelta(24*3600, format='sec') #start 1 day ago so we can always be moving forward in time in the calculations
        for secs in [1000, 100, 10, 1, 0.1]: #find the time of latitude crossing
            sat_lat = get_lat_long(t).lat.deg
            while (sat_lat <= lat): #get rise time
                t = t + TimeDelta(secs, format='sec')
                sat_lat = get_lat_long(t).lat.deg #the found crossing point
                print('TLL at ' + str(sat_lat))
            t = t - TimeDelta(secs, format='sec') #go back in time to the last position it was still up
        t.format='isot'
        print('previous latitude crossing time:')
        print(t)
        #find time of next TLL crossings and get latitude there (orbital period)
        for x in range(10):
            t = t + TimeDelta(x*orbital_period, format='sec')
            sat_lat = get_lat_long(t).lat.deg #the found crossing point
            sat_long = get_lat_long(t).lon.deg
            t.format='isot'
            print('TLL at ' + str("%.2f" % sat_lat) + ' Long error ' + str("%.2f" % (long-sat_long)) + ' ' + t.value)

        #check that the math is right by brute force
        t = Time.now() - TimeDelta(24*3600, format='sec') #start 1 day ago so we can always be moving forward in time in the calculations
        biglist = []
        for x in range(2592): #3 days of 100 seconds
            t.format='isot'
            t2 = t + TimeDelta(100*x, format='sec')
            aa = get_az_el(t2)
            t2.format='isot'
            biglist.append([aa.alt.deg, t2])
            print(str(x) + ' ' + t2.value)
        
        for element in biglist:
            if element[0] > 0:
                print(element)
                break
            print(element.el)

        #now make list of crossing times from yesterday's crossing until the next one that is above the horizon
        earth_rotation_rate = (360.9856)/(24*3600) #deg/day into deg/sec
        #longitudinal change at same latitude
        delta_long = earth_rotation_rate * orbital_period #degrees/sec * sec, change in longitude per orbit
        
        
        #t = t + TimeDelta(, format='sec')
        
        deltaV = get_lat_long(t).lon.deg - config.long
        
        
        
        N = (deltaV*n)/(360*earth_rotation_rate) #use degrees instead of radians in the 2.1 equations. Result is seconds to next TLL (target latitude location)
        for x in [1,2,3,4,5,6]:
            N = x*(deltaV*n)/(360*earth_rotation_rate)
            
    
    ########## END IN WORK ##########

    
    
    #get the in view times
    t = Time.now()
    t.format = 'isot'
    print('Current time:  ' + t.value)

    aa = get_az_el(t)
    if aa.alt<0: #object is currently below the horizon
        print('Currently set')
        
        t = Time.now()
        for secs in [1000, 100, 10, 1, 0.1]: #find the object's rise time
        #this is crude and can skip over small windows
        #direct calculation approach: https://digitalcommons.usu.edu/cgi/viewcontent.cgi?referer=&httpsredir=1&article=2075&context=smallsat
        #get times of closest approach from that
        #put them into Satrec and see if above horizon
        #if no, repeeat
#            print('steps of ' + str(secs) +' seconds')
            aa = get_az_el(t)
            while (aa.alt <= 0): #get rise time
                t = t + TimeDelta(secs, format='sec')
                aa = get_az_el(t)
#                print(aa.alt)
            t = t - TimeDelta(secs, format='sec') #go back in time to the last position it was still up
        rise_time = t
        t.format = 'isot'
        print('Rise time:     ' + t.value)
        
        #get the next set time for doing the trajectory calculations
        t = rise_time + TimeDelta(10, format='sec') #shift so don't accidentally hit the transition time
        for secs in [1000, 100, 10, 1, 0.1]: #find the object's rise time, this probably needs optimization to avoid edge cases?
#            print('steps of ' + str(secs) +' seconds')
            aa = get_az_el(t)
            while (aa.alt >= 0): #get set time
                t = t + TimeDelta(secs, format='sec')
                aa = get_az_el(t)
#                print(aa.alt)
            t = t - TimeDelta(secs, format='sec') #go back in time to the last position it was still up
        set_time = t
        t.format='isot'
        print('Next set time: ' + t.value)

        
        countdown = t-Time.now()
        # while countdown>0:
        countdown = t-Time.now()
        countdown.format = 'sec'
        countdown = countdown/60
        print(str("%.2f" % countdown.value) + ' minutes until risen')
            
    else: #object is above the horizon, get the time it sets
        rise_time = Time.now() #can't do stuff in the past so get calculations from current moment
        print('currently above the horizon')
        for secs in [1000, 100, 10, 1, 0.1]:
#            print('steps of ' + str(secs) +' seconds')
            aa = get_az_el(t)
            while (aa.alt >= 0): #get set time
                t = t + TimeDelta(secs, format='sec')
                aa = get_az_el(t)
#                print(aa.alt)
            t = t - TimeDelta(secs, format='sec') #go back in time to the last position it was still up
        set_time = t
        t.format = 'isot'
        print('Set time: ' + str(t.value))

        countdown = t-Time.now()
        countdown.format = 'sec'
        countdown = countdown/60
        print(str("%.2f" % countdown.value) + ' minutes until set')    
##### End in view time calculations
    rise_time = rise_time.to_datetime() #astropy Time class slow, avoid using
    set_time = set_time.to_datetime()

    rise_time = rise_time.replace(tzinfo=timezone.utc) #make UTC aware
    set_time = set_time.replace(tzinfo=timezone.utc) #make UTC aware
  
######Bug noticed: when re-running without restarting kernel the in view times report the n+1 rise window rather than the next one
#this might be if the rise is <1000 seconds away? Or if the 1000 sec step passes over an inview period based on its start point (I think this might be it)
# if time in view is in 1100-1300 seconds and we check at 0, 1000, 2000 seconds it will clip it

    
    ####### add yes/no check if in view time not in next 5 minutes (wait or try new target)
    config.coordinates = []
    t = rise_time
    my_loc = config.my_locs['EarthLocation']
    compass_dir = config.my_locs['compass_dir']   
    mag_declination = config.my_locs['mag_declination']

    coordinates = []
    target_time_az_el_list = []
    print('start getting az el coordinates to return')
    print('rise time         : ' + str(rise_time))
    print('set time          : ' + str(set_time))
    print('starting time     : ' + str(t))
    print('user ending time  : ' + str(timeend))

    print('pulling coordinate at time ' + str(t))
    aa = get_az_el(Time(t)) #get the coordinates        
    
#moved sun checks out of here        #sun = get_body("sun", Time(t), location = my_loc)
    #sunaltaz = sun.transform_to(AltAz(obstime=t, location=my_loc))
    #sun_az = sunaltaz.az.deg + compass_dir + mag_declination.d
    #sun_el = sunaltaz.alt.deg
    #coordinates.append({'az':aa.az.deg, 'el':aa.alt.deg, 'sun_az':sun_az, 'sun_el':sun_el, 'time':t})
    t = timestart #ignore the in view times for now, just pull coordinates for where it is at the specified time
    while t < timeend:
        aa = get_az_el(Time(t)) #get the coordinates        

        target_time_az_el_list.append([target_name_original, t, aa.az.deg, aa.alt.deg, 0, 0]) #the trailing 0 s are velocity placeholders, to be computed outside of this function
        print([target_name_original, t, aa.az.deg, aa.alt.deg])
        timespacing = timedelta(seconds=0.1)
        t = t + timespacing
        #timespacing might need to be dynamic, adjust until the minimum spacing between angles is the encoder resolution minimum angle or some small enough value
        #might also need to try to give it the next coodinate before it reaches the correct one or something to let it keep rolling smoothly
        
#    config.coordinates = coordinates #store in shared variable


#     while t < set_time: #compute the telescope time/az/el sets
#         #get az-el and time
#         if timeend is not None:
#             if t < timeend:
#                 pass #wait until time is in the specified zone
#             else:
#                 t = t + timespacing
#                 continue

#         else: #before the end time or no end time specified
#             print('pulling coordinate at time ' + str(t))
#             aa = get_az_el(Time(t)) #get the coordinates        
            
#     #moved sun checks out of here        #sun = get_body("sun", Time(t), location = my_loc)
#             #sunaltaz = sun.transform_to(AltAz(obstime=t, location=my_loc))
#             #sun_az = sunaltaz.az.deg + compass_dir + mag_declination.d
#             #sun_el = sunaltaz.alt.deg
#             #coordinates.append({'az':aa.az.deg, 'el':aa.alt.deg, 'sun_az':sun_az, 'sun_el':sun_el, 'time':t})
#             target_time_az_el_list.append([target_name, t, aa.az.deg, aa.alt.deg])

#             print([target_name, t, aa.az.deg, aa.alt.deg])
#             t = t + timespacing
# #    config.coordinates = coordinates #store in shared variable

    print('Completed trajectory calculations')
    mins= int((set_time - rise_time).seconds/60)
    secs = (set_time - rise_time).seconds
    print('In view time: ' + str(secs) +' seconds (' + str(mins) + ' min ' + str(secs-(60*mins)) + ' sec)')
    print(target_time_az_el_list[0])
    return target_time_az_el_list
    

#%% extra functionality below
def do_fancy_things():

    #check that the sun won't end up in the trajectory
    keep_out = []
    keep_out_angle = 5
    print('Check for sun angle intrusions and write to file for reference')
    with open(r'/Users/andrewmiller/telescope/track.txt', 'w') as f:

        for line in config.coordinates:
            f.write("%s\n" % line)     #save to file for review
            if abs(line['az']-line['sun_az'])<keep_out_angle and abs(line['el']-line['sun_el'])<keep_out_angle:
                print('Warning: trajectory approaches sun keep out: sun at ' + str(line['sun_az']) + ',' + str(line['sun_el']) + ' track at ' + str(line['az']) +','+str(line['el']))
                keep_out.append(line)
        print('Done writing to file')

    if len(keep_out) > 0:
        print('Halting due to intrusion into keep-out zone')
        #return
    else:
        config.tracking_ready = True # signal that the tracker is now producing tracking coordinates



#check for overspeeds (likely near keyhole events)
#coord(N+1) - coord(N) /time(N+1)-time(N) > system deg/sec capability
#probably try anyway as could be interesting and don't want to discard stuff that flies straight overhead in any case

# #Realtime        
#     time_old = datetime.now()
#     time_print = datetime.now()
    
#     while not config.end_program:
#         config.tracking_ready = True # signal that the tracker is now producing tracking coordinates
        
#         timing.append(datetime.now())
#         aa = get_az_el(Time.now())

#         if aa.alt > 0:
#             time1 = datetime.now()

#             countdown = set_time - datetime.now()
#             if datetime.now() > time_print + timedelta(seconds=10):
#                 print(str(countdown.seconds) + '.' + str(countdown.microseconds) + 's until out of view')
#                 print('current tracker generated angle az: ' + str(aa.az.deg) + ' el: ' + str(aa.alt.deg))

#                 time_print = datetime.now()

            
# #            print('ITRS: ' + str(Time.now()) + ' ' + str(aa.alt) + ' ' + str(aa.az))
#             #print('ITRS: ' + str(Time.now()) + ' ' + str(aa.alt.deg) + ' ' + str(aa.az.deg))
            
#             time1 = datetime.now()
#             send_commands(my_locs, aa.az.deg, aa.alt.deg, queues['telescope_q'])
#             #don't use astropy Time print((datetime.now() - time1).quantity_str + ' send commands')
# #            print(str((datetime.now() - time1).seconds) + '.' + str((datetime.now() - time1).microseconds) + 's send commands')
#             #print('deg az: ' + str(aa.az.deg) + ' deg el ' + str(aa.alt.deg))

#             if datetime.now() > time_old + timedelta(seconds=camera_period): #check if time to take a photo
#                 queues['camera_q'].put('az ' + str(aa.az.deg) + ' ' + 'el ' + str(aa.alt.deg),block=False)#take photo
#                 time_old = datetime.now()
#             #don't use astropy Time print((datetime.now() - time1).quantity_str + ' time to send')
# #            print(str((datetime.now() - time1).seconds) + '.' + str((datetime.now() - time1).microseconds) + 's time to send')

#         else:
#             countdown = datetime.now() - rise_time
#             #don't use astropy Time print(countdown.quantity_str + ' until in view')
#             print(str(countdown.seconds) + '.' + str(countdown.microseconds) + 's until in view')

#             pass
#     #print(timing)
#     config.x = timing
# #End realtime

    
    print('Terminated satellite tracking thread')

#use the above functions to build list of az/el commands
#then use the send_command to do the linking through lists and just update a target az/el variable in config that gets read by telescope control?
#send_commands does the real-time sun pointing safety checks but needs to do it fast enough

#needs fixing: doesn't check for sun during initial slews, maybe break those up into n degree steps?



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