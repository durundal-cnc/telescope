# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 10:52:25 2025

@author: AndrewMiller
"""


#point and shoot mode:
#aim at [center or top left] and raster scan as described
def point_and_shoot(my_loc, queues, FOV = 1, az_steps = 10, el_steps = 5): #units degrees
    import config # https://docs.python.org/3/faq/programming.html#how-do-i-share-global-variables-across-modules

    import time    
    
    print('Starting point and shoot thread')
    
    def in_position(target, actual, deadband = 1): #wait and see if settled for static shots
        if abs(target-actual) < deadband:
            return True
        else:
            return False
    
    #compute steps as degrees to send

    for el in range (0,el_steps): #ADD: need to relate steps to FOV for decent coverage with overlaps for stitching
        for az in range(0,az_steps):
            send_commands(my_loc, az, el, queues['telescope_q'])
            while not in_position(config.az_angle, az) and not in_position(config.el_angle, el):
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
    import config # https://docs.python.org/3/faq/programming.html#how-do-i-share-global-variables-across-modules

    import astropy.units as u
    from astropy.time import Time
    from datetime import datetime, date, timezone, timedelta
    from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body, get_sun
    from pygeomag import GeoMag

    print('Starting astronomy thread')

    
    compass_dir = my_locs['compass_dir']   
    mag_declination = my_locs['mag_declination']
    my_loc = my_locs['EarthLocation']

    x = 1
    time_old = datetime.now(timezone.utc)
    while not config.end_program:
   #     if     datetime.now(timezone.utc) > time_old + timedelta(0,5): #rate limiter, remove after testing

            #https://docs.astropy.org/en/latest/coordinates/example_gallery_plot_obs_planning.html
            telescope_time = datetime.now(timezone.utc) #have to make sure datetime is in utc for all the astro tools unless you specify it in them indivudally
            #end input
            
    
            #locate sky target
            try:
                target = SkyCoord.from_name(target_name) #sky objects (outside of solar system)
            except Exception as e:
                print('could not find target, trying body search')
                try:
                    target = get_body(target_name, Time(telescope_time), location = my_loc) #in solar system bodies
                except Exception as e:
                    print('could not find target or body')
                    print(e)
            #moon = get_body("moon", Time(telescope_time), location = my_loc)
            #sun = get_body("sun", Time(telescope_time), location = my_loc)
    
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


#Keep out zones (after getting desired commanded position)
def send_commands(my_locs, az, el, telescope_q):
    import config # https://docs.python.org/3/faq/programming.html#how-do-i-share-global-variables-across-modules

    import astropy.units as u
    from astropy.time import Time
    from datetime import datetime, date, timezone, timedelta
    from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body, get_sun
    from pygeomag import GeoMag

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
    
    #Check if will point into sun
    if (az-sun_az)**2 + (el - sun_el)**2 < keep_out_radius**2:
        send_command = False
        print("Keep out: pointing at sun")
    
    
    if send_command:
        #send_the_commands()
        telescope_q.put({'az':az, 'el':el, 'cmd':'move_az_el'})
    
    config.x = config.x + 1
    if config.x > 100:
        config.end_program = True
        

#tracking mode:
#def tracking:
    
    #realtime tracking of object in view at current moment
    #identify with selection or just pick thing in the center