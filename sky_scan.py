# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 10:28:56 2025

@author: AndrewMiller
"""
import astropy.units as u
from astropy.time import Time
from datetime import datetime, date, timezone, timedelta
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body, get_sun
from pygeomag import GeoMag
import multiprocessing as mp
import math
from threading import Thread
import queue
import sys
sys.path.append(r'C:\Users\AndrewMiller\OneDrive - Global Health Labs, Inc\Desktop')

#my functions
import config # https://docs.python.org/3/faq/programming.html#how-do-i-share-global-variables-across-modules
from telescope_control import telescope_control #actually drives the telescope
from camera_control import camera_control #takes photos
from image_save import image_save #saves photos
import tracking_modes #library of different pointing modes


config.end_program = False #reset globals for each run
config.az_angle = 0
config.el_angle = 0
config.camera_ready = False
config.x = 0 #use for quit debugging

#observation location details
#input location of observation
lat = 47.6205 #degrees
long = -122.3493 #degrees
altitude = 100 #meters
compass_dir = 0.000 #degrees, magnetic angle from 0 az on telescope to magnetic north
#compute magnetic declination 
geo_mag = GeoMag(coefficients_file="wmm/WMM_2025.COF")
telescope_time = datetime.now(timezone.utc) #have to make sure datetime is in utc for all the astro tools unless you specify it in them indivudally
mag_declination = geo_mag.calculate(glat=lat, glon=long, alt=altitude/1000, time=telescope_time.year+int(telescope_time.strftime('%j'))/1000) #altitude in km for geo_mag
#print('Magnetic declination: ' + str(mag_declination.d))
target_name = 'M33'
mode = 'astronomy'
camera_period = 10 #how many seconds between shots

my_locs = {'EarthLocation':EarthLocation(lat=lat * u.deg, lon = long * u.deg, height = altitude * u.m), 'compass_dir':compass_dir, 'mag_declination':mag_declination}

#end user inputs (this should live in GUI someday)

funct_dict = {
    "point_and_shoot": tracking_modes.point_and_shoot,
    "astronomy": tracking_modes.astronomy,
  #  "tracking": tracking_modes.tracking
}

camera_q = queue.Queue() #queue for taking a photo
image_save_q = queue.Queue()#  queue for saving photos
telescope_q = queue.Queue() #queue for telescope pointing
status_q = queue.Queue() #queue to return status items

queues = {'camera_q':camera_q, 'image_save_q':image_save_q, 'telescope_q':telescope_q, 'status_q':status_q}



# start the consumers
consumer_camera = Thread(target=camera_control, args=(queues,)) #You need to add a comma to create a tuple, so it would be args=(item,):
consumer_camera.start()

consumer_image_save = Thread(target=image_save, args=(queues,))
consumer_image_save.start()

consumer_telescope = Thread(target=telescope_control, args=(queues,))
consumer_telescope.start()

# start the producer
producer = Thread(target=funct_dict[mode], args=(my_locs, queues, target_name, camera_period)) #runs the mode in tracking_modes.py
producer.start()

print('All threads started')

# wait for all threads to finish
producer.join()
consumer_telescope.join()
consumer_image_save.join()
consumer_camera.join()

print('All threads finished')

# def point_and_shoot():
#     pass

# # def telescope_control(telescope_q):
# #     import time
# #     x = 0
# #     while(x<10):
# #         val = telescope_q.get()
# #         print('telescoepe control' + str(val))
# #         time.sleep(1) # Sleep for 3 seconds
# #         x = x + 1

# def camera_control(camera_q):
#     import time
#     x = 0
#     while(x<10):
#         print('camera control')
#         time.sleep(1) # Sleep for 3 seconds
#         x = x + 1
        
# def astronomy():
#     import time
#     x = 0
#     while(x<10):
#         print('astronomy')
#         telescope_q.put(x)
#         time.sleep(1) # Sleep for 3 seconds
#         x = x + 1
        


# #%%

# #point and shoot mode:
# #aim at [center or top left] and raster scan as described
# def point_and_shoot(FOV = 1000, az_steps = 10, el_steps = 5):
# # FOV = 1000 #microradians
# # az_steps = 10
# # el_steps = 5
#     for el in range (0,el_steps): #ADD: need to relate steps to FOV for decent coverage with overlaps for stitching
#         for az in range(0,az_steps):
#             done = False
#             while not done:
#                 in_postion = drive_gimbal(az,el)
#                 if in_position:
#                     take_picture() #filename as grid for stitching ease (0,0 to n,n)
#                     done = True
        

    
# #include way to scan by start/end lat/long coordinates
    
# #%%
    
# #astronomy mode:

# def astronomy(target_name = 'moon'):
#     x = 1
#     time = datetime.now(timezone.utc)
#     while x < 100:
#         if     datetime.now(timezone.utc) > time + timedelta(0,5):
            
#             #https://docs.astropy.org/en/latest/coordinates/example_gallery_plot_obs_planning.html
#             telescope_time = datetime.now(timezone.utc) #have to make sure datetime is in utc for all the astro tools unless you specify it in them indivudally
#             #end input
            
    
#             #locate sky target
#             my_loc = EarthLocation(lat=lat * u.deg, lon = long * u.deg, height = altitude * u.m)
#             try:
#                 target = SkyCoord.from_name(target_name) #sky objects (outside of solar system)
#             except Exception as e:
#                 print('could not find target, trying body search')
#                 try:
#                     get_body(target, Time(telescope_time), location = my_loc)
#                 except Exception as e:
#                     print('could not find target or body')
#             #moon = get_body("moon", Time(telescope_time), location = my_loc)
#             #sun = get_body("sun", Time(telescope_time), location = my_loc)
    
#             #need to convert moon to SkyCoord before proceeding to get alt/az
            
#             #altitude:  angle between the object and the observer's local horizon (0-90 deg)
#             #aziumuth: measured from true north and increasing eastward (N = 0, E = 90, S= 180, W = 270)
    
#             targetaltaz = target.transform_to(AltAz(obstime=telescope_time, location=my_loc))
#             #print(target_name + f"'s Altitude = {targetaltaz.alt:.2}")
#             #print(target_name + f"'s Azimuth = {targetaltaz.az:.2}")
            
#             #adjust azimuth for position of telescope base
#             az = targetaltaz.az.deg + compass_dir + mag_declination.d
#             el = targetaltaz.alt.deg
            
#             print('az: ' + str(az) + ' ' + 'el: ' + str(el))
#             x = x + 1
#     #print('el: ' + str(el))
# #%%
# #Keep out zones (after getting desired commanded position)
# def send_commands(az, el):
#     #Check if will point into ground
#     send_command = True
#     if el < 0:
#         send_command = False
#         print("Keep out: pointing into ground")
        
#     #get sun location
#     telescope_time = datetime.now(timezone.utc) #have to make sure datetime is in utc for all the astro tools unless you specify it in them indivudally
#     sun = get_body("sun", Time(telescope_time), location = my_loc)
#     sunaltaz = sun.transform_to(AltAz(obstime=telescope_time, location=my_loc))
#     sun_az = targetaltaz.az.deg + compass_dir + mag_declination.d
#     sun_el = targetaltaz.alt.deg
#     keep_out_radius = 5 #5 degree keep out radius to start
    
#     #Check if will point into sun
#     if (az-sun_az)^2 + (el - sun_el)^2 < keep_out_radius^2:
#         send_command = False
#         print("Keep out: pointing at sun")
    
    
#     if send_command:
#         #send_the_commands()
    
        
# #%%
# #tracking mode:
# def tracking:
#     #realtime tracking of object in view at current moment
#     #identify with selection or just pick thing in the center