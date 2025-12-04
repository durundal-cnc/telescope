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
import os
import csv
import statistics
import pandas as pd
import zipfile
import shutil
from collections import namedtuple

sys.path.append(r'/Users/andrewmiller/telescope/')

#my functions
import config # https://docs.python.org/3/faq/programming.html#how-do-i-share-global-variables-across-modules
from telescope_control import telescope_control #actually drives the telescope
from camera_control import camera_control #takes photos
from image_save import image_save #saves photos
import tracking_modes #library of different pointing modes

import sys
sys.path.append(r'/Users/andrewmiller/telescope/roboclaw_python')
from roboclaw_3 import Roboclaw
import math
import config # https://docs.python.org/3/faq/programming.html#how-do-i-share-global-variables-across-modules
from datetime import datetime, date, timezone, timedelta



#tool to align telescope images for stitching: https://astroalign.quatrope.org/en/latest/

def import_iphone_data(): #get location and orientation data from phone airdrop of CSV sensor logs
    #https://github.com/seemoo-lab/opendrop
    #and
    #https://www.tszheichoi.com/sensorlogger
    #get zip file via airdrop
    
    
    #unpack zip
    zip_file_loc = r'/Users/andrewmiller/telescope/14360_SE_Eastgate_Way-2025-07-18_17-03-55.zip'
    with zipfile.ZipFile(zip_file_loc, 'r') as zip_ref:
        unzip_dir = os.path.join(os.path.dirname(zip_file_loc),'iPhone_metadata')
        if os.path.exists(unzip_dir):
            shutil.rmtree(unzip_dir)
        os.mkdir(unzip_dir)
        zip_ref.extractall(unzip_dir)
    
    #import files
    sensor_files = r'/Users/andrewmiller/telescope/iPhone_metadata'
    files = [ name for name in os.listdir(sensor_files) if name[-4:] == '.csv' ]
    return_dict = {}
    for file in files:
        print('loading ' + file)
        if file == 'Compass.csv':
            #get stats on recorded values
            df = pd.read_csv(os.path.join(sensor_files,file))
            return_dict['compass'] = drop_outliers(df['magneticBearing'])
        if file == 'Location.csv':
            df = pd.read_csv(os.path.join(sensor_files,file))
            return_dict['altitude'] = drop_outliers(df['altitude']) #also altitudeAboveMeanSeaLevel
            return_dict['latitude'] = drop_outliers(df['latitude'])
            return_dict['longitude'] = drop_outliers(df['longitude'])
    return return_dict

def drop_outliers(input_list,num_stddevs=1):
    stddev = statistics.stdev(input_list)
    mean = statistics.mean(input_list)
    less_outliers = [item for item in input_list if item <= mean+num_stddevs*stddev or item >= mean-num_stddevs*stddev]
    avg_less_outliers = statistics.mean(less_outliers)
    print('dropped ' + str(len(input_list)-len(less_outliers)))
    return avg_less_outliers


#mount calibration
#goto given star, pick it out of an image, get error from predicted / actual based on deg/pixel in the image, repeat?
#or implement feature tracking in images

#make list of in-view satellites from location and list in/out view times for them as pre-processor to select targets
#use class for targets with associated data?

config.end_program = False #reset globals for each run
config.az_angle = 0
config.el_angle = 0
config.az_in_position = False
config.el_in_potiion = False
config.camera_ready = False
config.tracking_ready = False
config.camera_process_ready = False
config.image_process_ready = False
config.telescope_process_ready = False
config.coordinates = []



config.x = 0 #use for quit debugging


#observation location details
#input location of observation
#lat = 47.6205 #degrees
#long = -122.3493 #degrees
#altitude = 100 #meters, I think this is all WGS ellipsoid referenced but needs checking
#compass_dir = 0.000 #degrees, magnetic angle from 0 az on telescope to magnetic north

iphone = import_iphone_data()
lat = iphone['latitude']
long = iphone['longitude']
altitude = iphone['altitude']
compass_dir = iphone['compass']
config.lat = lat
config.long = long

#compute magnetic declination 
geo_mag = GeoMag(coefficients_file="wmm/WMM_2025.COF")
telescope_time = datetime.now(timezone.utc) #have to make sure datetime is in utc for all the astro tools unless you specify it in them indivudally
mag_declination = geo_mag.calculate(glat=lat, glon=long, alt=altitude/1000, time=telescope_time.year+int(telescope_time.strftime('%j'))/1000) #altitude in km for geo_mag
#print('Magnetic declination: ' + str(mag_declination.d))
target_name = 'ISS' #'Starlink-4727' #'M33' 'sun' 'moon 'SDO' 'ISS'
mode = 'satellite_tracking' #point and shoot, satellite tracking, astronomy
camera_period = 10 #how many seconds between shots
target_time_start = datetime.now()
target_time_end = datetime.now() + timedelta(hours=1)
timespacing = timedelta(seconds=1)

my_locs = {'EarthLocation':EarthLocation(lat=lat * u.deg, lon = long * u.deg, height = altitude * u.m), 'compass_dir':compass_dir, 'mag_declination':mag_declination}
config.my_locs = my_locs
config.mag_declination = mag_declination
config.compass_dir = compass_dir

#end user inputs (this should live in GUI someday)

target_track = namedtuple('target_track', ['target_name', 'time', 'az', 'el']) #think through how this should be organized - basic lists, dicts, other?
# from collections import namedtuple
# Point = namedtuple('Point', 'x y')
# pt1 = Point(1.0, 5.0)
# pt2 = Point(2.5, 1.5)


#for each target compute the time and coordinates:
    #first add the sun to be able to check for keep outs:
target_list = [['sun','astronomy'], ['moon','astronomy'], ['polaris','astronomy'], ['M42', 'astronomy'], ['ISS', 'satellite_tracking'] ]
#target_list = [['sun','astronomy'], ['moon','astronomy'], ['polaris','astronomy'], ['M42', 'astronomy'] ]

tracks = []

for targetpair in target_list:
    print('searching for target ' + str(targetpair))
    target = targetpair[0]
    mode = targetpair[1]
    match mode:
        case 'satellite_tracking':
            # Code to execute if subject matches pattern_1
            x = satellite_tracking(my_locs, target_name = 'ISS', timespacing = timedelta(seconds=1), timestart = datetime.now(timezone.utc))

#need to check that time end works properly when specified           #x = satellite_tracking(my_locs, target_name = 'ISS', timespacing = timedelta(seconds=1), timestart = datetime.now(timezone.utc), timeend = datetime.now(timezone.utc) + timedelta(seconds=60))
            tracks = tracks +[target_track(*i) for i in x]#(*i) unpacks the list as an input to the namedtuple

            print(mode)
        case 'astronomy':
            # Code to execute if subject matches pattern_2
            #x = astronomy(my_locs, target_name = target, timespacing = timedelta(seconds=1), timestart = datetime.now(timezone.utc), timeend = datetime.now(timezone.utc) + timedelta(seconds=60))

            x = astronomy(my_locs, target_name = target, timespacing = timedelta(seconds=1), timestart = datetime.now(timezone.utc))
            tracks = tracks +[target_track(*i) for i in x]#(*i) unpacks the list as an input to the namedtuple
    
            print(mode)
        case 'point_and_shoot':
            #
            print(mode)
        case _:
            print(mode)
            
#filter the output (optional)
specific_track = sorted(tracks, key=lambda x: x.target_name)
x = [record for record in tracks if record.target_name is 'sun']


#quick printout for checking
for x in tracks:
    print(  x.time.strftime("%Y-%m-%d %H:%M:%S")  + ' ' + str(x.az) + ' ' + str(x.el)  + ' ' + x.target_name)


#check for sun angle inclusion
#get time bounds of min-max of all times
all_times = [record.time for record in tracks]

observations_start = min(all_times)
observations_end   = max(all_times)

#compute sun location for all observations
sun = astronomy(my_locs, target_name = 'Sun', timespacing = timedelta(seconds=1), timestart = observations_start, timeend = observations_end)
sun =[target_track(*i) for i in sun]#(*i) unpacks the list as an input to the namedtuple

keepout = 5 #keep out degrees from the sun
sunintrusions = []
for track,sunloc in zip(tracks, sun):
    if abs(track.az-sunloc.az) < keepout and abs(track.el-sunloc.el) < keepout:
        sunintrusions.append([track,sunloc])
        print('intrusion into sun: ')    
        print(track)
            
print('Found ' + str(len(sunintrusions)) + ' sun keep out intrusions')
intrusions = set([x[0].target_name for x in sunintrusions])
print('Targets with sun intrustions: ' + str(intrusions))

#add stop point or drop any sun intrusion runs

#order runs by time of observation
observations = set([x.target_name for x in tracks])
individual_tracks = []
time_of_individual_tracks = []
for observation in observations:
    x = [record for record in tracks if record.target_name is observation]
    individual_tracks.append(x)
    time_of_individual_tracks.append(x[0].time)
time_of_individual_tracks, individual_tracks = zip(*sorted(zip(time_of_individual_tracks, individual_tracks )))
name_of_individual_tracks = [x[0].target_name for x in individual_tracks]
individual_tracks = dict(zip(name_of_individual_tracks, individual_tracks))

# funct_dict = {
#     "point_and_shoot": tracking_modes.point_and_shoot,
#     "astronomy": tracking_modes.astronomy,
#     "satellite_tracking": tracking_modes.satellite_tracking
# }

camera_q = queue.Queue() #queue for taking a photo
image_save_q = queue.Queue()#  queue for saving photos
# telescope_q = queue.Queue() #queue for telescope pointing
status_q = queue.Queue() #queue to return status items

queues = {'camera_q':camera_q, 'image_save_q':image_save_q, 'telescope_q':telescope_q, 'status_q':status_q}

#figure out max command rate roboclaw can receive, build table of locations to spit out at that rate (vice real-time computations)


# start the consumers
consumer_camera = Thread(target=camera_control, args=(queues,)) #You need to add a comma to create a tuple, so it would be args=(item,):
consumer_camera.start()
while not config.camera_process_ready:
    pass
print('consumer_camera started')


consumer_image_save = Thread(target=image_save, args=(queues,))
consumer_image_save.start()
while not config.image_process_ready:
    pass
print('consumer_image_save started')

# start the producer
# thread_tracker = Thread(target=funct_dict[mode], args=(my_locs, queues, target_name, camera_period)) #runs the mode in tracking_modes.py
# thread_tracker.start()

# while not config.tracking_ready:
#     pass
# print('thread_tracker started')

# consumer_telescope = Thread(target=telescope_control, args=(queues,))
# consumer_telescope.start()
# while not config.telescope_process_ready:
#     pass
# print('consumer_telescope started')

print('All threads started')
statuses = []
time = Time.now()

stats, rc = telescope_control(cmd='init') #initialize telescope controller

########## start the telescope control loop
x = 1
seq_num = 0
coordinate_loc = 0
state = 'idle'
last_print = datetime.now(timezone.utc)
time_print_spacing = 1 #seconds between prints for waiting for a target to come in view

while not config.end_program:
    #process user commands
    keyboard = ''

    if keyboard == 'q': #Time.now() > time + timedelta(seconds = 240):
        config.end_program = True
        print('ending program due to keyboard input')
    elif keyboard == 'new_targets':
        print('Input new target list ')
        #make the targeting a funtion where you input targets and times here
    elif keyboard == 'az':
        print(config.az_angle)
    elif keyboard == 'el':
        print(config.el_angle)   
    elif keyboard == 'new_target':
        print('List of available targets and times: ')
        '{0: <5}'.format('s')
        [print('{0: <20}'.format(str(x)) + ' ' +  str(y)) for x,y in zip(name_of_individual_tracks, time_of_individual_tracks)]
        target = input('input new target')
        active_track = individual_tracks[target]
        state == 'track'
        coordinate_loc = 0
    elif keyboard == 'run_seq':
        active_track = individual_tracks[seq_num]
        auto_next_track = True

    elif keyboard == 'manual':
        user_az = input('Input az angle')
        user_el = input('Input el angle')
        state = 'manual'
    elif keyboard == 'idle':
        state = 'idle'
    else:
        print('status check')
        #check for status every n seconds    
    
    #perform commands
    if state == 'idle':
        continue
    elif state == 'manual':
        #send command for manual move
        stats, rc = telescope_control(rc, cmd='move_az_el', coord=[datetime.now(timezone.utc), user_az, user_el]) #initialize telescope controller

    elif state == 'track':    #send the commands to the telescope for tracking
        if datetime.now(timezone.utc) < active_track[coordinate_loc].time: #this will keep the next coordinates from being sent from the caller
            if datetime.now(timezone.utc) > last_print:
                last_print = datetime.now(timezone.utc) + timedelta(seconds=time_print_spacing)
                print('Waiting for new coordinate in ' + str(active_track[coordinate_loc].time - datetime.now(timezone.utc) )+' H:M:S')
            continue #wait for coordinate time to occur
        else:
            if coordinate_loc > len(active_track): #if we have completed a track
                print('Completed track')
                if auto_next_track is True: #if we aren't in manual track selection
                    #move to next track
                    seq_num = seq_num + 1 #move to the next one
                    if seq_num <= len(individual_tracks): #make sure we have a next one to move to
                        print('selecting next track:')
                        print('{0: <20}'.format(str(name_of_individual_tracks[seq_num])) + ' ' +  str(time_of_individual_tracks[seq_num])) 
                        active_track = individual_tracks[seq_num]
                    else:
                        print('Ran out of tracks, moving to idle state')
                        state = 'idle'
                        continue
                else: #proceed to send coords
                    continue
            coords = active_track.iloc[coordinate_loc]
            stats, rc = telescope_control(rc, cmd='move_az_el', coords = [datetime.now(timezone.utc), 0,0]) #send next coordinate set to telescope
            coordinate_loc = coordinate_loc+1
    elif state == 'home':
        print('Sending homing sequence')
        #send the home command
    else:
        print('unknown state')
        
##### queue managing/shutdown below

    size = queues['status_q'].qsize() #check that an image is in the queue
    if size >= 1:
        #get status_q and display message
        status = status_q.get()
        statuses.append(status)
        #print(status)
    #sleep(1)
    #keyboard = input()
    if Time.now() > time + timedelta(seconds = 60):
        config.end_program = True
        print('ending program due to timer')

        
    
#    else:
#        print('End conditions: ' + str(size) + ' ' + str((time + timedelta(seconds = 20)-Time.now()).sec) + ' s remaining')

# wait for all threads to finish
thread_tracker.join()
print('thread_tracker join() complete')
consumer_telescope.join()
print('consumer_telescope join() complete')
# consumer_image_save.join()
# print('consumer_image_save join() complete')
# consumer_camera.join()
# print('consumer_camera join() complete')


print('All threads finished')



#y = [j - i for i,j in zip(config.x[:-1], config.x[1:])]

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