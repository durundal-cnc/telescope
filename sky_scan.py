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

#GUI tools, using Beeware
from time import time
from kivy.app import App
from os.path import dirname, join
from kivy.lang import Builder
from kivy.properties import (
    NumericProperty,
    StringProperty,
    BooleanProperty,
    ListProperty,
)
from kivy.clock import Clock
from kivy.animation import Animation
from kivy.uix.screenmanager import Screen
#end GUI tools







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

def initialize_config():
    config.end_program = False #reset globals for each run
    config.az_angle_SV = 0
    config.el_angle_SV = 0
    config.az_angle_PV = 0
    config.el_angle_PV = 0
    config.az_in_position = False
    config.el_in_potiion = False
    config.camera_ready = False
    config.camera_period = 10 #how many seconds between shots
    config.tracking_ready = False
    config.camera_process_ready = False
    config.image_process_ready = False
    config.telescope_process_ready = False
    config.coordinates = []
    config.roboclaw_stats = dict()
    
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

    my_locs = {'EarthLocation':EarthLocation(lat=lat * u.deg, lon = long * u.deg, height = altitude * u.m), 'compass_dir':compass_dir, 'mag_declination':mag_declination}
    config.my_locs = my_locs
    config.mag_declination = mag_declination
    config.compass_dir = compass_dir
    
    #end user inputs (this should live in GUI someday)



######## start the GUI



######## GUI running


######## Get the target coordinates
def compute_target_coords(target_list = [['sun','astronomy'], ['moon','astronomy'], ['polaris','astronomy'], ['M42', 'astronomy'], ['ISS', 'satellite_tracking'] ]):
    #for each target compute the time and coordinates:
        #first add the sun to be able to check for keep outs:
    #target_list = [['sun','astronomy'], ['moon','astronomy'], ['polaris','astronomy'], ['M42', 'astronomy'], ['ISS', 'satellite_tracking'] ]
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
                
    #filter the output (optional examples)
    #specific_track = sorted(tracks, key=lambda x: x.target_name)
    #x = [record for record in tracks if record.target_name is 'sun']
    
    
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
    
    keep_out_radius = 5 #keep out degrees from the sun
    sunintrusions = []
    for track,sunloc in zip(tracks, sun):
#        if abs(track.az-sunloc.az) < keepout and abs(track.el-sunloc.el) < keepout:
        if (track.az-sunloc.az)**2 + (track.el - sunloc.ell)**2 < keep_out_radius**2:
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
    return time_of_individual_tracks, name_of_individual_tracks, individual_tracks


target_track = namedtuple('target_track', ['target_name', 'time', 'az', 'el']) #think through how this should be organized - basic lists, dicts, other?
initialize_config()
# from collections import namedtuple
# Point = namedtuple('Point', 'x y')
# pt1 = Point(1.0, 5.0)
# pt2 = Point(2.5, 1.5)


camera_q = queue.Queue() #queue for taking a photo
image_save_q = queue.Queue()#  queue for saving photos
status_q = queue.Queue() #queue to return status items

queues = {'camera_q':camera_q, 'image_save_q':image_save_q, 'status_q':status_q}

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
print('All threads started')



######## start the hardware
statuses = []
stats, rc = telescope_control(cmd='init') #initialize telescope controller

########## start the telescope control loop
x = 1
seq_num = 0
coordinate_loc = 0
state = 'idle'
last_print = datetime.now(timezone.utc)
time_print_spacing = 1 #seconds between prints for waiting for a target to come in view
time_of_individual_tracks, name_of_individual_tracks, individual_tracks = compute_target_coords(target_list = [['sun','astronomy'], ['moon','astronomy'], ['polaris','astronomy'], ['M42', 'astronomy'], ['ISS', 'satellite_tracking'] ])

while not config.end_program:
    #process user commands
    keyboard = ''

#things to add for GUI
#manual +/- az/el control buttons and input numbers
#point and shoot select start point and end point
#display for images
#autofocus button
#some sort of progress bar
#debug function to display variables of interest (especially config variables)


#make separate functions bound to each button for each of these
    if keyboard == 'q': #Time.now() > time + timedelta(seconds = 240):
        config.end_program = True
        print('ending program due to keyboard input')
    elif keyboard == 'new_targets':
        print('Input new target list ')
        time_of_individual_tracks, name_of_individual_tracks, individual_tracks = compute_target_coords(target_list = [['sun','astronomy'], ['moon','astronomy'], ['polaris','astronomy'], ['M42', 'astronomy'], ['ISS', 'satellite_tracking'] ])

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


#make this a bound always-executing function
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

    size = queues['status_q'].qsize() #check that an image is in the queue (collects the camera and image save returns from their separate threads)
    if size >= 1:
        #get status_q and display message
        status = status_q.get()
        statuses.append(status)
        #print(status)
    #sleep(1)
    #keyboard = input()
    # if Time.now() > time + timedelta(seconds = 60):
    #     config.end_program = True
    #     print('ending program due to timer')

# wait for all threads to finish
# thread_tracker.join()
# print('thread_tracker join() complete')
# consumer_telescope.join()
# print('consumer_telescope join() complete')
consumer_image_save.join()
print('consumer_image_save join() complete')
consumer_camera.join()
print('consumer_camera join() complete')

print('All threads finished')