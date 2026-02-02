# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 10:28:56 2025

@author: AndrewMiller
"""
import astropy.units as u
from astropy.time import Time
from time import sleep
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
from tracking_modes import astronomy, satellite_tracking, point_and_shoot

import sys
sys.path.append(r'/Users/andrewmiller/telescope/roboclaw_python')
from roboclaw_3 import Roboclaw
import math
import config # https://docs.python.org/3/faq/programming.html#how-do-i-share-global-variables-across-modules
from datetime import datetime, date, timezone, timedelta

#GUI tools, using Beeware
from time import time
from kivy.app import App
from kivy.core.window import Window

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

from kivy.app import App, async_runTouchApp
from kivy.uix.gridlayout import GridLayout
from kivy.uix.label import Label
from kivy.uix.textinput import TextInput
from kivy.uix.button import Button
from kivy.uix.image import AsyncImage
from kivy.uix.dropdown import DropDown
from kivy.uix.checkbox import CheckBox
from kivy.metrics import dp
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.progressbar import ProgressBar
from kivy.uix.checkbox import CheckBox
from kivy.uix.scrollview import ScrollView
from kivy_garden.matplotlib.backend_kivyagg import FigureCanvasKivyAgg #python -m pip install https://github.com/kivy-garden/matplotlib/archive/master.zip
from kivy.clock import Clock


from functools import partial
import trio #for asynch running of the telescope thread



from kivy.properties import (
    NumericProperty,
    StringProperty,
    BooleanProperty,
    ListProperty,
)
import re
import sys

keep_running_app = True

#end GUI tools



#TO ADD
#DONE fix elevation behavior around 360 degrees (e.g. handle0 negative values)
#DONE to determine if lookup mode should be astronomy or satellite_tracking check to see if target is in NORAD csv, if so, use sat track, otherwise astronomy
#DONE after getting close to target read the velocity when target is reached and instead of letting the position command finish issue a velocity command while waiting for the next coordinate?
#DONE figure out why a second track isn't working right (wrong coordinates?) without rebooting program (coordinate times are in the future - only getting set on boot?)
#DONE figure out how to make program close gracefully
#Fine tune the pointing during tracking so it can converge on teh right value with the lookahead keeping velocity up (or just get close enough and proceed at constant velocity if the target will be in FOV?)
#Elevation SV value below 0 goes negative (no wraparound code), causes issue when halt telescope called (big swing to those coordinates)
#DONE blank the screen for compute_targets so obvisou when done calculating
#DONE add sun inclusion checking to the actual track checker
#Decide what to have the sun checker do if sun is detected (don't execute? Go up to sun angle and stop? Warn with a popup or something but proceed?)
#implement offline storage of all TLE for offline use if no internet (Celestrak has downloads?)
#DONE - this is a bug in the fake roboclaw logic somewhere - manual negative command not wrapping around properly now? make sure not a bug in the non-Roboclaw connected code
#investigate if driving to the first points at full speed until you get to near zero pointing error (it will be pulsey as it closes in) then engaging controlled speed as a binary works
#install limit switches and get home function working
#implement FITS format for photos (maybe? useful for stitching w/ metadata embedded?)
#handle can't-find-object case gracefully (crashes out now)
#idx, d2d, d3d = c.match_to_catalog_3d(catalog) #see if this can spit back the matched object in skycoordinates (since input may be partial match)
#function to add: for location in sky, find objects within https://eteq-astropy.readthedocs.io/en/latest/coordinates/matchsep.html
#adjust elevation from 0-180 to -90 to +90 for skycoordinate compatibility (? unsure how this is supposed to work)

#tool to align telescope images for stitching: https://astroalign.quatrope.org/en/latest/

def get_ra_dec(): #convert current alt-az angle to right ascension and declination
    current_time = Time(datetime.now(timezone.utc)) #have to make sure datetime is in utc for all the astro tools unless you specify it in them indivudally
    loc = config.my_locs['EarthLocation']

    sc = SkyCoord(alt=config.el_angle_SV*u.deg, az=config.az_angle_SV*u.deg, obstime = current_time, frame = 'altaz', location = loc)
   # sc = SkyCoord(alt=config.az_angle_PV*u.deg, az=config.el_angle_PV*u.deg, frame=altaz_frame) #need to put in the iPhone compass direction offsets
    icrs_coo = sc.transform_to('icrs')
   #print(f"RA/Dec: {icrs_coo.ra.deg} ({icrs_coo.ra.to_string(unit=u.hourangle, sep=':')}), {icrs_coo.dec.deg}")

    return icrs_coo.ra.deg, icrs_coo.dec.deg

def plot_nearby_stars(minimum_brightness_plot = 13, minimum_brightness_annotation = 6, radius = 10):
    from astropy.table import QTable
    from astroquery.gaia import Gaia
    Gaia.ROW_LIMIT = 10000  # Set the row limit for returned data
    
    import numpy as np
    from astropy.coordinates import SkyCoord, match_coordinates_sky
    from astropy import units as u
    
    from astroquery.simbad import Simbad
    
    
   #  #take earth location
   #  loc = config.my_locs['EarthLocation']
   #  current_time = Time(datetime.now(timezone.utc)) #have to make sure datetime is in utc for all the astro tools unless you specify it in them indivudally
   #  #take current az/el at that location
   #  altaz_frame = AltAz(obstime=current_time, location=loc)
   #  #error: thinks SkyCoord alt az should be -90 to +90 latitude??
   #  print('alt ' + str(config.el_angle_PV) + ' az ' + str(config.az_angle_PV))
    

   #  sc = SkyCoord(alt=config.el_angle_PV*u.deg, az=config.az_angle_PV*u.deg, obstime = current_time, frame = 'altaz', location = loc)
   # # sc = SkyCoord(alt=config.az_angle_PV*u.deg, az=config.el_angle_PV*u.deg, frame=altaz_frame) #need to put in the iPhone compass direction offsets
   #  icrs_coo = sc.transform_to('icrs')
   #  print(f"RA/Dec: {icrs_coo.ra.deg}, {icrs_coo.dec.deg}")
    
   #  #import GAIA database of objects
   #  #TODO: figure out how to store star catalog locally
   #  #get list of objects within radius of that skycoord 
   #  job = Gaia.cone_search_async(sc, radius=radius * u.deg)
   #  ngc188_table = job.get_results()

   #  # only keep stars brighter than G=19 magnitude
   #  ngc188_table = ngc188_table[ngc188_table["phot_g_mean_mag"] < minimum_brightness_plot * u.mag] #13 is limit of ~4" telescope
    
   #  # add all ids in the SIMBAD results
   #  Simbad.add_votable_fields('ids')
   #  named = Simbad.query_objects(ngc188_table['designation'])
    
    
   #  #for readability
   #  import pandas as pd
   #  p = ngc188_table.to_pandas()      #convert to Pandas dataframe
   #  q = named.to_pandas()
   #  #plot according to brightness
   #  #gaia_dist = Distance(parallax=ngc188_table_3d["parallax"].filled(np.nan)) #unused
   #  gaia_magnitude = ngc188_table["phot_g_mean_mag"].filled(np.nan)
    
    
    
    #plot and label so can match what seeing through eyepiece to what's on sky there
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.5, 5.2), constrained_layout=True)
    # cs = ax.scatter(
    #     ngc188_table['ra'],
    #     ngc188_table['dec'],
    #     c=gaia_magnitude,
    #     s=5, #marker size
    #     vmin=min(gaia_magnitude),
    #     vmax=max(gaia_magnitude),
    #     cmap="gray",
    # )
    # cb = fig.colorbar(cs)
    # cb.set_label(f"magnitude")
    
    # #restrict name plotting to brighter stars
    # for i, txt in enumerate(named['main_id']):
    #     if ngc188_table['phot_g_mean_mag'][i] < minimum_brightness_annotation:
    #         ax.annotate(txt, (ngc188_table['ra'][i], ngc188_table['dec'][i]), fontsize=8)
    
    ax.set_xlabel("RA [deg]")
    ax.set_ylabel("Dec [deg]")
    
    ax.set_title("Gaia DR2 sources near NGC 188", fontsize=18)
    return plt.gcf() 


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

def initialize_config():
    config.end_program = False #reset globals for each run

    config.camera_period = 10 #how many seconds between shots
    config.tracking_ready = False
    config.camera_process_ready = False
    config.image_process_ready = False
    config.telescope_process_ready = False
    config.coordinates = []
    config.current_coord_loc = '' #the current position in the target_track list
    config.selected_target = '' #name of the selected target
    config.selected_target_coords = [] #the list of target_track namedtuple target coords
    config.roboclaw_stats = dict()
    config.state = 'idle'
    config.selected_target_dropdown = 0
    
    config.az_angle_SV = 0 #degrees Azimuth is oriented East of North (i.e., N=0, E=90 degrees)
    config.el_angle_SV = 0 #degrees Altitude/elevation: Zero and 180º is your horizon, and 90º is straight above (zenith)
    config.az_angle_PV = 0 #degrees
    config.el_angle_PV = 0 #degrees
    config.az_in_position = False
    config.el_in_position = False
    config.camera_ready = False
    config.az_encoder_value = 0 #counts
    config.el_encoder_value = 0 #counts
    config.az_pointing_error = 0 #degrees
    config.el_pointing_error = 0 #degrees
    config.total_pointing_error = 0 #degrees
    config.Bit1_az = 0 #backward/forward indicator
    config.Bit1_el = 0 #backward/forward indicator
    config.az_speed = 0 #counts/sec
    config.el_speed = 0 #counts/sec
    config.az_current = 0 #current draw in 10 mA increments
    config.el_current = 0
    
    
    config.point_and_shoot_start = 0
    config.point_and_shoot_end = 0
    config.point_and_shoot_FOV = 1 #field of view in degreees
    

    config.x = 0 #use for quit debugging
    config.log = '' #store text for looking at performance
    
    
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
    config.mag_declination = 0# mag_declination #zero for debugging until can actually put phone on stage
    config.compass_dir = 0 #compass_dir #zero for debugging until can actually put phone on stage
    
    #end user inputs (this should live in GUI someday)



######## start the GUI



######## GUI running

######## Sun-keepout check function
def check_sun_keepout(tracks,observations_start,observations_end):
    #compute sun location for all observations
    sun = astronomy(config.my_locs, target_name = 'Sun', timespacing = timedelta(seconds=1), timestart = observations_start, timeend = observations_end)
    sun =[target_track(*i) for i in sun]#(*i) unpacks the list as an input to the namedtuple
    
    keep_out_radius = 5 #keep out degrees from the sun
    sunintrusions = []
    for track,sunloc in zip(tracks, sun):
#        if abs(track.az-sunloc.az) < keepout and abs(track.el-sunloc.el) < keepout:
        if (track.az-sunloc.az)**2 + (track.el - sunloc.el)**2 < keep_out_radius**2:
            sunintrusions.append([track,sunloc])
            print('intrusion into sun: ')    
            print(track)
                
    print('Found ' + str(len(sunintrusions)) + ' sun keep out intrusions')
    intrusions = set([x[0].target_name for x in sunintrusions])
    print('Targets with sun intrustions: ' + str(intrusions))
    return intrusions
    

######## Get the target coordinates
def compute_target_coords(target_list = [['sun','astronomy'], ['moon','astronomy'], ['polaris','astronomy'], ['M42', 'astronomy'], ['ISS', 'satellite_tracking'] ], coincident_times = False, duration=10):
   
    #for each target compute the time and coordinates:
        #first add the sun to be able to check for keep outs:
    #target_list = [['sun','astronomy'], ['moon','astronomy'], ['polaris','astronomy'], ['M42', 'astronomy'], ['ISS', 'satellite_tracking'] ]
    #target_list = [['sun','astronomy'], ['moon','astronomy'], ['polaris','astronomy'], ['M42', 'astronomy'] ]
    
    tracks = []
    timestart = datetime.now(timezone.utc)
    
    for targetpair in target_list:
        print('searching for target ' + str(targetpair))
        target = targetpair[0]
        mode = targetpair[1]
        match mode:
            case 'satellite_tracking':
                # Code to execute if subject matches pattern_1
                x = satellite_tracking(config.my_locs, target_name = target, timespacing = timedelta(seconds=1), timestart = timestart, timeend = timestart + timedelta(seconds=duration))
    
    #need to check that time end works properly when specified           #x = satellite_tracking(my_locs, target_name = 'ISS', timespacing = timedelta(seconds=1), timestart = datetime.now(timezone.utc), timeend = datetime.now(timezone.utc) + timedelta(seconds=60))
                tracks = tracks +[target_track(*i) for i in x]#(*i) unpacks the list as an input to the namedtuple
    
                print(mode)
            case 'astronomy':
                # Code to execute if subject matches pattern_2
                #x = astronomy(my_locs, target_name = target, timespacing = timedelta(seconds=1), timestart = datetime.now(timezone.utc), timeend = datetime.now(timezone.utc) + timedelta(seconds=60))
    
                x = astronomy(config.my_locs, target_name = target, timespacing = timedelta(seconds=1), timestart = timestart, timeend = timestart + timedelta(seconds=duration))
                tracks = tracks +[target_track(*i) for i in x]#(*i) unpacks the list as an input to the namedtuple
        
                print(mode)
            case 'point_and_shoot':
                x = point_and_shoot(config.my_locs, FOV = config.point_and_shoot_FOV, slew_speed = 5) #units degrees, degrees/sec
                tracks = tracks +[target_track(*i) for i in x]#(*i) unpacks the list as an input to the namedtuple

                print(mode)
            case _:
                print(mode)
        if coincident_times == False:
            timestart = tracks[-1].time #reset the start time to the last track's final time (otherwise it gets tracks for all objects coincidentally in time)
                
    #filter the output (optional examples)
    #specific_track = sorted(tracks, key=lambda x: x.target_name)
    #x = [record for record in tracks if record.target_name is 'sun']
    
    
    #quick printout for checking
    # for x in tracks:
    #     print(  x.time.strftime("%Y-%m-%d %H:%M:%S")  + ' ' + str(x.az) + ' ' + str(x.el)  + ' ' + x.target_name)
    
    
    #check for sun angle inclusion
    #get time bounds of min-max of all times
    all_times = [record.time for record in tracks]
    
    observations_start = min(all_times)
    observations_end   = max(all_times)
    
    
    #moved to helper function
    intrusions = check_sun_keepout(tracks,observations_start,observations_end)

#     #compute sun location for all observations
#     sun = astronomy(config.my_locs, target_name = 'Sun', timespacing = timedelta(seconds=1), timestart = observations_start, timeend = observations_end)
#     sun =[target_track(*i) for i in sun]#(*i) unpacks the list as an input to the namedtuple
#    
#     keep_out_radius = 5 #keep out degrees from the sun
#     sunintrusions = []
#     for track,sunloc in zip(tracks, sun):
# #        if abs(track.az-sunloc.az) < keepout and abs(track.el-sunloc.el) < keepout:
#         if (track.az-sunloc.az)**2 + (track.el - sunloc.el)**2 < keep_out_radius**2:
#             sunintrusions.append([track,sunloc])
#             print('intrusion into sun: ')    
#             print(track)
#                
#     print('Found ' + str(len(sunintrusions)) + ' sun keep out intrusions')
#     intrusions = set([x[0].target_name for x in sunintrusions])
#     print('Targets with sun intrustions: ' + str(intrusions))
    
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
    
    def compute_velocities(track):  #get the velocities between all the points
        if len(track)>1:
            first_coord = track[0]
            for idx, coord in enumerate(track[1:]):
                print(coord)
                delta_time = (coord.time-first_coord.time).total_seconds() #time between coords in seconds
                track[idx] = track[idx]._replace(az_vel = abs((coord.az  - first_coord.az)/delta_time)) #deg/sec
                track[idx] = track[idx]._replace(el_vel = abs((coord.el  - first_coord.el)/delta_time)) #deg/sec
            track[0] = track[0]._replace(az_vel = track[1].az_vel) #assume this probably won't matter at the start of the track as it will need time to get to location and speed anyway
            track[0] = track[0]._replace(el_vel = track[1].el_vel) #assume this probably won't matter at the start of the track as it will need time to get to location and speed anyway
            return track
        else:
            print('Error: track does not have any elements in it')
            return
    
    for key, value in individual_tracks.items(): #individual_tracks will have dict item for each track
        print(key)
        print(value)
        individual_tracks[key] = compute_velocities(value) #get the velocities between all the points
    
    return time_of_individual_tracks, name_of_individual_tracks, individual_tracks
#time_of_individual_tracks, name_of_individual_tracks, individual_tracks =  compute_target_coords(target_list = [['sun','astronomy'], ['moon','astronomy'], ['polaris','astronomy'], ['M42', 'astronomy'], ['ISS', 'satellite_tracking'] ])

target_track = namedtuple('target_track', ['target_name', 'time', 'az', 'el', 'az_vel', 'el_vel']) #think through how this should be organized - basic lists, dicts, other?
stats = namedtuple('stats', ['az_PV', 'el_PV', 'error_conds', ]) #think through how this should be organized - basic lists, dicts, other?

initialize_config() #set up location on earth and global variables
# from collections import namedtuple
# Point = namedtuple('Point', 'x y')
# pt1 = Point(1.0, 5.0)
# pt2 = Point(2.5, 1.5)

# camera_q = queue.Queue() #queue for taking a photo
# image_save_q = queue.Queue()#  queue for saving photos
# status_q = queue.Queue() #queue to return status items

# queues = {'camera_q':camera_q, 'image_save_q':image_save_q, 'status_q':status_q}

# #figure out max command rate roboclaw can receive, build table of locations to spit out at that rate (vice real-time computations)

# # start the consumers
# consumer_camera = Thread(target=camera_control, args=(queues,)) #You need to add a comma to create a tuple, so it would be args=(item,):
# consumer_camera.start()
# while not config.camera_process_ready:
#     pass
# print('consumer_camera started')


# consumer_image_save = Thread(target=image_save, args=(queues,))
# consumer_image_save.start()
# while not config.image_process_ready:
#     pass
# print('consumer_image_save started')
# print('All threads started')




########## start the telescope control loop

#time_of_individual_tracks, name_of_individual_tracks, individual_tracks = compute_target_coords(target_list = [['sun','astronomy'], ['moon','astronomy'], ['polaris','astronomy'], ['M42', 'astronomy'], ['ISS', 'satellite_tracking'] ])


async def run_app_GUI(root, nursery):
    '''This method, which runs Kivy, is run by trio as one of the coroutines.
    '''
    # trio needs to be set so that it'll be used for the event loop
    try:
        await async_runTouchApp(root, async_lib='trio')  # run Kivy, asynch_runTouchApp is kivy function
    finally:
        print('App done')
        # now cancel all the other tasks that may be running
        nursery.cancel_scope.cancel()



async def run_telescope(root):
    '''This method is also run by trio and periodically does something.'''
    last_print = datetime.now(timezone.utc)
    x = 1
    seq_num = 0
    coordinate_loc = 0
    time_print_spacing = 1 #seconds between prints for waiting for a target to come in view
    do_once = True

    ######## start the hardware
    stats, rc, address = telescope_control(cmd='init') #initialize telescope controller

    #the coordinates and targets are stored in root.time_of_individual_tracks, root.name_of_individual_tracks, root.individual_tracks 

    try:
        while keep_running_app:
            await trio.sleep(0.01) #checkpoint without blocking (e.g. GUI can operate now)

 #           root.console.text += 'Sitting on the beach'

 #           config.x = 'reset'
 #           root.console.text += config.x #this one gets slow fast

           # await trio.sleep(2)

            root.iPhone_stats.text = 'El vel' +str(config.el_speed)
            root.roboclaw_stats.text = 'Az vel' + str(config.az_speed)
            
            root.trackprogressbar.value = coordinate_loc
            root.trackprogressbar.max = max(1, len(config.selected_target_coords))
            
            ra,dec = get_ra_dec()
            root.ra_dec_display.text = str("%05.1f" % ra) + ' RA ' + str("%05.1f" % dec) + ' DEC'
            
            root.az_PV.text = str("%07.3f" % config.az_angle_PV) + ' (' + str("%08.3f" % config.az_pointing_error) + ') ' #update the angle display
            #print('######## el PV from config: ' + str(config.el_angle_PV))
            root.el_PV.text = str("%07.3f" % config.el_angle_PV) + ' (' + str("%08.3f" % config.el_pointing_error) + ')' #update the angle display

            #perform commands
            root.state_label.text = 'config.state = ' + config.state
            if config.state == 'idle':
                #update our stats
                root.current_target_stats.text='Current target: None'
                seq_num = 0
                coordinate_loc = 0
                stats, rc, address = telescope_control(rc, address = 0x80, cmd='noop', coord=[datetime.now(timezone.utc), config.az_angle_SV, config.el_angle_SV]) #update the PVs (done every telescope_control before command)
                #note: can't use the textbox values because if they are deleted they will be invalid when this command gets fired off in the background
                continue
            elif config.state == 'manual':
                #send command for manual move
                print('manual command issued')
                # config.az_angle_SV = float(root.manual_Az.text) #these don't get set anywhere else if a track isn't run first
                # config.el_angle_SV = float(root.manual_El.text)
                stats, rc, address = telescope_control(rc, address = 0x80, cmd='move_az_el', coord=[datetime.now(timezone.utc), float(root.manual_Az.text), float(root.manual_El.text)]) #initialize telescope controller
                config.state = 'idle' #don't keep hammering the roboclaw with the same position

#something fucky is going on here, need to get the right target list (not root.select_target as that's just the dropdown parent. Use config. instead?)
            elif config.state == 'track':    #send the commands to the telescope for tracking
                if do_once:
                    #set the starting track based on the dropdown
                    seq_num = config.selected_target_dropdown
                    print('root seq num ' + str(seq_num))
                    do_once = False
                    
                if len(config.selected_target_coords) < 1:
                    print('Need to compute some targets first')
                    config.state = 'idle'
                    continue

                if coordinate_loc >= len(config.selected_target_coords): #if we have completed a track
                    print('Completed track')
                    if root.automanual_checkbox.active is True: #if we aren't in manual track selection
                        #move to next track
                        seq_num = seq_num + 1 #move to the next one
                        print('seq num ' + str(seq_num))
                        print('len(root.individual_tracks)' + str(len(root.individual_tracks)))


                        if seq_num < len(root.individual_tracks): #make sure we have a next one to move to
                            seq_num_key = list(root.individual_tracks.keys())[seq_num] #get the key from the sequence number
                            print('seq_num_key ' + str(seq_num_key))

                            print('selecting next track:')
                            print('{0: <20}'.format(str(root.name_of_individual_tracks[seq_num])) + ' ' +  str(root.time_of_individual_tracks[seq_num])) 
                            print('seq_num_key ' + str(seq_num_key))

                            print(root.individual_tracks[seq_num_key])
                            config.selected_target_coords = root.individual_tracks[seq_num_key]
                            print('need to change the dropdowns for the new target selections here')
                            root.dropdown.value = seq_num
 ####mgiht be redundant                           root.current_target_stats.text='Current target: ' +root.individual_tracks[0].target_name

                            coordinate_loc = 0 #restart for the new track
                            continue
                        else:
                            print('Ran out of tracks, moving to idle state')
                            root.current_target_stats.text='Current target: None'
                            # seq_num = 0
                            # coordinate_loc = 0
                            config.state = 'idle'
                            do_once = True
                            continue
                    else:
                        config.state = 'idle'
                        do_once = True
                        # seq_num = 0
                        # coordinate_loc = 0
                        root.current_target_stats.text='Current target: None'

                        continue

                #need to check that we haven't run out of coordiantes first before moving to new coordinate
                config.current_coord_loc = coordinate_loc #update so others can see global variable
                config.selected_target_coords
                
                if datetime.now(timezone.utc) < config.selected_target_coords[coordinate_loc].time: #this will keep the next coordinates from being sent from the caller
                    stats, rc, address = telescope_control(rc, address=0x80,  cmd='noop', coord = [datetime.now(timezone.utc), float(root.manual_Az.text),float(root.manual_El.text)]) #update the scope PVs

                    if datetime.now(timezone.utc) > last_print:
                        last_print = datetime.now(timezone.utc) + timedelta(seconds=time_print_spacing)
                        print('Waiting for new coordinate in ' + str(config.selected_target_coords[coordinate_loc].time - datetime.now(timezone.utc) )+' H:M:S')

                    continue #wait for coordinate time to occur
                else:
                    root.current_target_stats.text='Current target: ' +config.selected_target_coords[0].target_name
                    #print('sending ' + str(coordinate_loc) + ' of ' + str(len(config.selected_target_coords)) + 'coords for ' + config.selected_target_coords[0].target_name)           

                    config.current_coord_loc = str(coordinate_loc)

                    config.coordinates = config.selected_target_coords[coordinate_loc] #store the coordinates in global variable
                    #set the displays
                    root.manual_Az.text = str(config.selected_target_coords[coordinate_loc].az)
                    root.manual_El.text = str(config.selected_target_coords[coordinate_loc].el)

                    #to do: compute az and el angular velocities to the future coordinate and pass those through to telescope_control too
                    #telescope control then commands a future distance n degrees away from the actual coordinate to go to, at the correct velocity
                    #this keeps the roboclaw from initiating deceleration once it gets close to the commanded point, but before the next one is sent
                    #add a bit in the telescope control to let it travel quickly if > 5 degrees away from the target or something
                    
                    
                    stats, rc, address = telescope_control(rc, address=0x80,  cmd='move_az_el', coord = [datetime.now(timezone.utc), float(root.manual_Az.text),float(root.manual_El.text), config.selected_target_coords[coordinate_loc].az_vel, config.selected_target_coords[coordinate_loc].el_vel], lookahead = False) #tell telescope to stop where it is
                    coordinate_loc = coordinate_loc+1
                    
                    #need to take the roboclaw stats to get PV
                    
            elif config.state == 'home':
                print('Sending homing sequence')

                #send the home command
                
            elif config.state == 'point_and_shoot':
                #bug: will not be able to access coord if still movign when start is triggered
                
                if config.az_in_position and  config.el_in_position:   #move to coordinate (first time will be in position because it isn't slewing somewhere)
                    sleep(1) #let any jitter settle
                    #take photo
                    if len(config.selected_target_coords)>0:    
                        coord = config.selected_target_coords.pop(0)
                        print('heading to ' + str(coord.az) + ',' + str(coord.el))
                        stats, rc, address = telescope_control(rc, address=0x80,  cmd='move_az_el', coord = [coord.time, coord.az,coord.el], lookahead = False) #tell telescope to stop where it is
                    else:
                        #remove first photo (taken at start arbitrary position)
                        config.state = 'idle'
                        print('Done with point and shoot')
                else:
                    try: #catch the first run where coord doesn't exist yet, but we want to not put in dummy starting points into the display
                        stats, rc, address = telescope_control(rc, address=0x80, coord = [datetime.now(timezone.utc), 0, 0, 0, 0], cmd='noop') #tell telescope to stop where it is
                        #print('Az error: '+ str(config.az_pointing_error) + 'el error: ' + str(config.el_pointing_error))
                        continue
                    except:
                        stats, rc, address = telescope_control(rc, address=0x80, coord = [datetime.now(timezone.utc), config.az_angle_PV, config.el_angle_PV], cmd='move_az_el') #stop the scope where it is
                        print('probably no coordinate due to unsettled start')
                        continue

                    
                
                # for coord in config.selected_target_coords:
                #     print('heading to coord az ' + str(coord[1]) + ' el ' + str(coord[2]))

                #     stats, rc, address = telescope_control(rc, address=0x80,  cmd='move_az_el', coord = coord, lookahead = False) #tell telescope to stop where it is
                #     if not config.az_in_position and not config.el_in_position:
                #         stats, rc, address = telescope_control(rc, address=0x80,  cmd='noop', coord = coord, lookahead = False) #tell telescope to stop where it is
                #         print('Az error: '+ str(config.az_pointing_error) + 'el error: ' + str(config.el_pointing_error))
                #         continue
                #     else :
                #         print('settled')
                #         sleep(1)#let any jitters settle
                #         #take photo here
                
                # config.state = 'idle'
                    
                    
                #time.sleep(1) # Sleep to let jitter settle
                #add photo request to camera queue
    #            queues['camera_q'].put('point_and_shoot_'+str(az)+'_'+str(el))
    #            #wait until queue size is zero (photo taken) before moving on
    #            while not config.camera_ready:
    #                pass #wait for camera to finish
                
                    
                
            else:
                print('unknown state')
                print(config.state)
            console_str = ('Al SV: ' + str(config.az_angle_SV) + ' PV: ' + str(config.az_angle_PV) + ' enc: ' + str(config.az_encoder_value) + '\n' +
                                'El SV: ' + str(config.el_angle_SV) + ' PV: ' + str(config.el_angle_PV) + ' enc: ' + str(config.el_encoder_value) +'\n')
            if len(config.selected_target_coords)>1 and coordinate_loc < len(config.selected_target_coords):
                console_str = console_str + str(config.selected_target_coords[coordinate_loc])
            root.console.text = console_str


    except trio.Cancelled as e:
        print('Wasting time was canceled', e)
        raise
    finally:
        # when canceled, print that it finished
        print('Done wasting time')

# async def read_log(root, nursery):
#     '''This method, which runs Kivy, is run by trio as one of the coroutines.
#     '''
#     try:
#         while True:
#             print('Reading console log')
#             await trio.sleep(2)
            
#             root.console.text = 

#             await trio.sleep(0) #checkpoint without blocking (e.g. GUI can operate now)
#     except trio.Cancelled as e:
#         print('Console log cancelled', e)
#     finally:
#         # when canceled, print that it finished
#         print('Done with console log')



class FloatInput(TextInput): #text input that only accepts numbers

    pat = re.compile('[^0-9]')
    def insert_text(self, substring, from_undo=False):
        pat = self.pat
        if '.' in self.text:
            s = re.sub(pat, '', substring)
        else:
            s = '.'.join(
                re.sub(pat, '', s)
                for s in substring.split('.', 1)
            )
        return super().insert_text(s, from_undo=from_undo)

class MainScreen(BoxLayout):
    
    
#size_hint: defines the size of a widget as a fraction of the parent’s size. Values are restricted to the range 0.0 - 1.0, e.g. 0.01 = 1/100th of the parent’s size (1%) and 1.0 = same size as the parent (100%).
#pos_hint: is used to place the widget relative to the parent.
#The size_hint and pos_hint are used to calculate a widget’s size and position only if the value(s) are not set to None. If you set these values to None, the layout will not position/size the widget and you can specify the values (x, y, width, height) directly in screen coordinates.

    def __init__(self, nursery = None, **kwargs):
        super(MainScreen, self).__init__(**kwargs)
        self.nursery = nursery #used to close window = stop code
        
        #layout
        #top level: controls (buttons)
        #mid level: text display (console)
        #bottom level: plot and images
        
        
        #self.cols = 1
        self.orientation = 'vertical' #set the direction of the MainScreen BoxLayout top level layout
        #self.top_level_layout = BoxLayout(orientation='vertical') #the overall window container
        #self.top_level_layout.add_widget(Label(text='top box'))
        #self.add_widget(self.top_level_layout)
        
        self.top_box = BoxLayout(orientation='vertical', size_hint_y = 0.375)    #the boxes that hold the top, middle and bottom sections
        self.mid_box = BoxLayout(orientation='horizontal', size_hint_y = 0.25) #two text fields (nominally track data and console output)
        self.bot_box = BoxLayout(orientation='horizontal', size_hint_y = 0.375) #two images, star field and image from camera

        self.add_widget(self.top_box)      #add the three sub-levels to the overall container
        self.add_widget(self.mid_box)
        self.add_widget(self.bot_box)

        self.top_box_row1 = BoxLayout(orientation='horizontal') #sub-boxes to hold button control to the top section
        self.top_box_row2 = BoxLayout(orientation='horizontal')
        self.top_box_row3 = BoxLayout(orientation='horizontal')
        self.top_box_row4 = BoxLayout(orientation='horizontal')
        self.top_box_row5 = BoxLayout(orientation='horizontal')
        self.top_box_row6 = BoxLayout(orientation='horizontal')
        self.top_box_row7 = BoxLayout(orientation='horizontal')
        #self.top_box_row8 = BoxLayout(orientation='horizontal')
        #self.top_box_row9 = BoxLayout(orientation='horizontal')

        #self.spacer = BoxLayout(orientation='horizontal')


        self.top_box.add_widget(self.top_box_row1)
        self.top_box.add_widget(self.top_box_row2)
        self.top_box.add_widget(self.top_box_row3)
        self.top_box.add_widget(self.top_box_row4)
        self.top_box.add_widget(self.top_box_row5)
        self.top_box.add_widget(self.top_box_row6)
        self.top_box.add_widget(self.top_box_row7)
        

        #self.top_box.add_widget(self.top_box_row8)
        #self.top_box.add_widget(self.top_box_row9)

        #self.top_box.add_widget(self.spacer)

        target_pairs = []
        self.time_of_individual_tracks = tuple()
        self.name_of_individual_tracks = []
        self.individual_tracks = dict() 

        
        #input_type is an OptionsProperty and defaults to ‘null’. Can be one of ‘null’, ‘text’, ‘number’, ‘url’, ‘mail’, ‘datetime’, ‘tel’ or ‘address’.
        #TextInput(hint_text='Kgs', input_filter = 'float', multiline=False, write_tab=False) #numeric only, tab moves to next object instead of writing \tab

# def on_focus(instance, value):
#     if value:
#         print('User focused', instance)
#     else:
#         print('User defocused', instance)

# textinput = TextInput()
# textinput.bind(focus=on_focus)


# def on_text(instance, value):
#     print('The widget', instance, 'have:', value)

# textinput = TextInput()
# textinput.bind(text=on_text)

        def manual_Az_focus_change(instance, value):
            if self.manual_Az.text.strip(): #check against empty inputs
                self.manual_Az.text = str( float(self.manual_Az.text) %360)
            else:
                self.manual_Az.text = '0.0'

        self.manual_Az = TextInput(multiline=False, hint_text='Manual Az (deg)', input_filter = 'float', write_tab = False) #input for numeric manual Az degrees
        self.manual_Az.text = '0.0'
        #self.manual_Az.bind(text = manual_Az_inputchange) #executes callback any time text is changed
        self.manual_Az.bind(focus = manual_Az_focus_change) #executes callback any time focus is changed


        def manual_El_focus_change(instance, value):
            if self.manual_El.text.strip(): #check against empty inputs
                self.manual_El.text = str( float(self.manual_El.text) %360) #Limit this to some bounds once limit switches/range of motion determined
                #self.manual_El.text = str(max(0, min(180, float(self.manual_El.text))))

            else:
                self.manual_El.text = '0.0'
            

        self.manual_El = TextInput(multiline=False, hint_text='Manual El (deg)', input_filter = 'float', write_tab = False) #input for numeric manual El degrees
        self.manual_El.text = '0.0'
        self.manual_El.bind(focus = manual_El_focus_change)


        def on_manual_AzEl_enter(instance): #button to accept manual typed Az/El
            #print('User pressed enter in', instance)
            #coerce output to 0-360
            self.manual_Az.text = str( float(self.manual_Az.text) %360)
            self.manual_El.text = str( float(self.manual_El.text) %360)

            print('Az : ' + self.manual_Az.text)
            print('El : ' + self.manual_El.text)
            config.state = 'manual'
        self.manual_AzEl = (Button(text='Move AzEl'))
        self.manual_AzEl.bind(on_press=on_manual_AzEl_enter)

#self.mytext.bind(text = self.calc)


    #     #Clock.schedule_interval(self.update, 1) #can use for timing updates
    #     return self.layout

    # def update(self, *args):
    #     self.name.text = str(self.current_i)
    #     self.current_i += 1
    #     if self.current_i >= 50:
    #         Clock.unschedule(self.update)

        #manual control buttons
        def plus_Az_callback(instance):
            self.manual_Az.text = str(float(self.manual_Az.text) + 1)
            print('The button <%s> is being pressed' % instance.text)
            self.manual_AzEl.trigger_action(0.1) # argument is for how long button should be pressed

        def minus_Az_callback(instance):
            self.manual_Az.text = str(float(self.manual_Az.text) - 1)
            print('The button <%s> is being pressed' % instance.text)
            self.manual_AzEl.trigger_action(0.1) # argument is for how long button should be pressed

        def plus_El_callback(instance):
            #self.manual_El.text = str(float(self.manual_El.text) + 1)
            self.manual_El.text = str( max(0, min(180, float(self.manual_El.text) + 1)) )

            print('The button <%s> is being pressed' % instance.text)
            self.manual_AzEl.trigger_action(0.1) # argument is for how long button should be pressed

        def minus_El_callback(instance):
            #self.manual_El.text = str(float(self.manual_El.text) - 1)
            self.manual_El.text = str( max(0, min(180, float(self.manual_El.text) - 1)) )

            print('The button <%s> is being pressed' % instance.text)
            self.manual_AzEl.trigger_action(0.1) # argument is for how long button should be pressed

        self.plus_Az_button = Button(text='+Az')
        self.plus_Az_button.bind(on_press=plus_Az_callback)

        self.minus_Az_button = Button(text='-Az')
        self.minus_Az_button.bind(on_press=minus_Az_callback)

        self.plus_El_button = Button(text='+El')
        self.plus_El_button.bind(on_press=plus_El_callback)

        self.minus_El_button = Button(text='-El')
        self.minus_El_button.bind(on_press=minus_El_callback)

        ImageUrl = 'https://kivy.org/doc/stable/_static/logo-kivy.png'

        self.iPhone_stats = Label(text='iPhone boot text', size_hint_x = 0.2)
        self.roboclaw_stats = Label(text='roboclaw boot text', size_hint_x = 0.2)
        self.az_PV = Label(text='0.0')
        self.el_PV = Label(text='0.0')
        self.current_target_stats = Label(text='Current target:', size_hint_x = 0.25)
        
        self.duration_grid_layout = BoxLayout(orientation='horizontal')#GridLayout(cols = 2)
        self.duration = TextInput(hint_text='Duration in sec per track', input_filter = 'float', multiline=False, write_tab=False, size_hint_x = 0.33) #numeric only, tab moves to next object instead of writing \tab
        self.duration.text = "1000"
        self.duration_label = Label(text='Track duration')
        self.duration_grid_layout.add_widget(self.duration)
        self.duration_grid_layout.add_widget(self.duration_label)
        
        
        

        # def target_list(instance, value):
        #     print(config.x)
        #     self.target_list.text = str( float(self.manual_El.text) %360)
        #     config.x = 'manual el config.x'


        self.target_list_textinput = TextInput(multiline=False, hint_text='List of targets', write_tab = False) #input for numeric manual El degrees
        self.target_list_textinput.text = 'sun moon iss'

    
        self.scroll_console = ScrollView(scroll_type = ['bars', 'content'], bar_width = 13, bar_margin = 5, size_hint_y = None)
        #dynamic heights her epossible, but tricky. Doesn't display scrolling behavior unless the text input is bigger than the scrollview
#        self.tracks_display = TextInput(multiline=True, hint_text='Start times for tracks', input_filter = 'float', write_tab = False) #input for numeric manual Az degrees
        self.console = TextInput(multiline=True, hint_text='console', write_tab = False) #input for numeric manual Az degrees
        self.console.text = 'console\n' + '\n'.join(['more' + str(i) +'\n' for i in range(1,50)]) #'console\nmore1\nmore2\nmore3\nmore4\n' #str([str(x) for x in self.individual_tracks])

        self.scroll_console.add_widget(self.console)

        
        
###works up until here (dropdown creation)

        self.dropdown = DropDown()
        def dropdown_callback(instance, value):
            print('inside dropdown callback')
            setattr(self.select_target, 'text', value)
            print('dropdown value' + value)
            s = ' '.join(value.split()[2:]) #ignore the Value 1 part of the dropdown menu
            config.selected_target = s #store for reference by the camera or other functions
            config.selected_target_dropdown = int(' '.join(value.split()[1]))
            print('dropdown value' + str(config.selected_target_dropdown))
            config.selected_target_coords = self.individual_tracks[s] #the list of target coordinates

            output_str = '' 
            for line in self.individual_tracks[s]:
                output_str += str([str(x) for x in line]) + '\n'
            self.console.text = output_str #display all the coordinates as a string for reference
        
        self.dropdown.bind(on_select=dropdown_callback)

        self.point_and_shoot = Button(text='Point and shoot')
        def point_and_shoot_callback(instance):
            self.time_of_individual_tracks, self.name_of_individual_tracks, self.individual_tracks = compute_target_coords(target_list=[['user_defined_pos', 'point_and_shoot']]) #default parameters used for point and shoot (start/end coords kept in config file)
            print('length of self.individual_tracks ' + str(len(self.individual_tracks)))
            print(self.individual_tracks)
            config.selected_target_coords = self.individual_tracks['point_and_shoot_target_name']
            config.state = 'point_and_shoot'
        self.point_and_shoot.bind(on_release=point_and_shoot_callback)
        
        self.point_and_shoot_start = Button(text='PaS start')
        def point_and_shoot_start_callback(instance):
            config.point_and_shoot_start = [config.az_angle_PV , config.el_angle_PV]
            self.console.text = 'Point and shoot start = ' + str(config.az_angle_PV) + ',' + str(config.el_angle_PV) +'\n'
        self.point_and_shoot_start.bind(on_release=point_and_shoot_start_callback)
        
        self.point_and_shoot_end = Button(text='PaS end')
        def point_and_shoot_end_callback(instance):
            config.point_and_shoot_end = [config.az_angle_PV , config.el_angle_PV ]
            self.console.text = self.console.text + 'Point and shoot end = ' + str(config.az_angle_PV) + ',' + str(config.el_angle_PV) + '\n'

        self.point_and_shoot_end.bind(on_release=point_and_shoot_end_callback)


        def point_and_shoot_FOV_callback(instance):
            config.point_and_shoot_FOV = float(self.point_and_shoot_FOV.text)
        self.point_and_shoot_FOV = TextInput(hint_text='FOV degrees', input_filter = 'float', multiline=False, write_tab=False) #numeric only, tab moves to next object instead of writing \tab
        self.point_and_shoot_FOV.text = "1"
        self.point_and_shoot_FOV.bind(on_release=point_and_shoot_FOV_callback)
        
        self.point_and_shoot_FOV_layout = BoxLayout(orientation = 'horizontal')
        self.point_and_shoot_FOV_layout.add_widget(self.point_and_shoot_FOV)
        self.point_and_shoot_FOV_layout.add_widget(Label(text = 'FOV deg'))



        def enable_buttons(disable = True): #turns buttons on/off for tasks that take a long time e.g. compute targets     
            self.manual_Az.disabled = disable
            self.manual_El.disabled = disable
            self.manual_AzEl.disabled = disable
     
            self.plus_Az_button.disabled = disable
            self.plus_El_button.disabled = disable
            self.checkbox_grid_layout.disabled = disable# automanual_checkbox
     
            self.minus_Az_button.disabled = disable
            self.minus_El_button.disabled = disable
     
     
            self.target_list_textinput.disabled = disable
            self.compute_targets.disabled = disable
            self.select_target.disabled = disable
            
            
            self.home.disabled = disable
            self.engage_track.disabled = disable
            self.halt_telescope.disabled = disable

        def compute_targets():
            print('computing target coordinates')
            self.console.text = 'computing target coordinates'

            target_list_str = self.target_list_textinput.text
            print(target_list_str)
            
            #check to see if astronomy mode or satellite tracking mode has been specified, if not, run them all as satellite_tracking 
            target_list = target_list_str.split() #= 'moon astronomy polaris astronomy sun astronomy ISS satellite_tracking'
            if len(target_list) == 1: #only one item
                target_pairs = [[target_list[0], 'satellite_tracking']]
                print('Only one target, presume it is satellite_tracking')
                print(target_list)
            elif len(target_list) ==  0:
                print('no tracks found! adding the moon so the rest of the program doesn''tbreak')
                target_list = ['moon', satellite_tracking]
            elif len(target_list) >= 2:
                if target_list[1] == 'satellite_tracking' or target_list[1] == 'astronomy':
                    #target lookup method was specified, continue
                    
                    print(target_list)
                    target_pairs = []
                    a = iter(target_list)
                    for x,y in zip(a, a):
                        target_pairs.append([x,y])
                    print(target_pairs)
                else:
                    target_pairs = [[x, 'satellite_tracking'] for x in target_list]
#            time_of_individual_tracks, name_of_individual_tracks, individual_tracks = compute_target_coords(target_list = target_pairs)
            #add to dropdown target selector
            
            #pass target_pairs to target computing function here
            self.time_of_individual_tracks, self.name_of_individual_tracks, self.individual_tracks = compute_target_coords(target_list = target_pairs, coincident_times=False, duration = float(self.duration.text))# [['sun','astronomy'], ['moon','astronomy'], ['polaris','astronomy'], ['M42', 'astronomy'], ['ISS', 'satellite_tracking'] ])
            s = ''
            for obstime,name in zip(self.time_of_individual_tracks, self.name_of_individual_tracks):
                s += name + ' ' +  obstime.strftime("%Y-%m-%d %H:%M:%S") + '\n'
                
            self.tracks_display.text = s
                        
            #create a dropdown with 10 buttons
            print('len of target pairs')
            print(len(target_pairs))
            self.console.text = self.console.text + '\ncomputing target coordinates' + '\n' + str(len(target_pairs)) + '\n' + 'Done computing targets'

            self.dropdown.clear_widgets() #remove the old target widgets from the dropdown

            #create a dropdown with 10 buttons
            for index in range(len(target_pairs)):
                # When adding widgets, we need to specify the height manually
                # (disabling the size_hint_y) so the dropdown can calculate
                # the area it needs.
            
                target = Button(text='Value %d ' % index + target_pairs[index][0] , size_hint_y=None)
            
                # for each button, attach a callback that will call the select() method
                # on the dropdown. We'll pass the text of the button as the data of the
                # selection.
                target.bind(on_release=lambda btn: self.dropdown.select(btn.text))
            
                # then add the button inside the dropdown
                self.dropdown.add_widget(target)

            enable_buttons(disable = False) #re-enable buttons after this is done running

                   #Clock.schedule_once(enable_buttons(disable = True), -1)
                   #Clock.schedule_once(lambda dt: enable_buttons(disable=True), -1)         
            
        def compute_targets_callback(instance):
            #have to do it this way because UI updates don't occur until the end of hte callback
            enable_buttons(disable = True) #disable buttons while this is running
            Clock.schedule_once(lambda dt: compute_targets(), 0)         

            
        self.compute_targets = Button(text='Compute target coords')#, size_hint_y = None)
        self.compute_targets.bind(on_press=compute_targets_callback)


        def halt_telescope_callback(instance):
            print('halting telescope')
            #coerce output to 0-360
            self.manual_Az.text = str(config.az_angle_PV)
            self.manual_El.text = str(config.el_angle_PV)
            config.state = 'manual'

        self.halt_telescope = Button(text='Halt telescope')
        self.halt_telescope.bind(on_press=halt_telescope_callback)

        #choose new target


        # def select_target_callback(instance):
            
        #     #create a dropdown with 10 buttons
        #     print('len of target pairs')
        #     print(len(target_pairs))
                

        self.select_target = Button(text='Select target')#, size_hint_y= None) #this is the main button for the dropdown (which contains other buttons and is hidden)
        self.select_target.bind(on_release=self.dropdown.open)

        self.trackprogressbar = ProgressBar(max=1000, size_hint_x = 0.75)#max=len(THE LENGTH OF THE TRACK))
        
        self.trackprogressbar.value = 750
        self.trackprogressbar.max = 10000
        
        def engage_track_callback(instance):
            config.state = 'track'

        self.engage_track = Button(text='Engage track') #this is the main button for the dropdown (which contains other buttons and is hidden)
        self.engage_track.bind(on_release=engage_track_callback)

        def home_callback(instance):
            #config.state = 'home'
            
            filename = datetime.now().strftime("%Y-%m-%d %H_%M_%S.%f")+".txt"
            folder = '/Users/andrewmiller/telescope/logs'
            with open(os.path.join(folder, filename), "w") as text_file:
                header = ('datetime,az_SV_input,el_SV_input,az_SV,el_SV,az_PV,el_PV,az_dest_counts,el_dest_counts,az_speed_SV,el_speed_SV,az_speed_SV_counts,el_speed_SV_counts,config.az_pointing_error,config.el_pointing_error,lookahead' + '\n')
                text_file.write(header)
                text_file.write(config.log)
            print('Done writing log')
        self.home = Button(text='Home') #this is the main button for the dropdown (which contains other buttons and is hidden)
        self.home.bind(on_release=home_callback)


        def ok_button(instance):
            print('ok pressed')
            enable_buttons(disable = True)
            #placeholder to update star plot
            self.remove_widget(self.plot_window)
            self.plot_window = FigureCanvasKivyAgg(figure = plot_nearby_stars(), size_hint_x = None, size_hint_y = None, height = 100)
            self.plot_window.height = 500 #for reasons unknown have to put the plot in, then adjust the height or it breaks (plot becomes NoneType)
            self.plot_window.width = 1000 #for reasons unknown have to put the plot in, then adjust the height or it breaks (plot becomes NoneType)
            self.add_widget(self.plot_window)

            #updating just the graph makes it resize (below)
                        # self.plot_window.figure.clf()
            # self.plot_window.figure = plot_nearby_stars()
            # self.plot_window.draw()


        def cancel_button(instance):
            print('cancel pressed')
            enable_buttons(disable = False)
####        self.ok_cancel_grid_layout = GridLayout(rows = 2, height = self.console.height)
        self.ok_cancel_grid_layout = BoxLayout(orientation = 'horizontal')

        self.ok = Button(text='OK') #this is the main button for the dropdown (which contains other buttons and is hidden)
        self.ok.bind(on_release=ok_button)
        self.cancel = Button(text='Cancel') #this is the main button for the dropdown (which contains other buttons and is hidden)
        self.cancel.bind(on_release=cancel_button)
        self.ok_cancel_grid_layout.add_widget(self.ok)
        self.ok_cancel_grid_layout.add_widget(self.cancel)


        def engage_track_callback(instance):
            config.state = 'track'

        self.engage_track = Button(text='Engage track') #this is the main button for the dropdown (which contains other buttons and is hidden)
        self.engage_track.bind(on_release=engage_track_callback)


        self.checkbox_grid_layout = GridLayout(cols = 2)
        def on_auto_checkbox_active(checkbox, value):
            if value:
                print('The checkbox', self.automanual_checkbox, 'is active')
            else:
                print('The checkbox', self.automanual_checkbox, 'is inactive')
        self.automanual_checkbox = CheckBox(active = True, size_hint_x = 0.25)
        self.automanual_checkbox.bind(active=on_auto_checkbox_active)
        self.automanual_label = Label(text='Auto next target', size_hint_x = 0.75)
        self.checkbox_grid_layout.add_widget(self.automanual_checkbox)
        self.checkbox_grid_layout.add_widget(self.automanual_label)

        #self.scroll_tracks_display = ScrollView(size_hint_y= None, scroll_type = ['bars', 'content'], bar_width = 13, bar_margin = 5, height = 200)
        self.scroll_tracks_display = ScrollView(scroll_type = ['bars', 'content'], bar_width = 13, bar_margin = 5, size_hint_y = None)

        #dynamic heights her epossible, but tricky. Doesn't display scrolling behavior unless the text input is bigger than the scrollview
#        self.tracks_display = TextInput(multiline=True, hint_text='Start times for tracks', write_tab = False, height =400, size_hint_y = None) #input for numeric manual Az degrees
        self.tracks_display = TextInput(multiline=True, hint_text='Start times for tracks', write_tab = False) #input for numeric manual Az degrees
        #self.tracks_display.text = 'tracks and times\nmore1\nmore2\nmore3\nmore4\n' #str([str(x) for x in self.individual_tracks])
        self.tracks_display.text = 'tracks and times\n' + ''.join(['track' + str(i) +'\n' for i in range(1,50)]) #'console\nmore1\nmore2\nmore3\nmore4\n' #str([str(x) for x in self.individual_tracks])

        self.scroll_tracks_display.add_widget(self.tracks_display)

        self.state_label = Label(text = config.state)

        #myplot = plot_nearby_stars()
        # if myplot is not None: 
        #     print('myplot not None')
        self.plot_window = FigureCanvasKivyAgg(figure = plot_nearby_stars(),  height = 100)
        #self.plot_window = FigureCanvasKivyAgg(figure = plot_nearby_stars(), size_hint_x = None, size_hint_y = None, height = 100)

        #self.plot_window.figure=plot_nearby_stars()
        self.plot_window.height = 500 #for reasons unknown have to put the plot in, then adjust the height or it breaks (plot becomes NoneType)
        self.plot_window.width = 1000 #for reasons unknown have to put the plot in, then adjust the height or it breaks (plot becomes NoneType)

        self.button_pressed = 0
        def button0_callback(instance):
            self.button_pressed = 0
        def button1_callback(instance):
            self.button_pressed = 1
        def button2_callback(instance):
            self.button_pressed = 2
        def button3_callback(instance):
            self.button_pressed = 3
        def button4_callback(instance):
            self.button_pressed = 4
        def button5_callback(instance):
            self.button_pressed = 5
        def button6_callback(instance):
            self.button_pressed = 6
        def button7_callback(instance):
            self.button_pressed = 7
        def button8_callback(instance):
            self.button_pressed = 8
        def button9_callback(instance):
            self.button_pressed = 9





#add the widgets in order for gridlayout (big downside to using that one really)

        self.ra_dec_display = Label(text = '', size_hint_x = 0.2)
        self.az_PV.size_hint_x = 0.2
        self.el_PV.size_hint_x = 0.2
        self.state_label.size_hint_x = 0.2
        self.point_and_shoot_FOV_layout.size_hint_x = 0.2
        
        self.top_box_row1.add_widget(self.az_PV)
        self.top_box_row1.add_widget(self.el_PV)
        self.top_box_row1.add_widget(self.ra_dec_display)
        self.top_box_row1.add_widget(self.point_and_shoot_FOV_layout)
        #self.top_box_row1.add_widget(Label(text='<--- Az, El PVs'))
        self.top_box_row1.add_widget(self.state_label)





        self.top_box_row2.add_widget(self.manual_Az)
        self.top_box_row2.add_widget(self.manual_El)
        self.top_box_row2.add_widget(self.manual_AzEl)
        self.top_box_row2.add_widget(self.point_and_shoot_start)
        self.top_box_row2.add_widget(self.home)


        
        

        self.top_box_row3.add_widget(self.plus_Az_button)
        self.top_box_row3.add_widget(self.plus_El_button)
        self.top_box_row3.add_widget(self.checkbox_grid_layout)# automanual_checkbox)
        self.top_box_row3.add_widget(self.point_and_shoot_end)
        self.top_box_row3.add_widget(self.engage_track)




        self.top_box_row4.add_widget(self.minus_Az_button)
        self.top_box_row4.add_widget(self.minus_El_button)
        self.top_box_row4.add_widget(self.duration_grid_layout)
        self.top_box_row4.add_widget(self.point_and_shoot)
        self.top_box_row4.add_widget(self.halt_telescope)


        self.target_list_textinput.size_hint_x = 0.4
        self.compute_targets.size_hint_x = 0.2
        self.select_target.size_hint_x = 0.2
        self.ok_cancel_grid_layout.size_hint_x = 0.2
        self.top_box_row5.add_widget(self.target_list_textinput)
        self.top_box_row5.add_widget(self.compute_targets)
        self.top_box_row5.add_widget(self.select_target)
        self.top_box_row5.add_widget(self.ok_cancel_grid_layout)

        
        self.top_box_row6.add_widget(self.current_target_stats)
        self.top_box_row6.add_widget(self.trackprogressbar) #(try to expand over 2 columns, labels as placeholders until then)

        # self.add_widget(Label(text='password'))
                
        self.top_box_row7.add_widget(self.roboclaw_stats)
        self.top_box_row7.add_widget(self.iPhone_stats)
        button_size_hint = (1-0.4)/10
        self.button0 = Button(text='0', size_hint_x = button_size_hint)
        self.button1 = Button(text='1', size_hint_x = button_size_hint)
        self.button2 = Button(text='2', size_hint_x = button_size_hint)
        self.button3 = Button(text='3', size_hint_x = button_size_hint) 
        self.button4 = Button(text='4', size_hint_x = button_size_hint)
        self.button5 = Button(text='5', size_hint_x = button_size_hint)
        self.button6 = Button(text='6', size_hint_x = button_size_hint)
        self.button7 = Button(text='7', size_hint_x = button_size_hint)
        self.button8 = Button(text='8', size_hint_x = button_size_hint)
        self.button9 = Button(text='9', size_hint_x = button_size_hint)
        self.top_box_row7.add_widget(self.button0)
        self.top_box_row7.add_widget(self.button1)
        self.top_box_row7.add_widget(self.button2)
        self.top_box_row7.add_widget(self.button3)
        self.top_box_row7.add_widget(self.button4)
        self.top_box_row7.add_widget(self.button5)
        self.top_box_row7.add_widget(self.button6)
        self.top_box_row7.add_widget(self.button7)
        self.top_box_row7.add_widget(self.button8)
        self.top_box_row7.add_widget(self.button9)      
        self.button0.bind(on_release=button0_callback)
        self.button1.bind(on_release=button1_callback)
        self.button2.bind(on_release=button2_callback)
        self.button3.bind(on_release=button3_callback)
        self.button4.bind(on_release=button4_callback)
        self.button5.bind(on_release=button5_callback)
        self.button6.bind(on_release=button6_callback)
        self.button7.bind(on_release=button7_callback)
        self.button8.bind(on_release=button8_callback)
        self.button9.bind(on_release=button9_callback)


        #configure the scrollbox and content heights
        self.console.size_hint_y = None
        self.console.height = max(self.console.minimum_height, self.scroll_console.height)
        self.tracks_display.size_hint_y = None
        self.tracks_display.height = max(self.tracks_display.minimum_height, self.scroll_tracks_display.height)
        print('self.console.minimum_height ' + str(self.console.minimum_height) + ' self.scroll_console.height'  + str(self.scroll_console.height))
        self.scroll_tracks_display.size_hint_y = 1
        self.scroll_console.size_hint_y = 1

        self.mid_box.add_widget(self.scroll_tracks_display)
        self.mid_box.add_widget(self.scroll_console)
        
        # self.add_widget(Label(text='Target list (name, method)'))
        # self.add_widget(Label(text='Manual El (deg)'))
        # self.add_widget(Label(text='Manual Az (deg)'))


        #self.add_widget(Label(text='placeholder')) #self.add_widget(self.automanual_label)

        # self.add_widget(self.username)
        # self.add_widget(Label(text='User Name'))
        # self.add_widget(self.password)


        

        self.bot_box.add_widget(self.plot_window)
        self.bot_box.add_widget(AsyncImage(source=ImageUrl)) #display camera images


        Window.bind(on_request_close=lambda *args: nursery.cancel_scope.cancel())
        
    def on_window_close(self, *args):
        if self.nursery:
            self.nursery.cancel_scope.cancel() #closes code running in nursery
        return False

#sun astronomy moon astronomy m43 astronomy venus astronomy

# class MyApp(App):
#     print('Doing some stuff')
#     def build(self):
#         tele = MainScreen()    
#         tele.manual_El.text = '123'

#         return tele
    
# if __name__ == '__main__':
#     MyApp().run()
    
        

if __name__ == '__main__':
    async def root_func():
        '''trio needs to run a function, so this is it. '''
        #root = MainScreen()    #Builder.load_string(kv)  # root widget

        async with trio.open_nursery() as nursery:
            root = MainScreen(nursery = nursery)    #Builder.load_string(kv)  # root widget

            '''In trio you create a nursery, in which you schedule async
            functions to be run by the nursery simultaneously as tasks.

            This will run all two methods starting in random order
            asynchronously and then block until they are finished or canceled
            at the `with` level. '''
            nursery.start_soon(run_app_GUI, root, nursery)
            nursery.start_soon(run_telescope, root)
        app = App.get_running_app()
        if app:
            app.stop()
        #root.close()
    trio.run(root_func)

        
##### queue managing/shutdown below

    # size = queues['status_q'].qsize() #check that an image is in the queue (collects the camera and image save returns from their separate threads)
    # if size >= 1:
    #     #get status_q and display message
    #     status = status_q.get()
    #     statuses.append(status)
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
# consumer_image_save.join()
# print('consumer_image_save join() complete')
# consumer_camera.join()
# print('consumer_camera join() complete')
keep_running_app = False
print('All threads finished')
#import time
#time.sleep(5)
#os._exit(0) #there is a bug in running kivy from spyder that keeps the window hung open forever when clicking it closed: https://github.com/spyder-ide/spyder/issues/19057