# -*- coding: utf-8 -*-
"""
Created on Fri Jun 13 11:29:37 2025

@author: AndrewMiller
"""

end_program = False #tells threads to exit
az_angle = 0
el_angle = 0
az_in_position = False
el_in_potiion = False
camera_ready = False
x = 0 #use for quit debugging
image = ''
tracking_ready = False
camera_process_ready = False
image_process_ready = False
telescope_process_ready = False
coordinates = []
my_locs = {}
mag_declination = 0
compass_dir = 0

lat = 0
long = 0

az_current = 0
el_current = 0