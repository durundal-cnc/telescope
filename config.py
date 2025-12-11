# -*- coding: utf-8 -*-
"""
Created on Fri Jun 13 11:29:37 2025

@author: AndrewMiller
"""

end_program = False #tells threads to exit

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

az_angle_PV = 0
el_angle_PV = 0
az_in_position = False
el_in_potiion = False
az_current = 0
el_current = 0
az_speed = 0
el_speed = 0
az_pointing_error = 0
el_pointing_error = 0
total_pointing_error = 0

error_conds = []

Bit1_az = 0 #testing only
Bit1_el = 0 #testing only
az_encoder_value = 0 #counts, replace with read function
el_encoder_value = 0 #counts, replace with read function

ReadMainBatteryVoltage = 0.0
ReadLogicBatteryVoltage = 0.0
ReadBuffers  = 0.0
ReadPWMs  = 0.0
ReadCurrents  = 0.0
ReadM1VelocityPID  = 0.0
ReadM2VelocityPID  = 0.0
ReadMinMaxMainVoltages  = 0.0#these are the set limits, doesn't appear to work? Only returns 0,0
ReadMinMaxLogicVoltages  = 0.0
ReadPinFunctions = 0.0
GetDeadBand = 0.0
ReadTemp  = 0.0
ReadTemp2 = 0.0
ReadError = 0.0
ReadEncoderModes = 0.0
GetConfig = 0.0
ReadM1MaxCurrent = 0.0
ReadM2MaxCurrent  = 0.0
ReadPWMMode = 0.0
