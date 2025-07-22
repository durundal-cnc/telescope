# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 11:17:53 2025

@author: AndrewMiller
"""
import sys
sys.path.append(r'C:\Users\AndrewMiller\OneDrive - Global Health Labs, Inc\Desktop\roboclaw_python')
import roboclaw_3 as roboclaw
from roboclaw import Roboclaw
import math
import config # https://docs.python.org/3/faq/programming.html#how-do-i-share-global-variables-across-modules

#All target/location variables are in degrees and convereted to counts for interacting with the motors

def telescope_control(queues):

    
    #Start new thread for this with queue to put in commanded positions and return current state
    #this is a consumer of the telescope queue


    deadband_counts = 10
    az_counts_per_rev = 10000
    el_counts_per_rev = 10000
    az_backlash = 66
    el_backlash = 66
    Bit1_az = 0 #testing only
    Bit1_el = 0 #testing only
    stats = {'test':'test'}


    # #Windows comport name
    # rc = Roboclaw("COM9",38400)
    # #Linux comport name
    # #rc = Roboclaw("/dev/ttyACM0",115200)
    # roboclaw_success = rc.Open()
    # if not roboclaw_success:
    #     print('Error opening roboclaw')
    #     return
    # address = 0x80
    # version = rc.ReadVersion(address)
    # print(version)
    
    # #initialize
    # az_backlash = 0 #counts for backlash
    # el_backlash = 0
    # az_counts_per_rev = 1000
    # el_counts_per_rev = 1000
    # # set pin funcs SetPinFunctions
    # Bit1 = 0
    # # get vers
    # # get config
    # # set M1 enc (to 4294967295/2 = 2147483647? or does zero just work too)
    # # set M2 enc
    # # set M1 max amps #stall current 5600mA, 10mA units for current limit
    # #SetDeadBand
    # # set M2 max amps 
    # # get M1POS PID
    # # set M1POS PID
    # # get M2POS PID
    # # set M2POS PID
    # # set M1VEL PID
    # # set M2VEL PID
    # # #set encoder size to equal stage max from home sensor
    
    
    # # #read functions
    # # get bufs
    # # get M1ENC
    # # get M2ENC
    # # get M1speed
    # # get M2speed
    # # get M1 current
    # # get M2 current
    
    # stats = 'this is the roboclaw stats pull'
    
    
    # # #write functions
    # # setM1speedacceldecelpos
    # # setM2speedacceldecelpos
    
    
        
    az_encoder_value = 0 #counts, replace with read function
    el_encoder_value = 0 #counts, replace with read function
    
    #queue commands
    run = True #send the command to the controller
    telescope_q = queues['telescope_q']
    status_q = queues['status_q']

    while not config.end_program:
        ###### Begin status updates

        
        #read state
        # az_encoder_value = rc.ReadEncM1(address)
        # el_encoder_value = rc.ReadEncM2(address)
        # az_speed = rc.ReadSpeedM1(address)
        # el_speed = rc.ReadSpeedM2(address)
        
        
        
        output_str = ''
        #read status
        # The status byte tracks counter underflow, direction and overflow. The byte value represents:
        # Bit0 - Counter Underflow (1= Underflow Occurred, Clear After Reading)
        # Bit1 - Direction (0 = Forward, 1 = Backwards)
        # Bit2 - Counter Overflow (1= Underflow Occurred, Clear After Reading)
        
        old_az_dir = Bit1_az
        if old_az_dir == 0:
            old_az_dir == -1 #used for backlash
            
        old_el_dir = Bit1_el
        if old_el_dir == 0:
            old_el_dir == -1 #used for backlash
        

        config.az_angle = 360*az_encoder_value/az_counts_per_rev #update global angle value
        config.el_angle = 360*el_encoder_value/el_counts_per_rev #update global angle value
        
        
        #pull queue
        try:
            queue_cmd = telescope_q.get(timeout = 10) #can never be filled if e.g. sun pointing or other keep out not met
            #debug print('1 received cmd az ' + str(queue_cmd['az']) + ' el ' + str(queue_cmd['el']) + ' in deg')

            az_dest_counts = queue_cmd['az'] * az_counts_per_rev / 360 #this will always be 0-360 degrees from the trajectory planner
            el_dest_counts = queue_cmd['el'] * el_counts_per_rev / 360 #this will always be 0-180 degrees from the trajectory planner
            cmd = queue_cmd['cmd']
            #debug print('2 az_dest_counts ' + str(az_dest_counts) + ' el_dest_counts ' + str(el_dest_counts) )
            
        except Exception as e:
            print('Exception in telescope queue pulling')
            print(e)
            cmd = 'hold'
        
        ###### End status updates
            
        
        ###### Begin command processing
        if cmd == 'home_motors':
            print('homing motors')
            #Home motors
                #rotate into sensors at constant velocity
                #record where sensor pulsed
                #move to pulse location
                #set encoder to 0
    
        elif cmd == 'move_az_el':
            #move to new az/el
            
            #determine right direction to move from current position (don't take long path around circle)    
            #don't forget about cable wrap limits!!
            positive_move = (az_dest_counts - az_encoder_value + az_counts_per_rev) % az_counts_per_rev
            if positive_move < (az_counts_per_rev/2):
                #go past the wrap around point
                az = az_encoder_value + positive_move
                new_az_dir = 1
            else:
                #negative move
                az =  az_counts_per_rev*math.floor(az_encoder_value / az_counts_per_rev) + az_dest_counts
                new_az_dir = -1
            #el can't wrap around so not an issue for it
            if el_dest_counts > el_encoder_value:
                new_el_dir = 1
            else:
                new_el_dir = -1
    
            #include backlash if reversing any axis
            if new_az_dir != old_az_dir:
                #include backlash
                az_dest_counts = az_dest_counts + new_az_dir*az_backlash #            need right direction
            if new_el_dir != old_el_dir:
                #include backlash
                el_dest_counts = el_dest_counts + new_el_dir*el_backlash #            need right direction
            
        elif cmd == 'hold': #stop and maintain current position
            az_dest_counts = az_encoder_value
            el_dest_counts = el_encoder_value
            
        elif cmd == 'status':
            print('status only')
            # #do nothing, just return status below
            # ReadMainBatteryVoltage = rc.ReadMainBatteryVoltage(address)
            # print('ReadMainBatteryVoltage ' + str(ReadMainBatteryVoltage))
            # ReadLogicBatteryVoltage = rc.ReadLogicBatteryVoltage(address)
            # print('ReadLogicBatteryVoltage ' + str(ReadLogicBatteryVoltage))

            # DutyM1 = rc.DutyM1(address)
            # print('DutyM1 ' + str(DutyM1))

            # DutyM2 = rc.DutyM2(address)
            # print('DutyM2 ' + str(DutyM2))

            # ReadBuffers = rc.ReadBuffers(address)
            # print('ReadBuffers ' + str(ReadBuffers))

            # ReadPWMs = rc.ReadPWMs(address)
            # print('ReadPWMs ' + str(ReadPWMs))

            # ReadCurrents = rc.ReadCurrents(address)
            # print('ReadCurrents ' + str(ReadCurrents))

            # ReadM1VelocityPID = rc.ReadM1VelocityPID(address)
            # print('ReadM1VelocityPID ' + str(ReadM1VelocityPID))

            # ReadM2VelocityPID = rc.ReadM2VelocityPID(address)
            # print('ReadM2VelocityPID ' + str(ReadM2VelocityPID))

            # ReadMinMaxMainVoltages = rc.ReadMinMaxMainVoltages(address)
            # print('ReadMinMaxMainVoltages ' + str(ReadMinMaxMainVoltages))

            # ReadMinMaxLogicVoltages = rc.ReadMinMaxLogicVoltages(address)
            # print('ReadMinMaxLogicVoltages ' + str(ReadMinMaxLogicVoltages))

            # ReadPinFunctions = rc.ReadPinFunctions(address)
            # print('ReadPinFunctions ' + str(ReadPinFunctions))

            # GetDeadBand = rc.GetDeadBand(address)
            # print('GetDeadBand ' + str(GetDeadBand))

            # ReadTemp = rc.ReadTemp(address)
            # print('ReadTemp ' + str(ReadTemp))

            # ReadTemp2 = rc.ReadTemp2(address)
            # print('ReadTemp2 ' + str(ReadTemp2))

            # ReadError = rc.ReadError(address)
            # print('ReadError ' + str(ReadError))

            # ReadEncoderModes = rc.ReadEncoderModes(address)
            # print('ReadEncoderModes ' + str(ReadEncoderModes))

            # GetConfig = rc.GetConfig(address)
            # print('GetConfig ' + str(GetConfig))

            # ReadM1MaxCurrent = rc.ReadM1MaxCurrent(address)
            # print('ReadM1MaxCurrent ' + str(ReadM1MaxCurrent))

            # ReadM2MaxCurrent = rc.ReadM2MaxCurrent(address)
            # print('ReadM2MaxCurrent ' + str(ReadM2MaxCurrent))

            # ReadPWMMode = rc.ReadPWMMode(address)
            # print('ReadPWMMode ' + str(ReadPWMMode))

            stats = {'test':'test'}

        else:
            print('command not recognized')
            output_str = output_str + '\n' + 'command not recognized'
    
        ###### End command processing
        
        ###### Begin checks and issue motor drive commands
    
        #do every time tasks
        #encoder runaway check
        #need to account for big slews onto target
            #break down slews into small segments in generator?
            #don't enable check until tracking?
            #only start if encoder counts not changing?
                #if encoder value not changing more than n between reads && motor current > some threshold
        # n = 100
        # if az_encoder_value < az_dest + n: #shut down:
        #     #run = False
        #     output_str = output_str + '\n' + 'encoder runaway on az'
        # if el_encoder_value < el_dest + n: #shut down:
        #     #run = False
        #     output_str = output_str + '\n' + 'encoder runaway on el'
            
        try:
            #check if we are where we are supposed to be close enough
            if abs(az_dest_counts - az_encoder_value) < 360*deadband_counts/az_counts_per_rev:
                config.az_in_position = True
            else:
                config.az_in_position = False
            if abs(el_dest_counts - el_encoder_value) < 360*deadband_counts/el_counts_per_rev:
                config.el_in_position = True
            else:
                config.el_in_position = False
        except Exception as e:
            print('No destination found for telescope control to slew to!')
            print(e)
                
    
        if run:
            #send the new command destinations
            #print('setM1speedacceldecelpos sent')
            #setM1speedacceldecelpos()
            #setM2speedacceldecelpos()
            
            #testing start (virtual controller)
            if az_dest_counts > az_encoder_value:
                az_encoder_value = az_encoder_value + 50
                Bit1_az = 0
            elif az_dest_counts < az_encoder_value:
                az_encoder_value = az_encoder_value - 50
                Bit1_az = 1

            else:
                continue #already in location
                
            if el_dest_counts > el_encoder_value:
                el_encoder_value = el_encoder_value + 50
                Bit1_el = 0
            elif el_dest_counts < el_encoder_value:
                el_encoder_value = el_encoder_value - 50
                Bit1_el = 1
            else:
                continue #already in location
            output_str = output_str + str(Bit1_az) + ' '  + str(Bit1_el)
            #debug print('3 Az current ' + str(az_encoder_value) + ' Az dest ' + str(az_dest_counts) + ' El current ' + str(el_encoder_value) + ' El dest ' + str(el_dest_counts))
            #testing end
            
        ###### End checks and motor drive commands
    
        #return status
            #include queue length to make sure not backing up
        status_q.put({'stats':stats, 'output':output_str})
        #print('Finished telescope loop, continuing')

    print('Terminated satellite tracking thread')

