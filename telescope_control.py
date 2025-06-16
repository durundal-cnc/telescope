# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 11:17:53 2025

@author: AndrewMiller
"""

def telescope_control(queues):
    import math
    import config # https://docs.python.org/3/faq/programming.html#how-do-i-share-global-variables-across-modules

    
    #Start new thread for this with queue to put in commanded positions and return current state
    #this is a consumer of the telescope queue
    import sys
    sys.path.append(r'C:\Users\AndrewMiller\OneDrive - Global Health Labs, Inc\Desktop\roboclaw_python')
    import roboclaw_3 as roboclaw
    from roboclaw import Roboclaw

    az_counts_per_rev = 1000
    el_counts_per_rev = 1000
    az_backlash = 100
    el_backlash = 100
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
    
    
        
    az_encoder_value = 0
    el_encoder_value = 0 #replace with read function
    
    #queue commands
    run = True #send the command to the controller
    telescope_q = queues['telescope_q']
    status_q = queues['status_q']

    while not config.end_program:
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
        queue_cmd = telescope_q.get()
        
        az_dest = queue_cmd['az']#this will always be 0-360 degrees (in counts) from the trajectory planner
        el_dest = queue_cmd['el'] #this will always be 0-180 degrees (in counts) from the trajectory planner
        
        cmd = queue_cmd['cmd']
        
        print('received cmd ' + cmd + ' az ' + str(az_dest) + ' el ' + str(el_dest))
        
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
            positive_move = (az_dest - az_encoder_value + az_counts_per_rev) % az_counts_per_rev
            if positive_move < (az_counts_per_rev/2):
                #go past the wrap around point
                az = az_encoder_value + positive_move
                new_az_dir = 1
            else:
                #negative move
                az =  az_counts_per_rev*math.floor(az_encoder_value / az_counts_per_rev) + az_dest
                new_az_dir = -1
            #el can't wrap around so not an issue for it
            if el_dest > el_encoder_value:
                new_el_dir = 1
            else:
                new_el_dir = -1
    
            #include backlash if reversing any axis
            if new_az_dir != old_az_dir:
                #include backlash
                az_dest = az_dest + new_az_dir*az_backlash #            need right direction
            if new_el_dir != old_el_dir:
                #include backlash
                el_dest = el_dest + new_el_dir*el_backlash #            need right direction
            
        elif cmd == 'hold': #stop and maintain current position
            az_dest = az_encoder_value
            el_dest = el_encoder_value
            
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
    
        #do every time tasks
        #encoder runaway check
        #need to account for big slews onto target
            #break down slews into small segments in generator?
            #don't enable check until tracking?
            #only start if encoder counts not changing?
                #if encoder value not changing more than n between reads && motor current > some threshold
        n = 100
        if az_encoder_value < az_dest + n: #shut down:
            #run = False
            output_str = output_str + '\n' + 'encoder runaway on az'
        if el_encoder_value < el_dest + n: #shut down:
            #run = False
            output_str = output_str + '\n' + 'encoder runaway on el'

    
        if run:
            #send the new command destinations
            print('setM1speedacceldecelpos sent')
            #setM1speedacceldecelpos()
            #setM2speedacceldecelpos()
            
            #testing start
            if az_dest > az_encoder_value:
                az_encoder_value = az_encoder_value + 10
                Bit1_az = 0
            elif az_dest < az_encoder_value:
                az_encoder_value = az_encoder_value - 10
                Bit1_az = 1

            else:
                continue #already in location
                
            if el_dest > el_encoder_value:
                el_encoder_value = el_encoder_value + 10
                Bit1_el = 0
            elif el_dest < el_encoder_value:
                el_encoder_value = el_encoder_value - 10
                Bit1_el = 1
            else:
                continue #already in location
            output_str = output_str + str(Bit1_az) + ' '  + str(Bit1_el)
            print('Az current ' + str(az_encoder_value) + ' El current ' + str(el_encoder_value))
            #testing end
    
        #return status
            #include queue length to make sure not backing up
        status_q.put({'stats':stats, 'output':output_str})


