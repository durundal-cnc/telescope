# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 11:17:53 2025

@author: AndrewMiller
"""
import sys
sys.path.append(r'/Users/andrewmiller/telescope/roboclaw_python')
from roboclaw_3 import Roboclaw
import math
import config # https://docs.python.org/3/faq/programming.html#how-do-i-share-global-variables-across-modules
from datetime import datetime, date, timezone, timedelta

#All target/location variables are in degrees and convereted to counts for interacting with the motors

def is_set(x, n): #get a specific bit from a byte
    # a more bitwise- and performance-friendly version:
    return x & 1 << n != 0


def telescope_control(queues):
    #Start new thread for this with queue to put in commanded positions and return current state
    #this is a consumer of the telescope queue

    #M1 = az
    #M2 = el

    az_encoder_cpr = 4000
    az_gearbox = 100 #75:1 on pololu website
    az_belt = 60/20 #number of teeth on large pulley/number of teeth on small pulley
    
    # el_encoder_cpr = 48
    # el_gearbox = 74.8317777777 #75:1 on pololu website
    # el_belt = 80/10
    
    el_encoder_cpr = 4000 #1000 pulse/rev with quadrature
    el_gearbox = 100 #harmonic drive
    el_belt = 60/20
    
    az_counts_per_rev = az_encoder_cpr * az_gearbox * az_belt
    az_arcsec_per_count = 60*60 * 360/az_counts_per_rev #for reference
    el_counts_per_rev = el_encoder_cpr * el_gearbox * el_belt
    el_arcsec_per_count = 60*60 * 360/el_counts_per_rev #for reference
    
    deadband_counts = 10
    az_backlash = 0 #needs measurement
    el_backlash = 0 #needs measurement
    az_accel = 100000 #counts/sec^2
    az_speed_SV = 40000#counts/sec #max is 240k
    az_deccel = 100000#counts/sec^2
    el_accel = 100000#counts/sec^2
    el_speed_SV = 40000#counts/sec
    el_deccel = 100000#counts/sec^2

    stats = {'test':'test'}
    
    robo_connected = True
    run = True #send the command to the controller
    
    if not robo_connected:
        Bit1_az = 0 #testing only
        Bit1_el = 0 #testing only
        az_encoder_value = 0 #counts, replace with read function
        el_encoder_value = 0 #counts, replace with read function
        
    if robo_connected:

        restore_defaults = True
        reconnect = True
        while reconnect:
            print('Opening roboclaw USB com port')
            #Windows comport name
            #rc = Roboclaw("COM4",38400)
            #Linux comport name
            rc = Roboclaw("/dev/tty.usbmodem1101",38400)
            roboclaw_success = rc.Open()

            address = 0x80 #decimal 128
        
            version = rc.ReadVersion(address)
            print(version)
            if not roboclaw_success or not version[0]:
                print('Error opening roboclaw')
               # return
                
            #initialize
            if restore_defaults:
                rc.RestoreDefaults(address) #need to re-establish serial connection once using this
                restore_defaults = False
            else:
                reconnect = False


        # # set M1 & M2 enc (to 4294967295 max range, but is signed))
        print('Setting encoder to 0')
        rc.SetEncM1(address,0)
        rc.SetEncM2(address,0)
        
        print('Setting encoder directions')
        rc.SetM1EncoderMode(address, int('00100000', 2)) #reverse motor relative direction direction int('00100000', 2))
        rc.SetM2EncoderMode(address, int('00000000', 2)) #may need to adjust depending on how final wiring goes (so positive encoder and positive angle match)
        
        # # set M1&M2 max amps #stall current 5600mA, 10mA units for current limit
        print('Setting max current')
        rc.SetM1MaxCurrent(address,5000) # Current value is in 10ma units. To calculate multiply current limit by 10.
        rc.SetM2MaxCurrent(address,5000)
        
        #voltage limit setting in python doesn't seem to work
        #rc.SetMainVoltages(address, 130, 110) #100mV increments
        #rc.SetMaxVoltageMainBattery(address, 66) #volts * 5.12 #this doesn't work in python but does work on their motion application
        #rc.SetMinVoltageMainBattery(address, 25) #(volts-6)*5=val
        
        #configure limit switches
        #S3 0x01 = E-stop
        #S4,S5:
        #home(user) = 0x62
        #home(auto) = 0xE2
        #Home(Auto)/Limit(Fwd) = 0xF2
        #Home(User)/Limit(Fwd) = 0x72
        #Limit(Rev) = 0x22
        #Limit(Fwd) = 0x12
        #Limit(Both) = 0x32
        print('read pin functions')
        rc.ReadPinFunctions(address)
        print('set pin functions')
        rc.SetPinFunctions(address, 0x01, 0x72, 0x72) #S3 is e-stop, home only for az, home and limits for el

        print('Set velocity PID to 0')
        rc.SetM1VelocityPID(address, 0x00000000, 0x00000000, 0x00000000, 256875) #QPPS set by measuring motor velocity open loop, PID here defaults
        rc.SetM2VelocityPID(address, 0x00000000, 0x00000000, 0x00000000, 256875) #QPPS set by measuring motor velocity open loop, PID here defaults
        #have to set the velocity PID terms to zero or it enables a cascaded PIV-D loop. Roboclaw defaults to having some terms in the velocity loop
#open loop speed seems to max out at 30000 but has significant variability there - down to 20k, averaging maybe 25k
#QPPS from PIV autotune = 102187

#these probably need to be re-evaluated after the drives have been run in a bit

        # # get M1POS PID [D(4 bytes), P(4 bytes), I(4 bytes), MaxI(4 bytes),Deadzone(4 bytes), MinPos(4 bytes), MaxPos(4 bytes)]
        print('Read position PID')
        rc.ReadM1PositionPID(address)
        # # set M1POS PID
        P = 39 #starting values from autotune in motion studio, for harmonic drives by their lonesome
        I = 2
        D = 540
        deadzone = 4
        max_I = 2617
        min_pos = 0
        max_pos = 2000000
        print('Set position PID')
        rc.SetM1PositionPID(address, P, I, D, max_I, deadzone, min_pos, max_pos)
        
        # # get M2POS PID
        print('Read position PID')
        print(rc.ReadM2PositionPID(address))
        # # set M2POS PID
        P = 39 #starting values from autotune in motion studio, for harmonic drives by their lonesome
        I = 2
        D = 540
        deadzone = 0
        max_I = 2617
        min_pos = 0 #negative values here seem to cause undefined behavior (though not in motion studio)
        max_pos = 2000000
        print('Set position PID')
        rc.SetM2PositionPID(address, P, I, D, max_I, deadzone, min_pos, max_pos)
        
        # stats = 'this is the roboclaw stats pull'
        print('Read errors in roboclaw registers')
        status = rc.ReadError(address)[0]
        status_list = [{'Normal' : status & 0x000000},
        {'E_Stop' : status & 0x000001}, #this always seems to be set?
        {'Temperature_Error' : status & 0x000002},
        {'Temperature_2_Error' : status & 0x000004},
        {'Main_Voltage_High_Error' : status & 0x000008},
        {'Logic_Voltage_High_Error' : status & 0x000010},
        {'Logic_Voltage_Low_Error' : status & 0x000020},
        {'M1_Driver_Fault_Error' : status & 0x000040},
        {'M2_Driver_Fault_Error' : status & 0x000080},
        {'M1_Speed_Error' : status & 0x000100},
        {'M2_Speed_Error' : status & 0x000200},
        {'M1_Position_Error' : status & 0x000400},
        {'M2_Position_Error' : status & 0x000800},
        {'M1_Current_Error' : status & 0x001000},
        {'M2_Current_Error' : status & 0x002000},
        {'M1_Over_Current_Warning' : status & 0x010000},
        {'M2_Over_Current_Warning' : status & 0x020000},
        {'Main_Voltage_High_Warning' : status & 0x040000},
        {'Main_Voltage_Low_Warning' : status & 0x080000},
        {'Temperature_Warning' : status & 0x100000},
        {'Temperature_2_Warning' : status & 0x200000},
        {'S4_Signal_Triggered' : status & 0x400000},
        {'S5_Signal_Triggered' : status & 0x800000},
        {'Speed_Error_Limit_Warning' : status & 0x01000000},
        {'Position_Error_Limit_Warning' : status & 0x02000000}]
        [print(stat) for stat in status_list if list(stat.values())[0] != 0]

#############testing

        # # rc.ReadEncoderModes(address)
        # # rc.DutyM1(address, 0)
        # # rc.DutyM2(address, 0)
        # # rc.DutyM1(address, 10000)
        # # rc.DutyM2(address, 10000)
        # # rc.SpeedAccelDeccelPositionM1(address, 100, 1000, 100, -20000, 1) #(address, accel, speed, deccel, position, buffer) #need to play with accel/decel in final system for smooth motion
        #  #rc.SpeedAccelDeccelPositionM2(address, 100, 1000, 100, 40000, 1) #(address, accel, speed, deccel, position, buffer)
        # el_dest_counts = 0
        # el_speed_SV = 100
        # buffer = 1
        # rc.SpeedAccelDeccelPositionM2(address, el_accel, el_speed_SV, el_deccel, el_dest_counts, buffer)
        # rc.GetDeadBand(address)
        # rc.GetConfig(address)
        # rc.ReadBuffers(address)
        # rc.ReadM2PositionPID(address)
        # rc.ReadM2VelocityPID(address)
        # rc.ReadPWMMode(address)
        
        # while True:
        #     az_encoder_value = rc.ReadEncM1(address) #Receive: [Enc1(4 bytes), Status, CRC(2 bytes)]
        #     Bit1_az = is_set(az_encoder_value[1], 1)
        #     az_encoder_value = az_encoder_value[1]
        #     el_encoder_value = rc.ReadEncM2(address)
        #     Bit1_el = is_set(el_encoder_value[1], 1)
        #     el_encoder_value = el_encoder_value[1]
        #     print(str(az_encoder_value) + ' ' + str(el_encoder_value))

##############testing

        #queue commands
#        telescope_q = queues['telescope_q']
#       status_q = queues['status_q']

    while not config.end_program:
        config.telescope_process_ready = True
        ###### Begin status updates

        if robo_connected:
            #read state
            az_encoder_value = rc.ReadEncM1(address) #Receive: [Enc1(4 bytes), Status, CRC(2 bytes)]
            Bit1_az = is_set(az_encoder_value[1], 1)
            az_encoder_value = az_encoder_value[1]
            print('az encoder is ' + str(az_encoder_value))

            el_encoder_value = rc.ReadEncM2(address)
            Bit1_el = is_set(el_encoder_value[1], 1)
            el_encoder_value = el_encoder_value[1]
            print('el encoder is ' + str(el_encoder_value))

            az_speed_PV = rc.ReadSpeedM1(address)[1] #status bit here also can tell speed (duplicate of ReadEncM1/2)
            el_speed_PV = rc.ReadSpeedM2(address)[1]
            
            # # get M1 current, get M2 current (10mA increments)
            az_current, el_current = rc.ReadCurrents(address)[1:]
            
            config.az_current = az_current
            config.el_current = el_current

            #check for sun intrusion (mostly during slews, since the start location won't be part of the trajectory that gets checked on creation)
            #get sun location at time from  config.coordinates.append({'az':aa.az.deg, 'el':aa.alt.deg, 'sun_az':sun_az, 'sun_el':sun_el, 'time':t})
            #if close to sun:
            #    run = False
        
        output_str = ''
        #read status
        # The status byte tracks counter underflow, direction and overflow. The byte value represents:
        # Bit0 - Counter Underflow (1= Underflow Occurred, Clear After Reading)
        # Bit1 - Direction (0 = Forward, 1 = Backwards)
        # Bit2 - Counter Overflow (1= Underflow Occurred, Clear After Reading)
        
        old_az_dir = Bit1_az #see if we have reversed the motors, and add backlash compensation if so
        if old_az_dir == 0:
            old_az_dir == -1 #used for backlash
            
        old_el_dir = Bit1_el
        if old_el_dir == 0:
            old_el_dir == -1 #used for backlash
        

        config.az_angle = 360*az_encoder_value/az_counts_per_rev #update global angle value
        config.el_angle = 360*el_encoder_value/el_counts_per_rev #update global angle value
        
        if len(config.coordinates) > 0:
            print('Reading ' + str(len(config.coordinates)) + ' coordinates')


        #pull queue
        try:
            # ####get command using queues
            # queue_cmd = telescope_q.get(timeout = 100) #can never be filled if e.g. sun pointing or other keep out not met
            # #debug print('1 received cmd az ' + str(queue_cmd['az']) + ' el ' + str(queue_cmd['el']) + ' in deg')

            # az_dest_counts = queue_cmd['az'] * az_counts_per_rev / 360 #this will always be 0-360 degrees from the trajectory planner
            # el_dest_counts = queue_cmd['el'] * el_counts_per_rev / 360 #this will always be 0-180 degrees from the trajectory planner
            # cmd = queue_cmd['cmd']
            # #debug print('2 az_dest_counts ' + str(az_dest_counts) + ' el_dest_counts ' + str(el_dest_counts) )
            # ####finish command using queues
            
            
            ####get command using config.coordinates
            time_print_spacing = 1
        
            if len(config.coordinates) > 0:
                print('Pulling coordinate')
                coord = config.coordinates.pop(0)
                print(coord['time'])
                print(datetime.now())
                last_print = datetime.now(timezone.utc)
                
                while coord['time'] < datetime.now(timezone.utc): #go through the list until we reach the present if not there already
                    coord = config.coordinates.pop(0)
    
                while datetime.now(timezone.utc) < coord['time']:
                    if datetime.now(timezone.utc) > last_print:
                        last_print = datetime.now(timezone.utc) + timedelta(seconds=time_print_spacing)
                        print('Waiting for new coordinate in ' + str(coord['time'] - datetime.now(timezone.utc) )+' H:M:S')
                    continue #wait for coordinate time to occur
                az_dest_counts = coord['az'] * az_counts_per_rev / 360 #this will always be 0-360 degrees from the trajectory planner
                el_dest_counts = coord['el'] * el_counts_per_rev / 360 #this will always be 0-180 degrees from the trajectory planner
                print('Issuing new coordinates to az_dest_counts '+str(az_dest_counts) + ' and el_dest_counts ' + str(el_dest_counts))

                cmd = 'move_az_el'
            else:
                cmd = 'noop'

            
            ####finish command using config.coordinates
        except Exception as e:
            print('Exception in pulling new coordinates in telescope_control!')
            print(e)
            cmd = 'hold'
            #return
        
        ###### End status updates
            
        
        ###### Begin command processing
        if cmd == 'home_motors':
            print('homing motors')
            #Home motors
                #rotate into sensors at constant velocity
                #record where sensor pulsed
                #move to pulse location
                #set encoder to 0
    
        elif cmd == 'noop': #no op
            pass
    
        elif cmd == 'move_az_el':
            #move to new az/el
            
            #determine right direction to move from current position (don't take long path around circle)    
            #don't forget about cable wrap limits!!
            positive_move = (az_dest_counts - az_encoder_value + az_counts_per_rev) % az_counts_per_rev
            if positive_move < (az_counts_per_rev/2):
                #go past the wrap around point
                az_dest_counts = az_encoder_value + positive_move
                new_az_dir = 1
            else:
                #negative move
                az_dest_counts =  az_counts_per_rev*math.floor(az_encoder_value / az_counts_per_rev) + az_dest_counts
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
            if robo_connected:
                #do nothing, just return status below
                ReadMainBatteryVoltage = rc.ReadMainBatteryVoltage(address)
                print('ReadMainBatteryVoltage ' + str(ReadMainBatteryVoltage[1]/10))
                ReadLogicBatteryVoltage = rc.ReadLogicBatteryVoltage(address)
                print('ReadLogicBatteryVoltage ' + str(ReadLogicBatteryVoltage[1]/10))
    
                ReadBuffers = rc.ReadBuffers(address)
                print('ReadBuffers ' + str(ReadBuffers)) #A return value of 0x80(128) indicates
                #the buffer is empty. A return value of 0 indiciates the last command sent is executing.
    
                ReadPWMs = rc.ReadPWMs(address)
                print('ReadPWMs ' + str(ReadPWMs))
    
                ReadCurrents = rc.ReadCurrents(address)
                print('ReadCurrents ' + str(ReadCurrents[1]/1000) + ' ' + str(ReadCurrents[2]/1000))
    
                ReadM1VelocityPID = rc.ReadM1VelocityPID(address)
                print('ReadM1VelocityPID ' + str(ReadM1VelocityPID))
    
                ReadM2VelocityPID = rc.ReadM2VelocityPID(address)
                print('ReadM2VelocityPID ' + str(ReadM2VelocityPID))
    
                ReadMinMaxMainVoltages = rc.ReadMinMaxMainVoltages(address) #these are the set limits, doesn't appear to work? Only returns 0,0
                print('ReadMinMaxMainVoltages ' + str(ReadMinMaxMainVoltages[1]/10) + ' ' + str(ReadMinMaxMainVoltages[2]/10))
                
                ReadMinMaxLogicVoltages = rc.ReadMinMaxLogicVoltages(address)
                print('ReadMinMaxLogicVoltages ' + str(ReadMinMaxLogicVoltages[1]/10) + ' ' + str(ReadMinMaxLogicVoltages[2]/10))
    
                ReadPinFunctions = rc.ReadPinFunctions(address)
                print('ReadPinFunctions ' + str(ReadPinFunctions))
    
                GetDeadBand = rc.GetDeadBand(address)
                print('GetDeadBand ' + str(GetDeadBand))
    
                ReadTemp = rc.ReadTemp(address)
                print('ReadTemp ' + str(ReadTemp[1]/10))
    
                ReadTemp2 = rc.ReadTemp2(address)
                print('ReadTemp2 ' + str(ReadTemp2[1]/10))
    
                ReadError = rc.ReadError(address)
                print('ReadError ' + str(ReadError))
    
                ReadEncoderModes = rc.ReadEncoderModes(address)
                print('ReadEncoderModes ' + str(ReadEncoderModes))
    
                GetConfig = rc.GetConfig(address)
                print('GetConfig ' + str(GetConfig))
    
                ReadM1MaxCurrent = rc.ReadM1MaxCurrent(address)
                print('ReadM1MaxCurrent ' + str(ReadM1MaxCurrent[1]/100))
    
                ReadM2MaxCurrent = rc.ReadM2MaxCurrent(address)
                print('ReadM2MaxCurrent ' + str(ReadM2MaxCurrent[1]/100))
    
                ReadPWMMode = rc.ReadPWMMode(address)
                print('ReadPWMMode ' + str(ReadPWMMode))

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
            try:
                #send the new command destinations
                
                if robo_connected:
                    #The Buffer argument can be set to a 1 or 0. If a value of 0 is used the command will be buffered
                    #and executed in the order sent. If a value of 1 is used the current running command is stopped,
                    #any other commands in the buffer are deleted and the new command is executed
                    buffer = 1
                    print('az_accel' + str(az_accel))
                    print('az_speed' + str(az_speed_SV))
                    print('az_deccel' + str(az_deccel))
                    print('az_dest_counts' + str(az_dest_counts))
                    az_dest_counts = round(az_dest_counts) #processes as a float but controller only does ints
                    print('az_dest_counts' + str(az_dest_counts))
                    el_dest_counts = round(el_dest_counts)


                    #if distance > some constant, slew fast
                    #then
                    #get current point and new point timestamps and distances
                    #set speed to match
                    #drive to new point
                    #remove old point from list
                    #repeat
                    

                    rc.SpeedAccelDeccelPositionM1(address, az_accel, az_speed_SV, az_deccel, az_dest_counts, buffer) #(address, accel, speed, deccel, position, buffer)
                    rc.SpeedAccelDeccelPositionM2(address, el_accel, el_speed_SV, el_deccel, el_dest_counts, buffer)
                else:
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
            except Exception as e:
                print(e)
            
        ###### End checks and motor drive commands
    
        #return status
            #include queue length to make sure not backing up
        #status_q.put({'stats':stats, 'output':output_str})
        #print('Finished telescope loop, continuing')

    print('Terminated satellite tracking thread')

