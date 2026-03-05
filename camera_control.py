# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 09:58:25 2025

@author: AndrewMiller
"""
import config # https://docs.python.org/3/faq/programming.html#how-do-i-share-global-variables-across-modules
from datetime import datetime, date, timezone, timedelta
from PIL import Image #used for testing only
import os
os.environ["OPENCV_VIDEOIO_MSMF_ENABLE_HW_TRANSFORMS"] = "0" #https://github.com/opencv/opencv/issues/17687 must go before import cv2
import cv2 
import numpy as np
import time
import sys
import trio
from kivy.core.image import Image as CoreImage
from PIL import Image
from io import BytesIO
from kivy.graphics.texture import Texture

#NOTE: the cv2 CAP_DSHOW and CAP_MSMF are windows specific, need to rewrite if using Linux or Mac

async def camera_control(root):

    #initialize camera
    print('starting camera thread')
    
    
    try:
        #import imageio as iio
       # import imageio.v3 as iio
    
        #import imageio_ffmpeg
        #import os, os.path
        #import matplotlib.pyplot as plt
        
        #get camera number
        camera = 1
        print('using camera # ' + str(camera))
        whitebalance = 'auto'
        exposure = -4
        resolution = 'HD'
        fps = 30
        config.focus = 'auto'
        
        if resolution == '4K':
            size = (4096,2160) #ffmpeg will re-scale image from cameras so need to get native resolutions to avoid this
        elif resolution == 'HD':
            size = (1920,1080)
        elif resolution == 'SD':
            size = (640,480)

        #for idx, frame in enumerate(iio.imiter(camera, fps = fps, size=(1920,1080))):
        #for idx, frame in enumerate(iio.imiter(camera, fps = fps, size=size, '-fflags=nobuffer')):
    
            
        # define a video capture object 
        
        
        #for the Logitech Brio 4K camera
        #CAP_MSMF supports setting the white balance but not auto white balance, CAP_DSHOW supports auto white balance but not setting it manually
        if whitebalance == 'auto':     
            vid = cv2.VideoCapture(camera, cv2.CAP_DSHOW   )  #CAP_MSMF  avoids delay loading the VideoCapture object
            wb = vid.set(cv2.CAP_PROP_AUTO_WB, 1)
        else:
            vid = cv2.VideoCapture(camera, cv2.CAP_MSMF   )  #CAP_MSMF  avoids delay loading the VideoCapture object
            wb = vid.set(cv2.CAP_PROP_AUTO_WB, 0)
            wb = vid.set(cv2.CAP_PROP_TEMPERATURE, whitebalance) #whitebalance is in Kelvin min=2800 max=6500
        print('White balance set: ' + whitebalance + '='+str(wb))
        if exposure == 'auto':        
            #note: auto-exposure can slow down camera considerably (at least for built in camera)
            ex = vid.set(cv2.CAP_PROP_AUTO_EXPOSURE, 1)
        else:
            ex = vid.set(cv2.CAP_PROP_AUTO_EXPOSURE, 0)
            ex = vid.set(cv2.CAP_PROP_EXPOSURE, exposure)#I believe the values are negative power of 2, e.g. -4 would be 1/16
        print('exposure set: ' + str(exposure) + ' success ' + str(ex))
    
        vid.set(cv2.CAP_PROP_BUFFERSIZE, 1) #keep only the most recent frame
        vid.set(cv2.CAP_PROP_FRAME_WIDTH, size[0])
        vid.set(cv2.CAP_PROP_FRAME_HEIGHT, size[1])
        vid.set(cv2.CAP_PROP_FPS, fps)
        vid.set(cv2.CAP_PROP_FOURCC, cv2.VideoWriter_fourcc('M', 'J', 'P', 'G')) #supposed to give good frame rates
        #vid.set(cv2.CAP_PROP_FOURCC, cv2.VideoWriter_fourcc(*"MJPG")) # add this line

    
        i = 0
        old_focus = '' #initialize
        old_exposure = ''
        old_whitebalance = ''
        last_photo = datetime.now(timezone.utc)
        photo = False
        attempt_focus = True #not all cameras support, so disable so timeout doesn't slow it down
        attempt_exposure = True
        attempt_whitebalance = True
        
        
        
        while(config.end_program is False):
            await trio.sleep(0.01) #checkpoint without blocking (e.g. GUI can operate now)

            config.camera_process_ready = True #signal loop is ready
            i = i + 1
            
            focus = config.focus #get the autofocus or manual focus value
            if  old_focus != focus and attempt_focus: #if change to focus setting
                if focus == 'auto':
                    focus_success = vid.set(cv2.CAP_PROP_AUTOFOCUS, 1) # turn the autofocus on
                else: 
                    focus_success1 = vid.set(cv2.CAP_PROP_AUTOFOCUS, 0) # turn the autofocus off
                    focus_success2 = vid.set(cv2.CAP_PROP_FOCUS , focus)  #set the focus,  min: 0, max: 255, increment:5
                    focus_success = focus_success1 and focus_success2
                if focus_success:
                    old_focus = focus
                else:
                    attempt_focus = False
            
            exposure = config.exposure #get the autoexposure or manual exposure value
            if  old_exposure != exposure and attempt_exposure: #if change to exposure setting
                if exposure == 'auto':
                    exposure_success = vid.set(cv2.CAP_PROP_AUTO_EXPOSURE, 1) # turn the autoexposure on
                else: 
                    exposure_success1 = vid.set(cv2.CAP_PROP_AUTO_EXPOSURE, 0) # turn the autoexposure off
                    exposure_success2 = vid.set(cv2.CAP_PROP_EXPOSURE , exposure)  #set the exposure,  min: 0, max: 255, increment:5
                    exposure_success = exposure_success1 and exposure_success2
                if exposure_success:
                    old_exposure = exposure
                else:
                    attempt_exposure = False
                    
                    
            whitebalance = config.whitebalance #get the autowhitebalance or manual whitebalance value
            if  old_whitebalance != whitebalance and attempt_whitebalance: #if change to whitebalance setting
                if whitebalance == 'auto':
                    whitebalance_success = vid.set(cv2.CAP_PROP_AUTO_WB, 1) # turn the autowhitebalance on
                else: 
                    whitebalance_success1 = vid.set(cv2.CAP_PROP_AUTO_WB, 0) # turn the autowhitebalance off
                    whitebalance_success2 = vid.set(cv2.CAP_PROP_TEMPERATURE , whitebalance)  #set the whitebalance,  min: 0, max: 255, increment:5
                    whitebalance_success = whitebalance_success1 and whitebalance_success2
                if whitebalance_success:
                    old_whitebalance = whitebalance
                else:
                    attempt_whitebalance = False
            
            
            if config.take_photo:
                photo = True
            
            if config.take_photos:
                if datetime.now(timezone.utc) >= last_photo + timedelta(seconds = config.sec_between_photos):
                    photo = True
                    last_photo = datetime.now(timezone.utc)
            #capture and display photos continuously on screen
            # print('camera read 1'+ str(round(time.time() * 1000)))
            ret, frame = vid.read() #double read is a hack to clear buffer (only for when doing long timelapse stuff where the stream isn't polled regularly)
            # print('camera read 2'+ str(round(time.time() * 1000)))
            #ret, frame = vid.read() 
            # print('camera read 3'+ str(round(time.time() * 1000)))

            # print('im start'+ str(round(time.time() * 1000)))
            # print('size:' + str(frame.shape) + '    ret: '+ str(ret) + '     ' + str(attempt_focus) + '       i=' + str(i)) #tells labview the photo size for scaling
            # im = Image.fromarray(frame) #display the image in Kivy
            # print('im end'+ str(round(time.time() * 1000)))
            
            # print('byte_io start'+ str(round(time.time() * 1000)))
            # byte_io = BytesIO()
            # print('byte_io end'+ str(round(time.time() * 1000)))
            
            # print('new process start'+ str(round(time.time() * 1000)))
            w, h, _ = frame.shape
            if camera == 0: #built in webcam on Dell laptop
                buf = bytes(np.flipud(frame))
            elif camera == 1: #Logi Brio
                buf = bytes(frame)

            texture = Texture.create(size=(h, w))
            texture.blit_buffer(buf, colorfmt='rgb', bufferfmt='ubyte')
            # print('new process end'+ str(round(time.time() * 1000)))
            
            
            
            # print('save start'+ str(round(time.time() * 1000)))
            # im.save(byte_io, format='PNG')
            # print('save end'+ str(round(time.time() * 1000)))
            
            # print('pointer start'+ str(round(time.time() * 1000)))
            # byte_io.seek(0)
            # print('pointer end'+ str(round(time.time() * 1000)))

            # print('texture start'+ str(round(time.time() * 1000)))
            # root.displayed_image.texture = CoreImage(byte_io, ext='png').texture
            root.displayed_image.texture = texture

            # print('texture end'+ str(round(time.time() * 1000)))

            # print('reload start'+ str(round(time.time() * 1000)))
            root.displayed_image.reload()
            # print('reload end'+ str(round(time.time() * 1000)))

            # print('old focus: ' + str(old_focus) + ' focus: ' + str(focus))

            #consider reducing size of displayed graphic to make faster? seems very slow and weird though)
            
            if photo: #take photo
                print('taking photo')
                photo = False
                config.take_photo = False #reset for single shots (for continuous the toggle has to be disabled by the user)
    
                filename = datetime.now().strftime('%Y-%m-%d_%H-%M-%S.%f')
                config.photos.append([filename, frame]) #add photo and timestamp to the list of photos for display (sky_scan) and storage (image_save)
                
        print('close windows and camera')
        # After the loop release the cap object 
        vid.release() 
        # Destroy all the windows 
        cv2.destroyAllWindows() 
        
        print('end')
        #plt.imshow(screenshot)
    except Exception as e:
        print('Error occured in camera_control: ')
        print(e)

# You might be requesting a camera output that is not native to your hardware, causing ffmpeg to work to encode to your request.

# Try this to get a list of your USB devices: ffmpeg -list_devices true -f dshow -i dummy

# Then use this with one of your cameras from the above output: ffmpeg -list_options true -f dshow -i video="Name of your camera here"

# This gives you the vcodec, video_size, framerate, and video:audio parameters to use in your ffmpeg command line.
    
    
#interleave the frame to make a thumbnail to export via stdout to labview
#for each color channel
# doonce = True
# for color in [0,1,2]:
#     col = frame[:,:,color][:,1::decimation] #interleave rows
#     col = col[1::decimation,:]#interleave columns
#     if doonce:
#         exportarray = np.empty([col.shape[0],col.shape[1],3],dtype=np.uint8)
#         doonce = False
#     exportarray[:,:,color] = col               
# #another thing to try, convert the PNG to hex and send that directly?
# import io
# import PIL
# from PIL import Image
# p = Image.fromarray(frame)

# f = io.BytesIO(frame)
# f.read(1000)
# byte_stream = io.BytesIO(frame)
# frames = iio.imwrite(byte_stream, frame, index=None)

# tmpFile = io.BytesIO()
# filename = '\\test.png'


#### png_encoded = iio.imwrite("<bytes>", frame, extension=".jpg")


# iio.imwrite(tmpFile+filename, frame)

# tmpFile.write(frame)
# f.read(1000)

# img_crop_pil = Image.fromarray(frame)
# img_crop_pil.save(tmpFile, format="JPG")


# savefig(tmpFile, format='png')
# filename = 'test.png'
# iio.imwrite(filename, frame)


        
    #close camera
    print('Terminating camera thread')