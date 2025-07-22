# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 09:58:25 2025

@author: AndrewMiller
"""
import config # https://docs.python.org/3/faq/programming.html#how-do-i-share-global-variables-across-modules
import datetime
from PIL import Image #used for testing only


def camera_control(queues):

    #initialize camera
    print(config.end_program)
    
    print('starting camera thread')
    #open camera
    

    while not config.end_program:
            
        size = queues['camera_q'].qsize()
        if size >= 1:
            print('taking image')
            config.camera_ready = False
    
            #pull item from queue
            camera = queues['camera_q'].get(block=False)
            
            #unpack camera parameters (to be implemented)
            #image = take picture
            #fake image for testing
            path = r'C:\Users\AndrewMiller\Telescope\dishwasher-mustard.png'
            image =Image.open(path)
            #image = config.image
           

            #date and coordinates as filename
            
            t = datetime.datetime.now()
            filename = t.strftime('%Y-%m-%d %H_%M_%S.%f') + camera
            img = {'img':image, 'filename':filename+'.jpg'}
            
            queues['image_save_q'].put(img, block=False)   #put picture into image processor queue along with filename
    
            config.camera_ready = True
            print('camera done taking image')
        
    #close camera
    print('Terminating camera thread')