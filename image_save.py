# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 09:58:25 2025

@author: AndrewMiller
"""
import os
import config # https://docs.python.org/3/faq/programming.html#how-do-i-share-global-variables-across-modules
import cv2 
import time

def image_save(queues):
    print('Starting image_save thread')
    
    if os.name == 'nt': #Windows 11
        savepath = r'C:\users\andym\telescope\images'
    elif os.name == 'posix':
        savepath = r'/Users/andrewmiller/telescope/images'
    
    while not config.end_program:
        config.image_process_ready = True #signal loop is ready
        
        if len(config.photos) > 0:
            frame, filename = config.photos.pop()
    
            print('JPEG start'+ str(round(time.time() * 1000)))
            is_success, jpg_encoded = cv2.imencode(".jpg", frame)
            print('JPEG end'+ str(round(time.time() * 1000)))
            
            print('JPEG print start'+ str(round(time.time() * 1000)))
            print(bytes(jpg_encoded).hex()) #need to add newline?
            print('JPEG print end '+ str(round(time.time() * 1000)))

            print('imwrite start' + str(round(time.time() * 1000)))
            cv2.imwrite(os.path.join(savepath, filename+'.jpg'), frame) 
            print('imwrite end ' + str(round(time.time() * 1000)))
        
    print('Terminating image_save thread')