# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 09:58:25 2025

@author: AndrewMiller
"""

def image_save(queues):
    from PIL import Image
    import os
    import config # https://docs.python.org/3/faq/programming.html#how-do-i-share-global-variables-across-modules
    savepath = r'C:\Users\AndrewMiller\Telescope\images'
    print('Starting image_save thread')
    
    while not config.end_program:
        size = queues['image_save_q'].qsize() #check that an image is in the queue
        if size >= 1:
            img_q = queues['image_save_q'].get(block=False) 
            filename = img_q['filename']
            img = img_q['img']
            print('saving ' + filename)
            #im = Image.fromarray(img)
            img.save(os.path.join(savepath, filename))
            print('done saving ' + filename)
    print('Terminating image_save thread')