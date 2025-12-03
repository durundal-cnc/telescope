# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 13:17:40 2025

@author: AndrewMiller
"""

import cv2
import math
from matplotlib import pyplot as plt
import os
from PIL import Image
import numpy as np

#acquire photo stack & focus location

#back out to start point

#compute optimal focus point

#drive to there


SHOW_STEP_IMAGES = True
SHOW_STEP_INFO = True
focus_scores: list[float] = []

def calculate_focus_score(image, blur, position):
    image_filtered = cv2.medianBlur(image, blur)
    laplacian = cv2.Laplacian(image_filtered, cv2.CV_64F)
    focus_score = laplacian.var()
   
    if SHOW_STEP_IMAGES:
        focus_scores.append(focus_score)
        grayscale_laplacian = cv2.convertScaleAbs(laplacian, alpha=50)
        fig = plt.figure(figsize=(10, 2))
        
        ax1 = fig.add_subplot(131)
        ax1.imshow(image_filtered)
        ax1.set_title("Filtered")
        ax1.set_xticks([])
        ax1.set_yticks([])
        
        ax2 = fig.add_subplot(132)
        ax2.imshow(grayscale_laplacian)
        ax2.set_title("Laplacian")
        ax2.set_xticks([])
        ax2.set_yticks([])
        
        ax3 = fig.add_subplot(133)
        ax3.plot(range(len(focus_scores)), focus_scores, color="red")
        ax3.set_title("Variance")
        ax3.set_xticks([])
        ax3.set_yticks([])
        
        #plt.savefig(figure_file_name(position))
        plt.close()
   
    return focus_score

best_focus_score = 0
best_focus_score = 0
best_focus_position = 0
end_mm = 10
start_mm = 0
step_size_mm = .1
blur = 1

file_loc = r'/Users/andrewmiller/telescope/001'
files = [ name for name in os.listdir(file_loc) if name[-4:] == '.jpg' ]

# How many steps to take to achieve the desired step size, +1 to check end_mm
steps = math.ceil((end_mm - start_mm)/step_size_mm) + 1
#for step in range(0, steps):
position = 0
for file in files:
    position = position + 1
    #position = min(start_mm + step * step_size_mm, end_mm)
    #z_axis.move_absolute(position, Units.LENGTH_MILLIMETRES)
    #image = get_image(cam)
    image = Image.open(os.path.join(file_loc, file)) #load file]
    image = np.array(image)

    focus_score = calculate_focus_score(image, blur, position)
    if focus_score > best_focus_score:
        best_focus_position = position
        best_focus_score = focus_score