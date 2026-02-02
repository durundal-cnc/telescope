#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 26 10:01:32 2026

@author: andrewmiller
"""

#doesn't work in spyder for unknown reasons (kernel just dies), maybe because not authorized by macos in security panel? Had to run in terminal to get the allow prompt

import cv2

frameWidth = 640
frameHeight = 480
cap = cv2.VideoCapture(0)
#when on Dell monitor dock:
#0 is Logi 4k
#1 is internal
#2 is iphone

#when on laptop only
#0 is iPhone camera
#1 is internal webcam
cap.set(3, frameWidth)
cap.set(4, frameHeight)
cap.set(10,150)

while True:
   success, img = cap.read()
   cv2.imshow("Result", img)
   if cv2.waitKey(1) & 0xFF == ord('q'):
       break



#%%
# from astropy.coordinates import SkyCoord
# from astropy import units as u
# ra1 = 10
# dec1 = 10
# ra2 = 10
# dec2 = 10
# distance1 = 1
# c = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree)  
# catalog = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)  
# idx, d2d, d3d = c.match_to_catalog_sky(catalog)  

# c = np.array([SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree)])

# idxc, idxcatalog, d2d, d3d = catalog.search_around_sky(c, 1*u.deg)  


# np.all(c[idxc].separation(catalog[idxcatalog]) == d2d)  
# np.all(c[idxc].separation_3d(catalog[idxcatalog]) == d3d)  
# print catalog_objectnames[idxcatalog]  
# ['NGC 1234' 'NGC 4567' ...]

# import numpy as np
# from astropy.coordinates import SkyCoord, match_coordinates_sky
# from astropy import units as u

# scatalog = SkyCoord(ra=np.linspace(0, 1, 10)*u.degree, dec=np.ones(10)*u.degree)
# pcatalog = SkyCoord(ra=np.linspace(0, 1, 100)*u.degree, dec=np.ones(100)*u.degree)
# idx, d2d, d3d = match_coordinates_sky(scatalog, pcatalog, nthneighbor=1)
# print(idx.shape)
# print(d2d.shape)
# #%%
# def plot_nearby_stars(minimum_brightness_plot = 13, minimum_brightness_annotation = 11):
#     from astropy.table import QTable
#     from astroquery.gaia import Gaia
#     Gaia.ROW_LIMIT = 10000  # Set the row limit for returned data
    
#     import numpy as np
#     from astropy.coordinates import SkyCoord, match_coordinates_sky
#     from astropy import units as u
    
#     from astroquery.simbad import Simbad
    
    
#     #take earth location
#     loc = config.my_locs['EarthLocation']
#     current_time = Time(datetime.now(timezone.utc)) #have to make sure datetime is in utc for all the astro tools unless you specify it in them indivudally
#     #take current az/el at that location
#     altaz_frame = AltAz(obstime=current_time, location=loc)
#     sc = SkyCoord(alt=config.az_angle_PV*u.deg, az=config.el_angle_PV*u.deg, frame=altaz_frame) #need to put in the iPhone compass direction offsets
#     icrs_coo = sc.transform_to('icrs')
#     print(f"RA/Dec: {icrs_coo.ra.deg}, {icrs_coo.dec.deg}")
    
#     #import GAIA database of objects
#     #TODO: figure out how to store star catalog locally
#     #get list of objects within radius of that skycoord 
#     job = Gaia.cone_search_async(sc, radius=0.5 * u.deg)
#     ngc188_table = job.get_results()

#     # only keep stars brighter than G=19 magnitude
#     ngc188_table = ngc188_table[ngc188_table["phot_g_mean_mag"] < minimum_brightness_plot * u.mag] #13 is limit of ~4" telescope
    
#     # add all ids in the SIMBAD results
#     Simbad.add_votable_fields('ids')
#     named = Simbad.query_objects(ngc188_table['designation'])
    
    
#     #for readability
#     import pandas as pd
#     p = ngc188_table.to_pandas()      #convert to Pandas dataframe
#     q = named.to_pandas()
#     #plot according to brightness
#     #gaia_dist = Distance(parallax=ngc188_table_3d["parallax"].filled(np.nan)) #unused
#     gaia_magnitude = ngc188_table["phot_g_mean_mag"].filled(np.nan)
    
    
    
#     #plot and label so can match what seeing through eyepiece to what's on sky there
#     import matplotlib.pyplot as plt
#     fig, ax = plt.subplots(figsize=(6.5, 5.2), constrained_layout=True)
#     cs = ax.scatter(
#         ngc188_table['ra'],
#         ngc188_table['dec'],
#         c=gaia_magnitude,
#         s=5,
#         vmin=min(gaia_magnitude),
#         vmax=max(gaia_magnitude),
#         cmap="gray",
#     )
#     cb = fig.colorbar(cs)
#     cb.set_label(f"magnitude")
    
#     #restrict name plotting to brighter stars
#     for i, txt in enumerate(named['main_id']):
#         if ngc188_table['phot_g_mean_mag'][i] < minimum_brightness_annotation:
#             ax.annotate(txt, (ngc188_table['ra'][i], ngc188_table['dec'][i]))
    
#     ax.set_xlabel("RA [deg]")
#     ax.set_ylabel("Dec [deg]")
    
#     ax.set_title("Gaia DR2 sources near NGC 188", fontsize=18)
#     return plt.gcf() 
