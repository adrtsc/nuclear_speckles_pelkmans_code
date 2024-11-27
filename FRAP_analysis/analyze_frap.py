#script to automate analysis of FRAP data from Leica .lif file (used for initial analyses and compared to manual analyses)
#need to add column to the output csv with Time information for later plotting 

from readlif.reader import LifFile
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import pandas as pd
from random import sample
from math import sqrt
from skimage.feature import blob_log
from skimage.measure import label, regionprops
from skimage.draw import disk 

import sys

liffi = LifFile('[fi path]') #provide path to .lif file from Leica
csv_out = 'automated_FRAP_analysis.csv' #provide path to save csv as output

def plot_roi(analysis_name,pre_bleach,post_bleach,diff,roi):
   fig, axes = plt.subplots(1, 3, figsize=(9, 3), sharex=True, sharey=True)
   ax = axes.ravel()
   ax[0].set_title('Pre-bleach')
   ax[0].imshow(pre_bleach)
   ax[0].set_axis_off()
   ax[1].set_title('Post-bleach 1')
   ax[1].imshow(post_bleach)
   ax[1].set_axis_off()
   ax[2].set_title('Laplacian of Gaussian')
   ax[2].imshow(diff)
   y,x,r = roi
   c = plt.Circle((x, y), r, color='red', linewidth=1, fill=False)
   ax[2].add_patch(c)
   ax[2].set_axis_off()
   plt.tight_layout()
   plt.savefig(analysis_name+'.pdf')

def plot_series(analysis_name,rois,img_set):
   fig, axes = plt.subplots(5,5,figsize=(10,10),sharex=True,sharey=True)

def get_max_intensity_roi(rois,img,dims):
   intensities = []
   for roi in rois:
      rr, cc = disk(roi[:2],roi[2],shape = dims)
      mask = np.zeros(dims)
      mask[rr, cc] = 1
      label_mask = label(mask)
      intens = measure_intensity(label_mask,img)
      intensities.append(intens[0])

   max_intensity = np.amax(intensities)
   roi1 = rois[np.where(intensities == np.amax(intensities))[0][0]]
   return(roi1,max_intensity)
 
def adjust_roi(roi_orig,roi0,img,pixel_shift,total_pixel_shift,dims):
   rois = []
   for i in range(-pixel_shift,pixel_shift+1): 
      for j in range(-pixel_shift,pixel_shift+1):
         new_roi = (roi0[0]+i,roi0[1]+j,roi0[2])
         if sum(abs(a-b) for a,b in zip(roi0[:2],new_roi[:2])) <= pixel_shift and \
           sum(abs(a-b) for a,b in zip(roi_orig[:2],new_roi[:2])) <= total_pixel_shift:
            rois.append(new_roi)
   return(get_max_intensity_roi(rois,img,dims))

def control_roi(roi_orig,img,dims):
   new_roi = (dims[0]*np.random.random(),dims[1]*np.random.random())
   rr,cc = disk(new_roi,roi_orig[2],shape = dims)
   mask = np.zeros(dims)
   mask[rr, cc] = 1
   label_mask = label(mask)
   intens = measure_intensity(label_mask,img)
   return(intens[0])

def measure_intensity(roi,img):
   roi_properties = regionprops(roi, img)[0]
   intensities = (roi_properties.mean_intensity,roi_properties.min_intensity,roi_properties.max_intensity)
   return(intensities)
 
img_names = {i.name:x for x,i in enumerate(liffi.get_iter_image())}
nfrap = 0
main_df = pd.DataFrame(columns = ['Measure', 'Mean', 'FRAP', 'condition'])

for img_name in img_names:
   print(img_name)
   if 'FRAP Pre Series' in img_name:
      nfrap+=1
      ind = int(img_name.split('FRAP Pre Series')[-1])
      pre = "{:02d}".format(ind)
      next = "{:02d}".format(ind+1)
      nextnext = "{:02d}".format(ind+2)

      pre_img = liffi.get_image(img_names[img_name])
      print(pre_img.path)
      img_set = [pre_img]
      post_img1_name = 'FRAP Pb1 Series'+next 
      post_img2_name = 'FRAP Pb2 Series'+nextnext
      try: 
         post_img1 = liffi.get_image(img_names[post_img1_name])
         img_set.append(post_img1)
      except KeyError:
         print('No images post-bleach for {}'.format(img_name))
         continue
      try:
         post_img2 = liffi.get_image(img_names[post_img2_name])
         img_set.append(post_img2)
      except KeyError:
         print('Single set of post-bleach images')
      last_pre_ind = pre_img.dims[3] - 1
      last_pre = np.array(pre_img.get_frame(z=0,t=last_pre_ind,c=0,m=0),dtype=np.int16)
      first_post = np.array(post_img1.get_frame(z=0,t=0,c=0,m=0),dtype=np.int16)
      diff = last_pre - first_post
      dims = (diff.shape[0], diff.shape[1])

      blobs_log = blob_log(diff, max_sigma=10, min_sigma = 2, num_sigma=40, threshold=.0001)
      blobs_log[:, 2] = blobs_log[:, 2] * sqrt(2)
      if len(blobs_log) > 0:
         roi_orig,detection_intensity = get_max_intensity_roi(blobs_log,diff,dims)
      else:
         roi_orig = (0,0,0)
         plot_roi('_'.join(img_name.split()),last_pre,first_post,diff,roi_orig)
         print('no blobs found for {}'.format(img_name))
         continue
      mask = np.zeros(dims)
      rr, cc = disk((blobs_log[0][0], blobs_log[0][1]), blobs_log[0][2],shape = dims)
      mask[rr, cc] = 1
      label_mask = label(mask)
      roi_orig = blobs_log[0]
      roi0 = tuple(blobs_log[0])

      #print(measure_intensity(label_mask,last_pre))      
      #print(measure_intensity(label_mask,first_post))

      intensities = []
      control_intensities = []
      try:
        for time in img_set:
         for t in range(time.dims[3]):
            frame = np.array(time.get_frame(z=0,t=t,c=0,m=0),dtype=np.int16)
            roi0,intensity = adjust_roi(roi_orig,roi0,frame,1,10, dims)
            control_intensities.append(control_roi(roi_orig,frame,dims)) 
            intensities.append(intensity)
      except IndexError:
         print('IndexError here')
         continue
      nmeasures = len(intensities)
      measures = range(1,nmeasures+4)
      frap = [nfrap]*(nmeasures+3)
      conditions = [img_name]*(nmeasures+3)      
      bg = [0]*nmeasures+[1]*3
      intensities = intensities + sample(sorted(control_intensities)[:int(nmeasures/5)],5)
      
      df = pd.DataFrame(list(zip(measures, intensities, frap, conditions,bg)),
               columns =['Measure', 'Mean', 'FRAP', 'condition','BG'])
      main_df = main_df.append(df)

main_df.to_csv(csv_out)
