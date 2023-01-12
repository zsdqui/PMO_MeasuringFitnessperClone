# Scripts to automate cellPose count 
# Saeed Alahmari

'''
Inputs:
1) path to an image. Image must be with .tif extension 
2) Type of images (cytoplasm, nucelus, mitochondria)


Outputs:
Two folders will be generated at the path to images folder.
the first folder, contains the coordinates of cell centers 
the second folder, contains the coordintates of entire segmeneted cells
'''
'''
HOW TO RUN THIS SCRIPT
First, make sure to create a conda enviroment and run the following command.
conda env create -f environment.yml
conda activate cellpose

Then run this command
python3 cellPose_getCount.py -dir <path2Image.tif> -t <typeofImages> -GPU <GPU number>

This code can be called as follows: python3 run_CellPose_3D_v2.py -imgPath mito.t.tif -t mitochondria -GPU 1

'''

import argparse
from cellpose import models
import cv2 
from cellpose.plot import mask_overlay
import numpy as np
import pandas as pd
import os
import sys 
from sys import platform
from moviepy.editor import ImageSequenceClip
from skimage import io, color
import skimage

np.random.seed(2020)


def visualize_mask(image,mask):
    colors = np.random.random((mask.max(),3))
    final_image = np.empty((image.shape[0],image.shape[1],image.shape[2],3),dtype=np.uint8)
    #print(final_image.dtype)
    #print(final_image.shape)
    for i in range(0,image.shape[0]):
        img = cv2.normalize(image[i,:,:], None, 0, 255, cv2.NORM_MINMAX, cv2.CV_8U)
        ColoredImg = color.label2rgb(mask[i,:,:],img,colors=colors,alpha=0.7, bg_label=0, bg_color=None)
        ColoredImg = cv2.normalize(ColoredImg, None, 0, 255, cv2.NORM_MINMAX, cv2.CV_8U)
        ColoredImg = ColoredImg.astype(np.uint8)
        #cv2.imwrite('image_overlay.tiff',ColoredImg)
        final_image[i,:,:,:] = ColoredImg
        #cv2.imwrite('image_overlay_first.tiff',final_image[0,:,:,:])
    return final_image 

def get_blob_prop(msk):
  contours,hierarchy = cv2.findContours(msk, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)
  for cnt in contours:
    try:
      m = cv2.moments(cnt)
      x = m['m10']/m['m00']
      y = m['m01'] /m['m00']
      area = cv2.contourArea(cnt)
      perimeter = cv2.arcLength(cnt,True)
    except:
      print('Trying to find moments of image')
      continue  
    return {'centroid':(x,y),'area':area,'perimeter':perimeter}

def get_fileName(image_path):
    if platform == 'linux' or platform == 'linux2':
        print('platform is linux')
        basepath,file_name = image_path.rsplit('/',1)[0],image_path.rsplit('/',1)[1]
    elif platform == 'darwin':
        print('platform is OS X')
        basepath,file_name = image_path.rsplit('/',1)[0], image_path.rsplit('/',1)[1]
    elif platform == 'win32':
        print('Windows')
        basepath,file_name = image_path.rsplit('\\',1)[0], image_path.rsplit('\\',1)[1]
    return basepath,file_name

#get count in pandas dataframe, to be saved to csv
def get_count2csv(list_of_cells_props):
    df = pd.DataFrame()
    df['centroid'] = [i['centroid'] for i in list_of_cells_props]
    df['area'] = [i['area'] for i in list_of_cells_props]
    df['perimeter'] = [i['perimeter'] for i in list_of_cells_props]
    return df

# get detections of masks overlay
def get_annotated_img(img,mask):
    overlay = mask_overlay(img,mask)
    return overlay
def get_center(maski):
  maski = maski.astype(np.uint8)
  M = cv2.moments(maski)
  try:
    cX = int(M["m10"] / M["m00"])
    cY = int(M["m01"] / M["m00"])
  except:
    cX = 0 
    cY = 0
  return (cX,cY)

# Save the cell centers in a file named with the mask name in the folder (Cells_center_coordinates)
def create_csv(images_path,locations_list,f):
  x = [i[0] for i in locations_list]
  y = [i[1] for i in locations_list]
  z = [i[2] for i in locations_list]

  df = pd.DataFrame()
  df['x'] = x
  df['y'] = y
  df['z'] = z

  if not os.path.isdir(os.path.join(images_path,'Cells_center_coordinates')):
    os.makedirs(os.path.join(images_path,'Cells_center_coordinates'))
  df.to_csv(os.path.join(images_path,'Cells_center_coordinates',f+'_Cells_Centers.csv'),index=False)

def euclidean_distance(p1,p2):
  distance = np.sqrt(np.power(p1[0] - p2[0],2) + np.power(p1[1] - p2[1],2))
  return distance
#Get the coordinates for all the cells. 
#Each cell coordinates will be written to a file with the mask name followed by the cell number in the folder (All_Cells_coordinates)
def get_all_coordinates(img,slice):
  indices = np.argwhere(img > 0)
  xyz_loc = []
  if len(indices > 0):
    indices=indices.tolist()
    xyz_loc = [j+[slice] for j in indices]
  return xyz_loc

def get_coordinates_of_cells(path2Images,mask,f):
        #Get center for each cell 
    print('Getting the coordinates of the image ....')
    if not os.path.isdir(os.path.join(path2Images,'All_Cells_coordinates')):
        os.makedirs(os.path.join(path2Images,'All_Cells_coordinates'))
    All_images_locations = []
    List_of_cell_Locations = []
    List_of_cell_Locations2 = []

    for i in range(1,mask.max()):
        maski = (mask == i)*255
        List_of_Images = [maski[i,:,:] for i in range(0,mask.shape[0])]
        list_of_centers = []
        xyz_loc_accum = []
        for ind,img in enumerate(List_of_Images):
            xy = get_center(img)
            xyz_loc = get_all_coordinates(img,ind)  #slice index (ind)
            if len(xyz_loc) > 0:
                if len(xyz_loc_accum) == 0:
                    xyz_loc_accum = xyz_loc 
                else:
                    xyz_loc_accum = xyz_loc_accum + xyz_loc 
            
            if xy[0] ==0 and xy[1] == 0:
                continue
            else:
                list_of_centers.append(xy)
        #print(len(list_of_centers))
        #Save all the location of a cell
        df = pd.DataFrame(xyz_loc_accum,columns=['y','x','z'])
        f_output = f.split('/')[-1].split('.tif')[0]
        if not os.path.exists(os.path.join(path2Images,'All_Cells_coordinates',f_output)):
            os.makedirs(os.path.join(path2Images,'All_Cells_coordinates',f_output))
        df.to_csv(os.path.join(path2Images,'All_Cells_coordinates',f_output,f_output+'_cell_'+str(i)+'_coordinates.csv'),index=False)

        center_index = np.floor(len(list_of_centers)/2)
        #print(center_index)
        center = list_of_centers[int(center_index)]
        center = tuple([center[0],center[1],int(center_index)])
        #  #print(center)
        List_of_cell_Locations.append(center)
    return List_of_cell_Locations



#call the model, and run the segmentations
def cellPose_getCount(images_path,images_type, machine_type='CPU'):

    # model_type='cyto' or model_type='nuclei'
    if machine_type.upper() == 'GPU': 
        print('Using GPU..')
        model = models.Cellpose(gpu=True, model_type='cyto')
    elif machine_type.upper() == 'CPU':
        print('Using CPU..')
        model = models.Cellpose(gpu=False, model_type='cyto')
    else:
        print('Wrong machine_type = {}'.format(machine_type))
        sys.exit()
        


    if not os.path.exists(images_path):
        print('Wrong image path = {}'.format(images_path))
        sys.exit()
    else:
        basepath,imageName = get_fileName(images_path) 
        print('basepath {}'.format(basepath))
        print('imageName {}'.format(imageName))
        if imageName.startswith('.'):
            print('Image {} must not start with dot (.)'.format(imageName)) 
            sys.exit()
        if not imageName.endswith('.tif'):
            print('Image {} must be .tif file'.format(imageName)) 
            sys.exit()
        print('Processing {} ...'.format(imageName))
        #imageName = get_fileName(image_path)
        #img = cv2.imread(os.path.join(images_path,imageName),-1)  # reading an image or stack of images 
        img = skimage.io.imread(images_path)
        print(img.shape)
        #img = img[:15,:,:]
        #model eval in 3D:
        #Option 1: use tif images as input. 
        # mask, flow, style, diam = model.eval(img, diameter=None, channels=[0,0], do_3D=True)
        #Option 2
        if images_type == 'cytoplasm':
            masks, flows, styles, diams = model.eval(img, diameter=30, channels=[0,0], do_3D=True)
        elif images_type == 'mitochondria':
            masks, flows, styles, diams = model.eval(img, diameter=15, channels=[0,0], do_3D=True, anisotropy=14)#  Run the code for 3D segmentation
        elif images_type == 'nucleus':
            masks, flows, styles, diams = model.eval(img, diameter=30, channels=[0,0], stitch_threshold=0.95)
        else:
            print('Sorry, was not able to run the code due to error on the image type')
    
        # Get coordinates of cells
        list_of_cell_centers = get_coordinates_of_cells(basepath,masks,imageName)

        #Save the centers coordinates
        print('Writing the coordinates for each image ....')
        f_output = imageName.split('.tif')[0]
        create_csv(basepath,list_of_cell_centers,f_output)

        final_image = visualize_mask(img,masks)
        list_of_images = [final_image[i,:,:,:] for i in range(0,final_image.shape[0])]
        clip = ImageSequenceClip(list(list_of_images), fps=5)
        clip.write_gif(os.path.join(basepath,f_output+'_visualized.gif'), fps=5)

def main(parser):
    args = parser.parse_args()
    image_type = args.t.lower()
    if image_type.startswith('cyto'):
        image_type = 'cytoplasm'
    elif image_type.startswith('nuc'):
        image_type = 'nucleus'
    elif image_type.startswith('mito'):
        image_type = 'mitochondria'
    else:
        print('-t input should be either cytoplasm, mitochondria, or nucleus')
        sys.exit()

    #Visible devices with GPU number 
    os.environ["CUDA_VISIBLE_DEVICES"] = str(args.GPU)
    cellPose_getCount(args.imgPath,image_type,machine_type='GPU')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Getting argument for running cellpose")
    parser.add_argument('-imgPath',required=True, help="Path to the an image (.tif) to run cellpose")
    parser.add_argument('-t',required=True,help="Type of images to run the cellpose for. The code support only (cytoplasm, mitochondria, or nucleus)")
    parser.add_argument('-GPU',required=False,help="GPU number to run the code on",default=0)
    main(parser)




