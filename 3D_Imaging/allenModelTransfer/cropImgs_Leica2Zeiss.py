'''
Crop images and convert to 16bits
Code by Saeed Alahmari , Jun 23rd, 2022

How to run the code
Type in the terminal the following command

python3 cropImgs.py -ImgsDir ./FoF1_220622_brightfield 

for help do

python3 cropImgs.py -h 


For questions contact Saeed at aalahmari.saeed@gmail.com
'''


import os 
import cv2
import numpy as np
import sys 
from tqdm import tqdm 

import argparse
def crop_and_convert2_16bit(path2Images,path2Save):
    if path2Save == '..':
        parentDir = os.path.dirname(path2Images)
        #print(parentDir)
        path2Save = os.path.join(parentDir,path2Images.rsplit('/',1)[-1]+'_croppedConverted2_16bits')
        if not os.path.exists(path2Save):
            os.makedirs(path2Save)
    print('Now cropping images from {} and converting to 16bits '.format(path2Images))
    print('Results are saved to {}'.format(path2Save))
    print('Please wait until this process complete ...')
    for image in tqdm(os.listdir(path2Images)):
        if image.startswith('.'):
            continue
        elif os.path.isdir(os.path.join(path2Images,image)):
            continue 
        elif not image.endswith('.tif'):
            continue
        #print(image)
        # Read Image ... 
        img = cv2.imread(os.path.join(path2Images,image),-1)
        # Crop the image, get first 1024 pixels from the x and y ....
        img = img[0:1024,0:1024,:]
        # convert image to uint16 bits .... 
        out_img = np.uint16(img)
        cv2.imwrite(os.path.join(path2Save,image),out_img)

def get_arguments():
   parser = argparse.ArgumentParser()
   parser.add_argument('-ImgsDir',help='directory to stack images')
   parser.add_argument('-saveDir',default='..',help='directory where images will be saved')
   args = parser.parse_args()
   crop_and_convert2_16bit(args.ImgsDir,args.saveDir)

if __name__ == "__main__":
    get_arguments()