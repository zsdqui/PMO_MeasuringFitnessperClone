
#Script to resize an image to 856x856px then pad the image with zeros (84 pixels in each side)
#Created by Saeed Alahmari 
# Jan 8th, 2023 


'''
How to run the code:
1) open the terminal and type: python3 preprocess.py -image <path2Image>
2) The results will be written to a folder at the same directory of the code, name of the directory is resized_cropped
3) If you want to specify a directory for the output, simply use -saveToDir <path2SaveTheOutputImage> 
'''
import cv2
import numpy as np 
from tqdm import tqdm 
import sys
import os 
import argparse 


def resize_pad(path2Image,path2SaveImage):
    imageName = path2Image.rsplit(os.sep,1)[-1]
    print(imageName)
    dsize = (856,856)
    padding_width = ((84, 84),(84, 84),(0,0))
    if not imageName.endswith('.tif'):
        print('Error, input image has to be a .tif image')
        print('Your input is {}'.format(path2Image))
        sys.exit()
    print('Reading an image')
    image = cv2.imread(path2Image,-1)
    print('resizing and padding the input image')
    new_image = cv2.resize(image,dsize,interpolation=cv2.INTER_AREA)
    if len(new_image.shape) == 3 and new_image.shape[-1] == 3:
        image_padded = np.pad(new_image, padding_width, 'constant', constant_values=(0, 0))
    elif len(new_image.shape) == 2:
        image_padded = np.pad(new_image, 84, 'constant', constant_values=(0, 0))
    else:
        print('Error, unsupported image, please check with the developer')
        print('Error description: image shape check')
    if not os.path.exists(path2SaveImage):
        os.makedirs(path2SaveImage)
    cv2.imwrite(os.path.join(path2SaveImage,path2Image),image_padded)
    print('Resized and Padded image is saved to {}'.format(path2SaveImage))


if __name__ =="__main__":
    print('In main')
    parser = argparse.ArgumentParser()
    parser.add_argument('-image',help='Path to an image',required=True)
    parser.add_argument('-saveToDir',help='Directory path where the resized and padded images is saved',required=False,default='./resized_cropped')
    args = parser.parse_args()
    args = parser.parse_args()
    resize_pad(args.image,args.saveToDir)

