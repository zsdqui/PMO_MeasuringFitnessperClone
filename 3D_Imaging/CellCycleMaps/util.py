# Util file
import numpy as np
from enum import unique
import os
import sys
import cv2
import sys
from tqdm import tqdm
from sklearn.preprocessing import LabelBinarizer
import pandas as pd
from sklearn.utils import shuffle
import json
import pickle
import tensorflow.keras as K
import tensorflow as tf
import tifffile


def normalize8(I):
    # print(I.dtype)
    # print(type(I))
    if isinstance(I, tf.Tensor):
        I = I.numpy()
    mn = I.min()
    mx = I.max()
    mx -= mn
    I = ((I - mn) / mx) * 255
    return I.astype(np.uint8)


def load_data(train_dir, dataframe):
    images = []
    labels = []
    print("Loading data...")
    for imagePath, label in tqdm(
        zip(dataframe["image"].tolist(), dataframe["label"].tolist())
    ):
        imagePath = imagePath.lstrip(os.sep)
        if imagePath.endswith(".png") or imagePath.endswith(".jpg"):
            image = cv2.imread(os.path.join(train_dir, imagePath), -1)
        else:
            image = tifffile.imread(os.path.join(train_dir, imagePath))
        # img_temp = image[0,:,:]
        # image = img_temp
        # print(imagePath,label)
        # print(image.shape)
        if image.shape[0] == 3:
            img_resized = np.zeros((3, 64, 64))
            img_resized[0, :, :] = cv2.resize(
                image[0, :, :], (64, 64), interpolation=cv2.INTER_CUBIC
            )
            img_resized[1, :, :] = cv2.resize(
                image[1, :, :], (64, 64), interpolation=cv2.INTER_CUBIC
            )
            img_resized[2, :, :] = cv2.resize(
                image[2, :, :], (64, 64), interpolation=cv2.INTER_CUBIC
            )

            img_transposed = img_resized.transpose(1, 2, 0)
            image = img_transposed
            # image = cv2.resize(image,(3,64,64),interpolation=cv2.INTER_CUBIC)
        else:
            image = cv2.resize(image, (64, 64), interpolation=cv2.INTER_CUBIC)

        if not image.dtype == np.uint8:
            image = normalize8(image)
        if image.shape[-1] == 3:
            pass
        else:
            image = np.stack([image, image, image], -1)
        images.append(image)
        labels.append(label)
    return images, labels


def Binarize_labels(label):
    label = np.array(label)
    mlb = LabelBinarizer()
    new_labels = mlb.fit_transform(label)
    return new_labels, mlb.classes_, mlb


def split_data_val_test(df):
    # Here the validation data is a two image with all correspondings z planes patches for each class
    # Here the validation data is a two image with all correspondings z planes patches for each class
    # Validation images:  FoF2_231005_fluorescent.nucleus  and FoF3001_220523_brightfield
    # Testing images: FoF2002004_221018_brightfield and FoF3_231005_fluorescent.nucleus
    print(df.shape)
    val_df = df(
        df["image"].str.startswith("FoF2_231005_fluorescent.nucleus")
        or df["image"].str.startswith("FoF3001_220523_brightfield")
    )
    print(val_df.shape)


def load_cellCycleData(path2Data):
    list_of_names = []
    list_of_labels = []
    df = pd.DataFrame()
    for folder in os.listdir(path2Data):
        if folder.startswith("._") or not os.path.isdir(
            os.path.join(path2Data, folder)
        ):
            continue
        for cellCycleImage in os.listdir(
            os.path.join(path2Data, folder, "images_padded")
        ):
            if not cellCycleImage.startswith("FoF"):
                continue
            imagePath = os.path.join(folder, "images_padded", cellCycleImage)
            image_label = int(cellCycleImage.split("cellCyle")[1].split("_")[0])
            list_of_names.append(imagePath)
            list_of_labels.append(image_label - 1)

    df["image"] = list_of_names
    df["label"] = list_of_labels
    return df


# Load the keras models
def load_keras_model(filePath):
    if filePath.endswith(".keras") and not filePath.startswith("._"):
        model = K.models.load_model(filePath)
        return model

    model_names = [
        file
        for file in os.listdir(filePath)
        if file.endswith(".keras") and not file.startswith("._")
    ]
    if len(model_names) > 1:
        print(
            "There should be a single model to use for fine-tuning, please specify the model name in -path2Model argument"
        )
        sys.exit()
    else:
        model = K.models.load_model(os.path.join(filePath, model_names[0]))
        return model


# path2Data = '/Volumes/WD Element/Collaboration/Moffitt_Noemi/BioinformaticsPaper/code/patches'
# data = load_cellCycleData(path2Data)
