
from random import shuffle
import tensorflow as tf 
import tensorflow.keras as k
from tensorflow.keras.preprocessing import image 
import os 
import tensorflow.keras as K
import sys
import numpy as np
import pandas as pd 
from util import load_data,Binarize_labels,load_cellCycleData,load_keras_model
from sklearn.metrics import accuracy_score
import random as rn
import argparse
from sklearn.model_selection import train_test_split
from DataGenerator import DataGenerator 
from custom_resnet import create_resnet
from custom_cnn import buildModel
from tensorflow.keras.losses import categorical_crossentropy
from pre_augment import AugmentImages


import tensorflow.keras as k
k.backend.set_image_data_format('channels_last')

'''
Note: Due to this error  GPU MaxPool gradient ops do not yet have a deterministic XLA implementation. Error occurred when finalizing GeneratorDataset iterator: 
FAILED_PRECONDITION: Python interpreter state is not initialized. The process may be terminated. I disabled the deterministic behaviour of tensorflow
'''
def init_seeds(seed=2022):
    np.random.seed(seed)
    rn.seed(seed)
    session_conf = tf.compat.v1.ConfigProto()
    session_conf = tf.compat.v1.ConfigProto(intra_op_parallelism_threads=1,
                              inter_op_parallelism_threads=1)
    #session_conf.gpu_options.visible_device_list = gpu
    os.environ['TF_CUDNN_DETERMINISTIC'] ='true'
    os.environ['TF_DETERMINISTIC_OPS'] = 'true'
    import tensorflow.keras as k
    k.backend.set_image_data_format('channels_last')
    tf.random.set_seed(seed)
    sess = tf.compat.v1.Session(graph=tf.compat.v1.get_default_graph(),config=session_conf)
    #tf.compat.v1.keras.backend.set_session(sess)
    return sess

def ResNet(model,input_shape=(64,64,3)): 
  input = k.layers.Input(input_shape)
  layer = model(input)
  layer = k.layers.BatchNormalization()(layer)
  layer = tf.keras.layers.Activation('relu')(layer)
  layer = k.layers.GlobalAveragePooling2D()(layer)
  layer = k.layers.Dense(units=512, activation='relu')(layer)
  layer = k.layers.Dense(units=4,activation='softmax')(layer)
  DenseNet_model = k.models.Model(inputs=input,outputs=layer)
  #print(DenseNet_model.summary()) 
  DenseNet_model.compile(loss=categorical_crossentropy,optimizer=k.optimizers.Adam(),metrics=['accuracy'])
  return DenseNet_model

def train(train_dir,val_dir,nlabels,nepochs,expName,arch,batch_size,path2PretrainedModel=None,pre_augment=False):
    try:
        sess.close()
        k.clear_session()
    except:
        pass
    #sess = init_seeds(seed=2022)

    if not os.path.exists('models'):
        os.makedirs('models')

    
    train_df = load_cellCycleData(train_dir)
    if not os.path.exists(os.path.join(train_dir,'train_data.csv')) and not os.path.exists(os.path.join(train_dir,'val_data.csv')) and not os.path.exists(os.path.join(train_dir,'test_data.csv')):
        print('Didnt find csv files, splitting the data...')
        train_df, temp_df = train_test_split(train_df, random_state=2024,test_size=0.4, stratify=train_df['label'])
        val_df, test_df = train_test_split(temp_df, random_state = 2024, test_size=0.5, stratify=temp_df['label'])
        train_df.to_csv(train_dir+'/train_data.csv',index=False)
        val_df.to_csv(train_dir+'/val_data.csv',index=False)
        test_df.to_csv(train_dir+'/test_data.csv',index=False)
    else:
        train_df = pd.read_csv(train_dir+'/train_data.csv')
        val_df = pd.read_csv(train_dir+'/val_data.csv')
    #Pre-augment option goes here and updating the df_train file
    if pre_augment == True:
        if train_dir.endswith('/'): # remove the most right '/' in the path to images.
            train_dir = train_dir.rsplit('/',1)[0]
        if not os.path.exists(train_dir+'-aug'):
            print('Starting pre-augmentation ...')
            Augmentor = AugmentImages(train_dir,train_dir+'-aug')
            Augmentor.iterate_df()
        train_df = load_cellCycleData(train_dir+'-aug')
        training_generator = DataGenerator(train_df,num_classes = 4, batch_size=batch_size,augment=False)
        val_generator = DataGenerator(val_df,num_classes = 4, batch_size=batch_size,augment=False)
    else:
        #val_df = pd.read_csv('./'+train_dir+'/val_data.csv')
        training_generator = DataGenerator(train_df,num_classes = 4, batch_size=batch_size,augment=True)
        val_generator = DataGenerator(val_df,num_classes = 4, batch_size=batch_size,augment=True)

    model_checkpoint = tf.keras.callbacks.ModelCheckpoint(os.path.join('models',expName+'_'+arch+'.keras'),monitor="val_loss",verbose=0,
        save_best_only=True,
        save_weights_only=False,
        mode="auto",
        save_freq="epoch")
    
    #model = ResNet(ResNet_model,(64,64,3))
    if arch == 'custom_cnn':
        # hyperparameters
        drop = 0.3 # dropout rate
        hidden_layers = 2 # number of hidden convolutional/max pooling layer sets
        batch = 32 # batch size
        neurons = 28 # number of convolution kernels
        model = buildModel(neurons,drop,hidden_layers,nb_classes=4)
        if path2PretrainedModel is not None: # fine-tune model path available.
            # load the model from the path. 
            print('Loading pre-trained model, please wait...')
            model = load_keras_model(path2PretrainedModel)
    elif arch == 'custom_resnet':
        model = create_resnet()
        if path2PretrainedModel is not None: # fine-tune model path available.
            # load the model from the path. 
            print('Loading pre-trained model, please wait...')
            model = load_keras_model(path2PretrainedModel)
    model.summary()
    print('Fitting the model with images of one label')
    model.fit(training_generator,
        validation_data= val_generator,
        batch_size=batch_size,steps_per_epoch= int(training_generator.__len__()/batch_size), epochs=nepochs, callbacks= [model_checkpoint])

def get_test_accuracy_per_label(df):
    df_results_summary = pd.DataFrame()
    accuracy_list = []
    unique_labels = set(df['true_labels'].tolist())
    count_list = []
    for l in unique_labels:
        df_subset = df[df['true_labels'] == l]
        accuracy = accuracy_score(df_subset['true_labels'].tolist(),df_subset['pred'].tolist())*100
        print('Accuracy for label {} for {} instance is {}'.format(l,str(df_subset.shape[0]),accuracy))
        accuracy_list.append(accuracy)
        count_list.append(df_subset.shape[0])

    df_results_summary['true_label'] = list(unique_labels)
    df_results_summary['Accuracy'] = accuracy_list
    df_results_summary['Count'] = count_list
    return df_results_summary

        
def test(train_dir, test_dir,exp,arch):
    #test_datagen = image.ImageDataGenerator(rescale=1./255)
    df = pd.read_csv(os.path.join(test_dir + 'test_data.csv'))
    img_list,labels = load_data(train_dir,df)
    #test_datagen.fit(img_list_test)
    print(len(img_list))
    #img_test = np.concatenate(img_list, axis=0)
    img_test = np.array(img_list)
    print(img_test.shape)
    img_test = img_test / 255.

    #CellCycle_true,_ = Binarize_labels(labels)
    CellCycle_true = K.utils.to_categorical(labels,4)

    model = tf.keras.models.load_model(
        os.path.join('models',exp+'_'+arch+'.keras'), custom_objects=None, compile=True)

    bottle_model =  k.models.Model(inputs=model.input, outputs=model.layers[-3].output)
    print(f"Bottleneck layer (before dropout) is: {model.layers[-3].name}\n")
    
    #save bottleneck features to dataframe
    feats = bottle_model.predict(img_test)
    df_feats = pd.DataFrame(feats)

    #results = model.evaluate(test_datagen.flow(img_list_test,y=s_label_numeric_test, batch_size = 32, shuffle=False))
    cellCycle = model.predict(img_test)
    print(cellCycle.shape)
 

    pred = np.argmax(cellCycle,axis=1)
    labels_true = np.argmax(CellCycle_true,axis=1)
    print('state_pred {}'.format(pred.shape))
    print('state_true {}'.format(labels_true.shape))

    df['pred'] = pred
    df['true_labels'] = labels_true
    #print bottleneck features to csv with labels
    combined_df = pd.concat([df['true_labels'],df['pred'],df_feats],axis=1)
    combined_df.to_csv(os.path.join('./Results', exp + '_bottleneck_features.csv'), index=False)

    accuracy = accuracy_score(labels_true,pred)
    if not os.path.exists('./Results'):
        os.makedirs('./Results')
    with open(os.path.join('Results',exp+'.txt'),"w") as txt:
        txt.write('accuracy for {} is \n {}\n\n'.format(exp,accuracy))
    print('Accuracy is {}'.format(accuracy))
        
    df_results = get_test_accuracy_per_label(df)
    df_results.to_csv(os.path.join('Results',exp+'_test.txt'), header=True, index=None, sep=' ', mode='a')
    
def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-dataDir',help='Train directory',required=False)
    parser.add_argument('-valDir',help='Val directory',required=False)
    parser.add_argument('-epochs',help='training epochs', required=False)
    parser.add_argument('-arch',help='name of architecture to use, custom_resnet is for the method implement in the paper titled'
                        '(An Imbalanced Image Classification Method for the Cell Cycle Phase, whereas the option custom_cnn is for the paper'
                        'Robust classification of cell cycle phase and biological feature extraction by image-based deep learning',required=True)
    parser.add_argument('-exp',help='experiment name',required=True)
    parser.add_argument('-m',help='mode either train or test',required = True)
    parser.add_argument('-nlabel',help='number of label for the data ex. 1 or 2', required= False)
    parser.add_argument('-testDir',help='Test directory containint test.csv',required=False,default='')
    parser.add_argument('-batch_size',help='Batch size',required=False,default=32)
    parser.add_argument('--finetune',help='Finetune a pre-trained model for epochs',action='store_true')
    parser.add_argument('-path2Model',help='path to model to fine-tune',required=False,default=None)
    parser.add_argument('-finetune_epochs',help='Number of epochs to fine-tune',required=False,default=100)
    parser.add_argument('--pre_augment',help='Pre-augment images before start training, if False, online-augmentation per batch is done',action='store_true')
    args = parser.parse_args()
    if args.m == 'train':
        if args.finetune and args.path2Model is None:
            print('Please specify the path to pre-trained model for finetuning')
            sys.exit()
        if args.epochs is None:
            print('Please specify argments')
            sys.exit()
        elif args.finetune and args.path2Model is not None:
            print('Started finetuning of {} model from path {}'.format(args.arch,args.path2Model))
            train(args.dataDir,args.valDir,args.nlabel,int(args.finetune_epochs),args.exp+'_fine-tuned',args.arch,int(args.batch_size),args.path2Model,args.pre_augment)
        elif not args.finetune:
            print('Started training  {} model '.format(args.exp))
            train(args.dataDir,args.valDir,args.nlabel,int(args.epochs),args.exp,args.arch,int(args.batch_size),None,args.pre_augment)
    else:
        if args.testDir == '':
            print('Error, please enter test directory')
            sys.exit()
        else:
            print('Started testing {} model'.format(args.exp))
            if args.finetune:
                test(args.dataDir, args.testDir, args.exp+'_fine-tuned',args.arch) #testing finetuned model
            else:
                test(args.dataDir, args.testDir,args.exp,args.arch) #testing 


if __name__ =="__main__":
    main()
