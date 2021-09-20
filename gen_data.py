import numpy as np
import pandas as pd
import os
import pickle
from sklearn import preprocessing  



def unpickle(file):
    with open(file, 'rb') as fo:
        dict = pickle.load(fo)
    return dict

def load_databatch(data_folder, idx, img_size=32):
    data_file = os.path.join(data_folder, 'train_data_batch_')

    d = unpickle(data_file + str(idx))
    x = d['data']

    return x

def gen_imagedata():
    color_his = []
    for k in range(1,11):

        # change filepath to your image file folder
        X = load_databatch("E:\\dataset_LIMS\\Imagenet32_train", k)
        lenth = len(X[0])
        for i in range(len(X)):
            hist, bins = np.histogram(X[i], bins=32)
            color_his.append(hist/lenth)
    data = np.array(color_his)

    print(data.dtype)
    print(data.shape)
    min_max_scaler = preprocessing.MinMaxScaler()
    X_scaled = min_max_scaler.fit_transform(data) 
    print(X_scaled.dtype)
    print(X_scaled.shape)

    # change output file path to your path
    np.savetxt('E:\\dataset_LIMS\\Imagenet32_train\\color_32.txt',X_scaled,fmt="%.5f",delimiter=',', newline='\n', header='', footer='')
    np.save('E:\\dataset_LIMS\\Imagenet32_train\\color_32.npy',X_scaled)


def gen_forestdata():
    # change input file path to your forest cover type data storage path
    df = pd.read_csv("E:\\dataset_LISHD\\forest\\test.csv", error_bad_lines=False)
    print(df.head())
    print(df.shape)
    print(df.isna().sum())


    data = df[['Elevation','Aspect','Slope']].copy()

    print(data.head())
    print(data.shape)

    data['dis'] = df['Horizontal_Distance_To_Hydrology'].apply(lambda x: x * x) + df['Vertical_Distance_To_Hydrology'].apply(lambda x: x * x)
    data['dis'] = data['dis'].apply(lambda x: np.sqrt(x))
    data['Horizontal_Distance_To_Roadways'] = df['Horizontal_Distance_To_Roadways']
    data['time'] = (df['Hillshade_9am'] + df['Hillshade_Noon'] + df['Hillshade_3pm'])/3
    print(data.head())
    print(data.shape)
    print(data.isna().sum())

    print(data['Elevation'].max(axis=0))
    print(data['Aspect'].max(axis=0))
    print(data['Slope'].max(axis=0))
    print(data['dis'].max(axis=0))
    print(data['Horizontal_Distance_To_Roadways'].max(axis=0))
    print(data['time'].max(axis=0))

    data = np.array(data,dtype=float)

    max_val = 0.0
    for i in range(6):
        if data[:, i].max(axis=0) > max_val:
            max_val = data[:, i].max(axis=0)

    print(max_val)

    for i in range(6):
        data[:, i] = data[:, i] / max_val
        noise = np.random.randint(0,1000,size=[len(data),1])
        noise = noise.astype(float)
        noise[:, 0] = noise[:, 0] / 100000000
        data[:, i] = data[:, i] + noise[:, 0]

    # change output file path to your path
    np.savetxt('E:\\dataset_LISHD\\forest\\forest.txt',data,fmt="%.8f",delimiter=',', newline='\n', header='', footer='')
    np.save('E:\\dataset_LISHD\\forest\\forest.npy',data)

if __name__ == '__main__':
    gen_imagedata()
    gen_forestdata()