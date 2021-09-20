# LIMS

LIMS is a novel index structure that integrates compact-partitioning methods, pivot-based techniques and learned indexes to support efficient similarity query processing in metric spaces. More details about LIMS can be found in our paper.

## Getting Started

### 1. Prepare datasets

ImageNet : https://image-net.org/download-images.php

Forest Cover Type : https://www.kaggle.com/c/forest-cover-type-prediction/data

```bash
python gen_data.py
```

### 2. Change path and parameter

change parameter for LIMS in ./common/Constants.cpp to fit your own dataset

change dataset path and query file path in main.cpp 

### 3. Run

we use [FFT](https://github.com/ZJU-DBL/PSAMS) to aggregate data and we do the query based on clusterd data and their pivot points.

you need to change the clustered data path, cluster center path and other pivot points path to your own path in main.cpp. Data is separated by ','. One example can be found in ./inputFiles.

If you don't know how many clusters of data should be partitioned, you can use our evaluation indicators as a reference. We can score each cluster number from the start clusters number to the end of clusters number, the score will help you quickly determine the appropriate cluster number.

```bash
make clean
make all
./main [number of pivot point] [number of cluster] [dimention of point] [evaluate the clustering or not(0 /1 : [min number of cluster] [max number of cluster] [step size])]
```