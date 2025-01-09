## Cell Cycle Prediction with CNN
These scripts predict the cell cycle assignment of specially formatted imaging data where each file contains a 3D brightfield image of a cell

## 1) Install requirements (tested with conda package manager, python 3.10)
```bash
pip install -r requirements.txt
```
## 2) (optional) train model:
```bash
python runModel.py -m train -epochs N_epochs -dataDir dataPATH
```

## 3) test model/generate bottleneck features (reproduce results in Alahmari et al.):
```bash
python runModel.py -dataDir dataPATH -m test -testDir ./ -arch custom_cnn -exp CNN_exp1_merged_data_ch1
```
