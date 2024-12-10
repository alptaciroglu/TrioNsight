# TrioNsight

## Table of Contents
- [About](#-about)
- [Installation](#-Installation)
  - [Base](#base)
  - [VariPred](#varipred)
  - [SAAPdap](#saapdap)
  - [R packages](#r-packages)
- [Run TrioNsight](#run-trionsight)

## 🚀 About

**TrioNsight** - TrioN-like Dbl-homology domain specific variant effect predictor

## 📝 Installation 

### <ins>Base</ins>

```
git clone git@github.com:alptaciroglu/TrioNsight.git
cd TrioNsight

conda env create -f environment.yml
conda activate trionsight
```

### <ins>VariPred</ins>

For the VariPred repo please visit https://github.com/wlin16/VariPred

***Not mandatory if you'd like to use precalculated VariPred output used in the training of the TrioNsight -> Go to [SAAPdap](#saapdap)***

***If you need to run VariPred yourself, continue with VariPred installation directly below***

```
git clone https://github.com/wlin16/VariPred.git
cd VariPred
pip install -r requirements.txt
pip install torch==2.0.0 torchvision==0.15.1 torchaudio==2.0.1 --index-url https://download.pytorch.org/whl/cu118
pip install fair-esm

cd ../
cp -f VariPred_files/predict.sh VariPred/VariPred/predict.sh
# File contains adjusted paths only
cp -f VariPred_files/main.py VariPred/VariPred/main.py
# File contains adjusted paths only
```

### <ins>SAAPdap</ins>

***Not mandatory if you'd like to use precalculated SAAPdap output used in the training of the TrioNsight -> Go to Run TrioNsight [Run TrioNsight](#run-trionsight)***

***If you need to run SAAPdap yourself, continue with SAAPdap installation directly below***

https://github.com/ACRMGroup/SAAP3

Please go to SAAP3 repository and follow the installation instructions.
SAAPdap should be installed under TrioNsight/structure_analysis/. Alternatively, modify the paths under TrioNsight/scripts/SAAPdap.R wrapper function based on chosen SAAPdap installation directory.

### <ins>R packages</ins>

```
R
source('scripts/helper.R')
```

## <ins>Run TrioNsight</ins>
