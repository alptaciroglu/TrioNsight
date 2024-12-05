## TrioNsight

## Table of Contents
- [About](#-about)
- [Installation](#-step-1)

## 🚀 About

**TrioNsight** - TrioN-like Dbl-homology domain specific variant effect predictor

## 📝 Step 1: Installation

```
git clone git@github.com:alptaciroglu/TrioNsight.git
cd TrioNsight

conda env create -f VariPred_files/environment.yml
conda activate trionsight
```
### Step 1.1: Necessary R libraries

```
R
source('scripts/helper.R')
```

### Step 1.2: VariPred

***Not mandatory if you'd like to use precalculated VariPred output used in the training of the TrioNsight -> Go to Step 1.3 [here](#step-1.3) ***

***If you need to run VariPred yourself, continue with VariPred installation directly below***

https://github.com/wlin16/VariPred
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

### Step 1.3: SAAPdap

***Not mandatory if you'd like to use precalculated SAAPdap output used in the training of the TrioNsight -> Go to Step 2***

***If you need to run SAAPdap yourself, continue with SAAPdap installation directly below***

https://github.com/ACRMGroup/SAAP3

Please go to SAAP3 repository and follow the installation instructions.
SAAPdap should be installed under TrioNsight/structure_analysis/. Alternatively, modify the paths under TrioNsight/scripts/SAAPdap.R wrapper function based on chosen SAAPdap installation directory.



### Step 2: Run TrioNsight
