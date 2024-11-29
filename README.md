# TrioNsight
Codes and results for TrioNsight - TrioN-like Dbl-homology domain specific variant effect predictor

## Step 1: Installation

```
git clone git@github.com:alptaciroglu/TrioNsight.git
cd TrioNsight

conda env create -f VariPred_files/environment.yml
conda activate trionsight
```

## Step 1.1: VariPred
***Not mandatory if you'd like to use precalculated VariPred output used in the training of the TrioNsight -> Go to Step 1.2***
***If you need to run VariPred yourself, continue directly from below***

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

## Step 1.2: SAAPdap
***Not mandatory if you'd like to use precalculated SAAPdap output used in the training of the TrioNsight -> Go to Step 1.3***
***If you need to run SAAPdap yourself, continue directly from below***

Please go to SAAP3 repository and follow the installation instructions.
SAAPdap should be installed under TrioNsight/StructureAnalysis or modify the paths under TrioNsight/Scripts/SAAPdap.R wrapper function based on SAAPdap installation directory.

https://github.com/ACRMGroup/SAAP3

## Step 1.3: Necessary R libraries

```
R
source('Scripts/helper.R')
```

## Step 2: Run TrioNsight
