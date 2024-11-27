# TrioNsight
Codes and results for TrioNsight - TrioN-like Dbl-homology domain specific variant effect predictor

## Step 1: Installation
```
git clone git@github.com:alptaciroglu/TrioNsight.git
cd TrioNsight

conda env create -f environment.yml
conda activate trionsight
```
## Step 1.1: Include VariPred installation
*** Not mandatory if you'd like to use precalculated VariPred output used in the training of the TrioNsight -> Go to Step 2. ***
*** If you need to use VariPred, continue here. ***
```
pip install -r varipred_requirements.txt
pip install torch==2.0.0 torchvision==0.15.1 torchaudio==2.0.1 --index-url https://download.pytorch.org/whl/cu118
pip install fair-esm
```
## Step 2:
