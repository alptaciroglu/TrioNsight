# TrioNsight

<h4 align="center">A specialised variant effect predictor for TrioN-like Dbl-homology domains.</h4>

<p align="center">
    <a href="https://github.com/alptaciroglu/TrioNsight/commits/master">
    <img src="https://img.shields.io/github/last-commit/alptaciroglu/TrioNsight.svg?style=flat-square&logo=github&logoColor=white"
         alt="GitHub last commit">
    <a href="https://github.com/alptaciroglu/TrioNsight/issues">
    <img src="https://img.shields.io/github/issues-raw/alptaciroglu/TrioNsight.svg?style=flat-square&logo=github&logoColor=white"
         alt="GitHub issues">
    <a href="https://github.com/alptaciroglu/TrioNsight/pulls">
    <img src="https://img.shields.io/github/issues-pr-raw/alptaciroglu/TrioNsight.svg?style=flat-square&logo=github&logoColor=white"
         alt="GitHub pull requests">
</p>
      
<p align="center">
  <a href="#-1-installation">1. Installation</a> •
    <a href="#11-base">1.1. Base</a> •
    <a href="#12-varipred">1.2. VariPred</a> •
    <a href="#13-saapdap">1.3. SAAPdap</a> •
    <a href="#14-r-packages">1.4. R Packages</a> •
  <a href="#-2-run-trionsight">2. Run TrioNsight</a> •
</p>

---

<table>
<tr>
<td>
  
**TrioNsight** is a variant effect predictor trained on TrioN-like DH domains - N-terminal Dbl-homology domain of TRIO (TrioN) and Dbl-homology domains that are highly similar to TrioN.


![workflow](workflow.png)
<p align="right">
<sub>(Preview)</sub>
</p>

</td>
</tr>
</table>

## ⚙️ 1. Installation 

### <ins>1.1. Base</ins>

```
git clone git@github.com:alptaciroglu/TrioNsight.git
cd TrioNsight

conda env create -f environment.yml
conda activate trionsight
```

### <ins>1.2. VariPred</ins>

For the VariPred repo please visit https://github.com/wlin16/VariPred

> [!NOTE]  
> Not mandatory if you'd like to use precalculated VariPred output used in the training of the TrioNsight -> Go to [SAAPdap](#saapdap)
>
> If you need to run VariPred yourself, continue with VariPred installation directly below

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

### <ins>1.3. SAAPdap</ins>

> [!NOTE]  
> Not mandatory if you'd like to use precalculated SAAPdap output used in the training of the TrioNsight -> Go to Run TrioNsight [Run TrioNsight](#run-trionsight)
>
> If you need to run SAAPdap yourself, continue with SAAPdap installation directly below

https://github.com/ACRMGroup/SAAP3

Please go to SAAP3 repository and follow the installation instructions.
SAAPdap should be installed under TrioNsight/structure_analysis/. Alternatively, modify the paths under TrioNsight/scripts/SAAPdap.R wrapper function based on chosen SAAPdap installation directory.

### <ins>1.4. R packages</ins>

```
R
source('scripts/helper.R')
```

## 🚀 2. Run TrioNsight
