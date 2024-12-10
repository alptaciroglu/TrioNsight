# TrioNsight

<h4 align="center">A specialised variant effect predictor for TrioN-like Dbl-homology domains.</h4>

<p align="center">
    <a href="https://github.com/alptaciroglu/TrioNsight/commits/master">
    <img src="https://img.shields.io/github/last-commit/alptaciroglu/TrioNsight.svg?style=flat-square&logo=github&logoColor=white"
         alt="GitHub last commit">
    <a href="https://github.com/alptaciroglu/TrioNsight/issues">
    <img src="https://img.shields.io/github/issues-raw/alptaciroglu/TrioNsight.svg?style=flat-square&logo=github&logoColor=white"
         alt="GitHub issues">
</p>
      
<p align="center">
  <a href="#-installation">🔧 Installation</a>
  <a href="#-run-trionsight">🚀 Run TrioNsight</a>
  <a href="#-input-data">📊 Input Data</a>
  <a href="#-citing-this-work">❛❛❞ Citing this work</a>
  <a href="#-data-availability">✔ Data availability</a>
  <a href="#-acknowledgement">🤝 Acknowledgement</a>
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

## 🔧 Installation 

### <ins>Base</ins>

```
git clone git@github.com:alptaciroglu/TrioNsight.git
cd TrioNsight

conda env create -f environment.yml
conda activate trionsight
```

### <ins>VariPred</ins>

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

### <ins>SAAPdap</ins>

For the SAAPdap repo please visit https://github.com/ACRMGroup/SAAP3

> [!NOTE]  
> Not mandatory if you'd like to use precalculated SAAPdap output used in the training of the TrioNsight -> Go to [Run TrioNsight](#run-trionsight)
>
> If you need to run SAAPdap yourself, continue with SAAPdap installation directly below

Please go to SAAP3 repository and follow the installation instructions.
SAAPdap should be installed under TrioNsight/structure_analysis/. Alternatively, modify the paths under TrioNsight/scripts/SAAPdap.R wrapper function based on chosen SAAPdap installation directory.

### <ins>R packages</ins>

```
R
source('scripts/helper.R')
```

## 🚀 Run TrioNsight

## 📊 Input Data

## ❛❛❞ Citing this work

## ✔ Data availability

## 🤝 Acknowledgement
