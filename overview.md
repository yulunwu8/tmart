# Overview

## T-Mart (Topography-adjusted Monte-carlo Adjacency-effect Radiative Transfer code)

T-Mart solves radiative transfer in a 3D surface-atmosphere system. This can be used in modelling and correcting for the adjacency effect. Although initially designed for aquatic environments, T-Mart can also handle simple terrestrial applications.


## AE correction

See Tab <a href="https://tmart-rtm.github.io/ins_aec.html" target="_blank">Instruction - Adjacency-Effect Correction</a> for details.


## AE modelling 

See Tab <a href="https://tmart-rtm.github.io/ins_basic.html" target="_blank">Instruction - Basic Modelling</a> and <a href="https://tmart-rtm.github.io/ins_advanced.html" target="_blank">Instruction - Advanced Modelling</a> for details.


## Installation 

1 - Create a conda environment and activate it: 

```bash
conda create --name tmart python=3.9
conda activate tmart
```

2 - Install dependencies: 

```bash
conda install -c conda-forge Py6S rasterio==1.3.9
```

3 - Install tmart: 

```bash
pip3 install tmart
```


---

Home page: <a href="https://github.com/yulunwu8/tmart" target="_blank">https://github.com/yulunwu8/tmart</a>

Yulun Wu | February 28, 2025 | [yulunwu8@gmail.com](mailto:yulunwu8@gmail.com)







