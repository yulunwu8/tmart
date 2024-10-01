# Overview

## T-Mart (Topography-adjusted Monte-carlo Adjacency-effect Radiative Transfer code)

T-Mart solves radiative transfer in a 3D surface-atmosphere system. This can be used in modelling and correcting for the adjacency effect (AE). Although initially designed for aquatic environments, T-Mart can also handle simple terrestrial applications.


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
conda install -c conda-forge Py6S
```

3 - Install tmart: 

```bash
pip3 install tmart
```

## Publications

Wu, Y., Knudby, A., & Lapen, D. (2023). Topography-Adjusted Monte Carlo Simulation of the Adjacency Effect in Remote Sensing of Coastal and Inland Waters. *Journal of Quantitative Spectroscopy and Radiative Transfer*, 108589. <a href="https://doi.org/10.1016/j.jqsrt.2023.108589" target="_blank">https://doi.org/10.1016/j.jqsrt.2023.108589</a>

Wu, Y., Knudby, A., Pahlevan, N., Lapen, D., & Zeng, C. (2024). Sensor-generic adjacency-effect correction for remote sensing of coastal and inland waters. *Remote Sensing of Environment*, 315, 114433. <a href="https://doi.org/10.1016/j.rse.2024.114433" target="_blank">https://doi.org/10.1016/j.rse.2024.114433</a>

--

Home page: <a href="https://github.com/yulunwu8/tmart" target="_blank">https://github.com/yulunwu8/tmart</a>

Yulun Wu | September 29, 2024 | [yulunwu8@gmail.com](mailto:yulunwu8@gmail.com)







