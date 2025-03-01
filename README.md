# T-Mart: Topography-adjusted Monte-carlo Adjacency-effect Radiative Transfer Code

## Description 

T-Mart solves radiative transfer in a 3D surface-atmosphere system. It supports customizable surface models and enables simulation and correction for the adjacency effect (AE) in optical aquatic remote sensing. AE correction substantially improves satellite-based retrieval of water-leaving reflectance in nearshore environments (Wu et al., 2024). 


## Links


Home page: <a href="https://github.com/yulunwu8/tmart" target="_blank">https://github.com/yulunwu8/tmart</a>

User guide: <a href="https://tmart-rtm.github.io" target="_blank">https://tmart-rtm.github.io</a>

## Publications

Wu, Y., Knudby, A., & Lapen, D. (2023). Topography-adjusted Monte Carlo simulation of the adjacency effect in remote sensing of coastal and inland waters. *Journal of Quantitative Spectroscopy and Radiative Transfer*, 108589. <a href="https://doi.org/10.1016/j.jqsrt.2023.108589" target="_blank">https://doi.org/10.1016/j.jqsrt.2023.108589</a>

Wu, Y., Knudby, A., Pahlevan, N., Lapen, D., & Zeng, C. (2024). Sensor-generic adjacency-effect correction for remote sensing of coastal and inland waters. *Remote Sensing of Environment*, 315, 114433. <a href="https://doi.org/10.1016/j.rse.2024.114433" target="_blank">https://doi.org/10.1016/j.rse.2024.114433</a>


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

## Quick start: adjacency-effect correction 

T-Mart supports AE correction for Sentinel-2 MSI and Landsat 8/9 OLI/OLI-2 products. Correction is performed directly on level-1 products and can be followed by any amtospheric correction tools. 

Minimal input: 

```python
import tmart
file = 'user/test/S2A_MSIL1C_20160812T143752_N0204_R096_T20MKB_20160812T143749.SAFE'

# NASA EarthData Credentials, OB.DAAC Data Access needs to be approved
username = 'abcdef'
password = '123456'

# T-Mart uses multiprocessing, which needs to be wrapped in 'if __name__ == "__main__":' for Windows systems. This is optional for Unix-based systems
if __name__ == "__main__":
    tmart.AEC.run(file, username, password)
```

The tool takes approximately 20 min to process a Landsat 8/9 scene and 30 min for a Sentinel-2 scene on an eight-core personal computer. See <a href="https://tmart-rtm.github.io/ins_aec.html" target="_blank">Instruction - Adjacency-Effect Correction</a> for detailed instructions.


## Others

T-Mart can calculate reflectances of various units, see Table 1 in Wu et al. (2023) for examples. 

For questions and suggestions (which I'm always open to!), please open an issue or email Yulun at [yulunwu8@gmail.com](mailto:yulunwu8@gmail.com)

