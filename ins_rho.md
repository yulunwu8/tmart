# Instruction - Modelling Sea-Surface Reflectance factor

The sea-surface reflectance factor, *ρ*, is used to remove skyglint and sunglint from above-water measurements of water's reflectance. 


To calculate *ρ* with a solar zenith angle of 30°, viewing zenith angle of 40° and a relative azimuth angle of 135°: 

```python
import tmart
rho = tmart.surface_rho.calculate(wl=400, viewing_zenith=40, solar_zenith=30, relative_azimuth=135)
print(rho)
```

Note that Windows users need to wrap the code within \<\< if \_\_name\_\_ == "\_\_main\_\_" \>\> because it is required to call functions that use the *multiprocessing* library in Windows. 

To change the default AOT<sub>550</sub> of 0.05 and wind speed of 3 m/s: 

```python
rho = tmart.surface_rho.calculate(wl=400, viewing_zenith=40, solar_zenith=30, relative_azimuth=135, aot550=0.05, wind_speed = 8)
print(rho)
```

To calculate the values of *ρ* from 400 to 800 nm at an interval of 100 nm: 

```python
rho = tmart.surface_rho.calculate(wl=[400,800,100], viewing_zenith=40, solar_zenith=30, relative_azimuth=135, aot550=0.05, wind_speed = 8)
print(rho)
```

Other parameters of the function include atmosphere type, aerosol type, number of photons for computation. See full parameters and their descriptions <a href="https://tmart-rtm.github.io/tmart.html#module-tmart.surface_rho.calculate" target="_blank">HERE</a>.








