# Eddy-Covariance-Partitioning
This is an example code to partition net ecosystem CO2 exchange from eddy covariance into gross primary productivity (GPP) 
and ecosystem respiration (Reco) using by fitting a light response curve with vapor pressure deficit (VPD) limitation (Lasslop et al. 2010). 
For this, we used daytime data (observed values between 08:00 – 19:00, no gap-filled data was included), 
using 7-day moving windows (Equation 1). The model intercept was used as an estimate of Reco.

### Equation 1:
```math
NEE=VPD>{VPD}_0,\ -\beta\exp(-k\left(VPD-{VPD}_0\right)\left(\frac{1-exp{\left(-\alpha PAR\right)}}{\beta\exp(-k(VPD-{VPD}_0)}\right)+\gamma\ ,
```
```math
VPD\le{VPD}_0,\ -\left(\beta+\gamma\right)(1-exp\left(\frac{aPAR}{\beta+\gamma}\right)+\gamma)
```

Whereby, 
$\beta\$ is the maximum CO2 uptake rate of the canopy at light saturation, 
$\alpha\$ is canopy light utilization which represents the initial slope of the light–response curve, 
$\gamma\$ is ecosystem respiration and PAR is photosynthetically active radiation. 

Parameter k was set to a constant value, which was derived from a separate VPD response curve. The VPD response curve was fitted to the whole 
dataset but using only light-saturated (photosynthetically active radiation, PAR>1200) 
and VPD limited (>VPD0) conditions, assuming aVPD0 threshold of 10hPa (Körner, 1995)
A constant value for \alpha was applied which was the slope of the linear regression between NEE and PAR under 
low light conditions (<200 μmol) (Xu et al., 2019).


### References
Ch. Körner, “Leaf diffusive conductances in the major vegetation types of the globe” in Ecophysiology of Photosynthesis, 100th 
Ed., E. D. Schulze, M. M. M.Caldwell, Eds. (Springer Study Edition, 1995), pp. 463–490.

G. Lasslop, et al., Separation of net ecosystem exchange into assimilation and respiration using a light response curve approach: 
critical issues and global evaluation. Glob Chang Biol 16, 187–208 (2010).

J. Xu, et al., A general non-rectangular hyperbola equation for photosynthetic light response curve of rice at various 
leaf ages. Sci Rep 9, 1–8 (2019)
