# diagnostics of the Q-G omega equation 
This contains a set of Fortran codes used to diagnose the contributions of three processes to the long-term trends in vertical motion over the region of interest, based on the traditional Q-G omega equation (Eq. 1: slightly modified from Eq. 6.34 in Holton 2004, and the formulation from Nigam et al. 2000) and ERA5 daily (or 6-hourly) data from 1980 to 2020:

<img width="1680" height="204" alt="image" src="https://github.com/user-attachments/assets/ff507f2e-47e0-48fc-8eb8-2c419f2aebef" />

 
omega.term1.f is a Fortran code to calcualte term1 on a monthly basis ( Vg∙∇p (ζg+f) in our Eq.1 ) using daily ERA5 data. 
More information on the calculation procedure is provided in the code comments.
 
omega.term2.f is a Fortran code to calcualte term2 on a monthly basis (Vg∙∇p T in our Eq.1) using daily ERA5 data. 
More information on the calculation procedure is provided in the code comments.


To calculate the diabatic heating term Q ( term 3 in Eq.1), we follow the equation described in Nigam et al. (2000), which is adapted here.
 


<img width="300" height="111" alt="image" src="https://github.com/user-attachments/assets/e3d041b0-0da5-4e67-8541-591c9038d7de" />



Q = term A + term B + term C + term D + term E, where the overbar denotes the monthly average, and the prime indicates the deviation of the 6-h data from the monthly average.

diabatic.term3.f is a Fortran code to calcualte monthly Q fields using 6-h ERA5 data.
More information on the calculation procedure is provided in the code comments.

Based on these monthly fields of terms 1 to 3, we further calculate their linear trends in the annual mean fields using various conventional approaches, which are not detailed here.

References

Holton, J. R., 2004: An Introduction to Dynamic Meteorology (4th ed.). Academic Press.

Nigam, S., C. Chung, and E. DeWeaver, 2000: ENSO Diabatic Heating in ECMWF and NCEP–NCAR Reanalyses, and NCAR CCM3 Simulation. J. Climate, 13, 3152–3171,

