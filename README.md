# diagnostics of the Q-G omega equation 
This contains a set of Fortran codes used to diagnose the contributions of three processes to the long-term trends in vertical motion over the region of interest, based on the traditional Q-G omega equation (Eq. 1 and the equation from Nigam et al. 2000) and ERA5 daily (or 6-hourly) data from 1980 to 2020:

<img width="1680" height="204" alt="image" src="https://github.com/user-attachments/assets/ff507f2e-47e0-48fc-8eb8-2c419f2aebef" />

 
omega.term1.f is a Fortran code to calcualte term1 ( Vg∙∇p (ζg+f) in our Eq.1 ) using daily ERA5 data. 
More information on the calculation procedure is provided in the code comments.
 
omega.term2.f is a Fortran code to calcualte term2 (Vg∙∇p T in our Eq.1) using daily ERA5 data. 
More information on the calculation procedure is provided in the code comments.


To calculate the diabatic heating term Q ( term 3 in Eq.1), we follow the equation described in Nigam et al. (2000), which is adapted here.
 


<img width="300" height="111" alt="image" src="https://github.com/user-attachments/assets/e3d041b0-0da5-4e67-8541-591c9038d7de" />



Q = term A + term B + term C + term D + term E, where the overbar denotes the monthly average, and the prime indicates the deviation of the 6-h data from the monthly average.

diabatic.term3.f is a Fortran code to calcualte Q using 6-h ERA5 data.
More information on the calculation procedure is provided in the code comments.


lineartrend.f is used to calculate linear trends in our study.


Reference
Nigam, S., C. Chung, and E. DeWeaver, 2000: ENSO Diabatic Heating in ECMWF and NCEP–NCAR Reanalyses, and NCAR CCM3 Simulation. J. Climate, 13, 3152–3171,

