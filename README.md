# diagnostics of the Q-G omega equation
This contains a set of Fortran codes used to diagnose the contributions of three processes to the long-term trends in vertical motion over the region of interest, based on the traditional Q-G omega equation (Eq1) and ERA5 daily (or 6-hourly) data

<img width="1680" height="204" alt="image" src="https://github.com/user-attachments/assets/ff507f2e-47e0-48fc-8eb8-2c419f2aebef" />

 
 omega.term1.f is a Fortran code to calcualte term1 (Vg∙∇p (ζg+f)) in our Eq.1. more informaiton of the calcualtioon procedue are commented in the code
 omega.term2.f is a Fortran code to calcualte term2 (Vg∙∇p T) in our Eq.1. more informaiton of the calcaultion procedue are commented in the code

To clacuate the diabtic heating term ( term 3 in Eq.1), we follow oow the equation descibe in Nigam et al. 2010., whcih is adpated here 
 


<img width="300" height="111" alt="image" src="https://github.com/user-attachments/assets/e3d041b0-0da5-4e67-8541-591c9038d7de" />
Q=term A+ term B + term C + term D + term E
diabtic.f is a Fortran code to calcualte Q.  more informaiton of the calcaultion procedue are commented in the code

list.f is a Fortan code to calcualte linear trend , whci extneve used in our study
