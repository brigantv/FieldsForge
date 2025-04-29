module parameters
implicit none

 double precision, parameter     :: PI=4.D0*DATAN(1.D0)
 double precision, parameter     :: A_to_B= 1.8897259886d0
 double precision, parameter     :: Har_to_Kc=627.503d0
 double precision, parameter     :: F_conv=Har_to_Kc*A_to_B
 double precision, parameter     :: boltz=3.11811E-06
 double precision, parameter    :: amu_to_emass=1822.89d0

end module parameters

