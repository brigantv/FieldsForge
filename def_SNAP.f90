module SNAP_class

use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
use lammps_class
use max_class
use parameters
use lapack_inverse
implicit none

contains

subroutine fit_energy_forces(set,num_bisp,energies,DFT_forces,lambda,ezero_constr,set_type,Y,flag_energy,flag_forces,weight)
implicit none

double precision, dimension(:,:), allocatable                           :: A,C
double precision, dimension(:,:),allocatable                            :: F
double precision, dimension(:),allocatable                              :: b,forces_predicted
double precision, dimension(:),allocatable,intent(inout)                :: energies,DFT_forces
character(len=5),intent(in)                                             :: set_type
double precision                                                        :: ave_atom
double precision,intent(in)                                             :: weight
integer                                                                 :: tot_kinds,tot_atom,comp
integer                                                                 :: i,j,k,l,m
integer,intent(in)                                                      :: num_bisp
integer                                                                 :: size_ref
integer                                                                 :: counter,count_kinds_ezero,ez_cons_rows
integer                                                                 :: tot_atoms
type(lammps_obj),dimension(:),allocatable,intent(in)                    :: set
logical,dimension(:),allocatable,intent(in)                             :: ezero_constr
integer                                                                 :: start_cycle,end_cycle,help_counter
integer                                                                 :: start_snap_force
character(len=1)                                                        :: TRANS
integer                                                                 :: MF,N,NRHS,LDA,LDB,LWORK,INFO
double precision,dimension(:),allocatable                               :: WORK
double precision,dimension(:),allocatable,intent(out)                   :: Y
double precision,intent(in)                                             :: lambda
double precision, dimension(:,:), allocatable                           :: D,D_t,test,test_copy
logical,intent(in)                                                      :: flag_forces,flag_energy

call get_tot_kinds(set,tot_kinds)
call get_ave_atoms(set,ave_atom,tot_atom)

ez_cons_rows = count(ezero_constr)

if ((flag_energy).and.(.not.flag_forces)) then

        size_ref=size(set) + ez_cons_rows

end if

if ((flag_forces).and.(.not.flag_energy)) then
        
        size_ref=3*tot_atom + ez_cons_rows

end if

if ((flag_forces).and.(flag_energy)) then

        size_ref=size(set)+3*tot_atom+ez_cons_rows

end if

if (lambda.ne.0.0) then

  size_ref=size_ref+(num_bisp-1)*tot_kinds

end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(A(size_ref,num_bisp*tot_kinds))
A=0.0

allocate(b(size_ref))
b=0.0

if ((flag_energy).and.(.not.flag_forces)) then

 b(1:size(set))=energies

end if 

if ((flag_forces).and.(.not.flag_energy)) then

 b(1:3*tot_atom)=DFT_forces

end if

if ((flag_forces).and.(flag_energy)) then
        
        b(1:size(set))=weight*energies
        b(size(set)+1:size(set)+3*tot_atom)=DFT_forces

end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SOME VARIABLES FOR REGULARIZATION 

if ((flag_energy).and.(.not.flag_forces)) then

 start_cycle=size(set)+1       
 end_cycle=size_ref-ez_cons_rows
 help_counter=size(set)

end if

if ((flag_forces).and.(.not.flag_energy)) then

 start_cycle=3*tot_atom+1
 end_cycle=size_ref-ez_cons_rows
 help_counter=3*tot_atom
end if

if ((flag_forces).and.(flag_energy)) then

 start_cycle=size(set)+3*tot_atom+1       
 end_cycle=size_ref-ez_cons_rows
 help_counter=size(set)+3*tot_atom

end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!REGULARIZATION

if (lambda.ne.0.0) then

 counter=0

  do j=start_cycle,end_cycle


   if (mod(j-help_counter,num_bisp-1)==1) then

    counter=counter+1
  
   end if

   A(j,j-help_counter+counter) = sqrt(lambda)

 end do

end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!BUILDING THE SNAP MATRIX FOR ENERGY

if ((flag_energy).or.((flag_energy).and.(flag_forces))) then

do j=1,size(set)

 do i=1,set(j)%nats

  do k=1,num_bisp

   A(j,(set(j)%kind(i)-1)*num_bisp+k) = A(j,(set(j)%kind(i)-1)*num_bisp+k) + set(j)%at_desc(i)%desc(k)

  end do

 end do

end do

allocate(C(size(set),num_bisp*tot_kinds))

C=A(1:size(set),:)
A(1:size(set),:)=A(1:size(set),:)*weight
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FLAG E0 TO REGULARIZE THE FIRST COEFFICIENTS OF BISPECTRUM

counter=0

do j=1,tot_kinds

 if (ezero_constr(j)) then
     
   counter=counter+1
   A(size_ref-ez_cons_rows+counter,(j-1)*num_bisp+1)=1
 
 end if

end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CONSTRUCTING THE SNAP MATRIX OF FORCES
if ((flag_forces).and.(.not.flag_energy)) then
        start_snap_force=0
else if ((flag_forces).and.(flag_energy)) then
        start_snap_force=size(set)
end if
 
if (flag_forces) then

 do j=1,size(set)

  do i=1,set(j)%nkinds

   do l=1,set(j)%nats
                        
    do comp=1,3

     do k=1,num_bisp-1 
   
        A(start_snap_force+(j-1)*set(j)%nats*3+(l-1)*3+comp,num_bisp*(i-1)+1+k)=&
                set(j)%bisp_der((i-1)*((num_bisp-1)*3)+(num_bisp-1)*(comp-1)+k,l)

     end do

    end do
  
   end do

  end do

 end do

end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (flag_forces) then
        
        allocate(F(3*tot_atom,num_bisp*tot_kinds))
        F=A(start_snap_force+1:start_snap_force+3*tot_atom,:)

end if

!!!!DEBUG
!open(11,file='SNAP_matrix',action='write')
! do i=1,size(A,1)
!  write(11,*) A(i,:)
! end do
!close(11)

!open(11,file='target_values',action='write')
! do i=1,size(b,1)
!  write(11,*) b(i)
! end do
!close(11)

!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SOLVING THE LINEAR LEAST SQUARES PROBLEM
if (set_type=='TRAIN') then
        
TRANS='N'
M=size_ref
LDA=size_ref
N=tot_kinds*num_bisp
NRHS=1
LDB=max(M,N)
LWORK=2*min(M,N)
allocate(WORK(LWORK))

call dgels(TRANS,M,N,NRHS,A,LDA,b,LDB,WORK,LWORK,INFO)
deallocate(A)

if (info.ne.0) then

 write(*,*) 'Convergence issues: could not solve the linear least square problem (def_SNAP.F90/fit_energy)'

end if

open(11,file='snapcoeff_energy',action='write')

 do i=1,num_bisp*tot_kinds

  write(11,*) b(i)

 end do

close(11)

b=b(1:num_bisp*tot_kinds)

else

deallocate(b)
allocate(b(num_bisp*tot_kinds))

open(11,file='snapcoeff_energy',action='read')

 do i=1,tot_kinds*num_bisp
  
  read(11,*) b(i)

 end do

close(11)

end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CALCULATING THE PREDICTED ENERGIES

if (flag_energy) then

allocate(Y(size(set)))
Y=matmul(C,b)
deallocate(C)

end if

if (flag_forces) then
        
  allocate(forces_predicted(3*tot_atom))
  forces_predicted=matmul(F,b)
  deallocate(F)
end if  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PRINTING TO FILES COEFFICIENTS OF FIT AND ENERGIES PREDICTED

if (flag_energy) then

open(11,file='energy_rms_fpp.dat')

write(11,*) '#','      ','Num config','       ','ML energy','             ','DFT energy','           ','Error'

do i=1,size(set)
 write(11,*) ' ', i, Y(i), energies(i),Y(i)-energies(i)
end do

write(11,*)'# RMS=',(1/ave_atom)*(sqrt(sum((Y-energies)**2)/size(set))), 'kcal/mol/atom'

close(11)

end if

if (flag_forces) then

open(11,file='forces_rms_fpp.dat')

write(11,*) '#','      ','Num config','       ','ML force','             ','DFT force','           ','Error'

do i=1,3*tot_atom
 write(11,*) ' ', i, forces_predicted(i), DFT_forces(i),forces_predicted(i)-DFT_forces(i)
end do

write(11,*)'# RMS=',(sqrt(sum((forces_predicted(:)-DFT_forces(:))**2)&
        /(dble(3*tot_atom)))), 'kcal/mol/angstrom'

close(11)

end if

deallocate(b)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine fit_energy_forces

subroutine fit_dipoles(set,num_bisp,dipoles_file,lambda,ezero_constr,set_type,tot_charge,len_file_dip)

type(lammps_obj),dimension(:),allocatable,intent(inout)   :: set
integer, intent(in)                                       :: num_bisp
character(len=*),optional,intent(in)                      :: dipoles_file
character(len=150)                                        :: mulliken
character(len=150)                                        :: dip_file
character(len=5),intent(in)                               :: set_type
double precision,dimension(:),allocatable,intent(in)      :: tot_charge
double precision, dimension(:,:), allocatable             :: A,C
logical,dimension(:),allocatable,intent(in)               :: ezero_constr
integer                                                   :: tot_kinds
integer                                                   :: i,j,k,l,m
integer                                                   :: tot_atoms
double precision, dimension(:), allocatable               :: b
double precision,dimension(:),allocatable                 :: dipoles
integer                                                   :: size_ref
integer,intent(in)                                        :: len_file_dip
double precision                                          :: ave_atom
integer                                                   :: counter
integer                                                   :: ezero_rows
!variables to call function dgels
character(len=1)                                          :: TRANS
integer                                                   :: MF,N,NRHS,LDA,LDB,LWORK,INFO
double precision,dimension(:),allocatable                 :: WORK
!variables to call dgemv
double precision                                          :: ALPHA,BETA
integer                                                   :: INCX,INCY
double precision,dimension(:),allocatable                 :: Y
double precision,intent(in)                               :: lambda
double precision                                          :: charge_mol_1,charge_mol_2


dip_file=dipoles_file(1:len_file_dip)

call get_tot_kinds(set,tot_kinds)
call get_ave_atoms(set,ave_atom,tot_atoms)

ezero_rows=count(ezero_constr)


if (lambda.ne.0.0) then

  size_ref=3*size(set)+(num_bisp-1)*tot_kinds+size(set)+ezero_rows

else

 size_ref=3*size(set)+size(set)+ezero_rows

end if

allocate(A(size_ref,num_bisp*tot_kinds))
A=0.0

allocate(b(max(size_ref,num_bisp*tot_kinds)))
b=0.0

open(10,file=dip_file)

!initialize the array b with the values of dipoles

do l=1,3*size(set)

 read(10,*) b(l)

end do

close(10)

counter=0

do l=size_ref-size(set)-ezero_rows+1,size_ref-ezero_rows
 
 counter=counter+1
 b(l)=tot_charge(counter)

end do

close(10)

allocate(dipoles(3*size(set)))
dipoles=b(1:3*size(set))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!BUILDING THE SNAP MATRIX

do j=1,size(set)

 do m=1,3

  do i=1,set(j)%nats

   do k=1,num_bisp

    A(3*(j-1)+m,(set(j)%kind(i)-1)*num_bisp+k) = A(3*(j-1)+ m,(set(j)%kind(i)-1)*num_bisp + k)&
           + set(j)%at_desc(i)%desc(k)*set(j)%x(i,m)

   end do

  end do

 end do

end do

A=A_to_B*A

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (lambda.ne.0.0) then

counter=0

 do j=3*size(set)+1,size_ref-size(set)-ezero_rows

  if (mod(j-3*size(set),num_bisp-1)==1) then
    counter=counter+1
  end if

   A(j,j-3*size(set)+counter) = sqrt(lambda)

 end do

end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CHARGE CONSERVATION
counter=0

do j=size_ref-size(set)-ezero_rows+1,size_ref-ezero_rows
 
 counter=counter+1
 do i=1,set(counter)%nats

  do k=1,num_bisp
  
   A(j,(set(counter)%kind(i)-1)*num_bisp+k) = A(j,(set(counter)%kind(i)-1)*num_bisp+k) + set(counter)%at_desc(i)%desc(k)

  end do

 end do

end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FLAG EO FOR DIPOLES

counter=0

do j=1,tot_kinds

 if (ezero_constr(j)) then
  counter=counter+1
  A(size_ref - ezero_rows + counter,(j-1)*num_bisp + 1)= 1

end if

end do


allocate(C(3*size(set),num_bisp*tot_kinds))

C=A(1:3*size(set),:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!DEBUGGER
!open(222,file='electro_snap_matrix',position="append")

!do j=1,size_ref

! write(222,*) A(j,:)

!end do

!close(222)

!open(222,file='target_snap_matrix',position="append")

!do j=1,size_ref

! write(222,*) b(j)

!end do

!close(222)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SOLVING THE LINEAR LEAST SQUARES PROBLEM
if (set_type=='TRAIN') then

 TRANS='N'
 M=size_ref
 LDA=size_ref
 N=tot_kinds*num_bisp
 NRHS=1
 LDB=max(M,N)
 LWORK=2*min(M,N)
 allocate(WORK(LWORK))

 call dgels(TRANS,M,N,NRHS,A,LDA,b,LDB,WORK,LWORK,INFO)
 if (info.ne.0) then
  write(*,*) 'Convergence issues: could not solve the linear least square problem for the charges'
 end if
 deallocate(A)
 open(11,file='snapcoeff_dipoles',action='write')

 do i=1,num_bisp*tot_kinds

  write(11,*) b(i)

 end do

 close(11)


else

deallocate(A)

b=0.0

 open(11,file='snapcoeff_dipoles',action='read')

 do i=1,num_bisp*tot_kinds

  read(11,*) b(i)
 
 end do

 close(11)

end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CALCULATING THE PREDICTED DIPOLES


b=b(1:num_bisp*tot_kinds)
allocate(Y(3*size(set)))
Y=matmul(C,b)
deallocate(C)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PRINTING TO FILES COEFFICIENTS FOR DEBUGGING OF FIT AND DIPOLES PREDICTED
open(11,file='dipoles_rms_fpp.dat',action='write')

write(11,*) '#','      ','Num config','       ','ML dipole','             ','DFT dipole','           ','Error'

do i=1,3*size(set)
 write(11,*) ' ', i, Y(i), dipoles(i),Y(i)-dipoles(i)
end do

write(11,*)'# RMS=',(1/ave_atom)*sqrt(sum((Y-dipoles)**2)/(3*size(set))), 'a.u./atom (dipoles)'

close(11)

deallocate(dipoles,Y)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CALCULATING CHARGES FROM DIPOLES

do i=1,size(set)
  
if (.not. allocated(set(i)%charge)) allocate(set(i)%charge(set(i)%nats))

 set(i)%charge=0.0

 do j=1,set(i)%nats
  
  do m=1,num_bisp

   set(i)%charge(j)=set(i)%charge(j)+set(i)%at_desc(j)%desc(m)*b((set(i)%kind(j)-1)*num_bisp+m)
  
  end do
 
 end do

end do

deallocate(b)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DEBUGGER CHARGES

!open(222,file='snap_charges',action='write',position="append")

!do j=1,size(set)
 
! write(222,*) 'Total charge:',sum(set(j)%charge)

!end do

!close(222)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine fit_dipoles

subroutine SNAP_coulomb_energy(set,R_screen,coul_energy)
double precision, intent(in)                                 :: R_screen
integer                                                      :: i,j,k
type(lammps_obj),dimension(:),allocatable,intent(in)         :: set
double precision, dimension(:),allocatable,intent(out)       :: coul_energy
double precision                                             :: screen
double precision,dimension(:,:),allocatable                  :: rel_dist

allocate(coul_energy(size(set)))
coul_energy=0.0

do i=1,size(set)
 
 allocate(rel_dist(set(i)%nats,set(i)%nats))
 rel_dist=0.0

 call set(i)%dist_ij 

 do j=1,set(i)%nats
 
  do k=1,set(i)%nats

   if (k>j) then
    
    rel_dist(j,k)=A_to_B*set(i)%dist(j,1,k,1)

    if (rel_dist(j,k) <= R_screen) then 
            
     screen=0.5*(1-dcos(PI*rel_dist(j,k)/R_screen))
    
    else
    
     screen=1.0

    end if     

    coul_energy(i)=coul_energy(i)+ screen*(set(i)%charge(j)*set(i)%charge(k)/rel_dist(j,k))
    
   end if

  end do

 end do

 deallocate(rel_dist)

end do

coul_energy=coul_energy*Har_to_Kc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!DEBUGGER ELECTROSTATIC ENERGY
!open(222,file='snap_electro_energy')

! do j=1,size(set)
  
 ! write(222,*) coul_energy(j)
 
 !end do

!close(222)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine SNAP_coulomb_energy

end module SNAP_class

