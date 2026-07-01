module force_field_min_class

  use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int,C_char
  use LAMMPS
  !use MPI
  use SNAP_class
  use target_functions_class
  use atoms_class
  use lammps_class
  use parameters
  use lapack_inverse
  !use rescale
  implicit none

  type, extends(target_function)   :: force_field_min
  type(lammps_obj)                 :: object_lammps
  type(lammps_obj),allocatable     :: set(:)
  type(lammps_obj),allocatable     :: set_AL(:)
  double precision                 :: R_screen
  double precision                 :: cutoff_dip
  double precision                 :: cutoff_en
  integer                          :: twojmax_dip,twojmax_en
  integer                          :: num_bisp_en,num_bisp_dip
  integer                          :: nconfig_AL
  integer,dimension(4)             :: iseed
  double precision                 :: timestep
  double precision                 :: temperature,temperature_in,temperature_final
  double precision                 :: lambda_dip    
  double precision                 :: weight
  double precision                 :: factor_thresh  
  double precision                 :: sigma_AL,thresh_AL
  double precision                 :: k_AL,C_M0
  double precision                                              :: lambda_en
  double precision                                          :: s_z
  double precision                 :: error
  double precision,allocatable                 :: error_forces(:)
  double precision,dimension(:),allocatable                 :: charge
  logical,dimension(:),allocatable :: coeff_mask_en,coeff_mask_dip
  double precision,dimension(:),allocatable :: tot_charge
  double precision, dimension(:,:),allocatable :: SNAP_prediction_matrix,SNAP_matrix_A
  character(len=120)               :: dipoles_file,geometry_file,energy_file,record_file
  logical                          :: VdW_flag, coul_flag
  logical                          :: minim_flag, md_flag, rampa_flag
  logical                          :: active_learn, shift_flag,post_AL
  logical                          :: flag_energy,flag_forces
  logical                          :: periodic_flag
  integer                          :: counter,nsteps
   double precision               :: eps=1.0e-7
  double precision               :: beta1=0.9d0
  double precision               :: beta2=0.999d0
   integer                        :: nval
  double precision               :: ener
  double precision, allocatable  :: val(:)
  double precision, allocatable  :: grad(:)
  double precision, allocatable  :: loc_lr(:)
  double precision               :: max_grad=0.01e5
  double precision               :: lr=0.0001d0
  integer                        :: max_iter=3000
  logical                        :: print_grad=.false.
  integer                        :: print_grad_io=22
  logical                        :: print_val=.false.
  integer                        :: print_val_io=23
  

  contains
  procedure                        :: get_fval  => get_energy
  procedure                        :: get_fgrad => calc_force
  procedure                        :: init_vel
  procedure                        :: propagate_md
  procedure                        :: control_structure
  procedure                        :: minimize_adam
  procedure                        :: get_prediction_err
  procedure                        :: shift_CM
  procedure                        :: update_temperature
 end type force_field_min

  contains

  subroutine get_energy(this,vec,val)
  implicit none
  class(force_field_min)             :: this
  double precision                   :: val
  double precision, allocatable      :: vec(:),beta(:),coul_energy(:)
  integer                            :: i,j,tot_kinds

   call get_tot_kinds(this%set,tot_kinds)
   allocate(beta(this%num_bisp_en*tot_kinds))

  open(222,file='snapcoeff_energy',action='read')
  do i=1,this%num_bisp_en*tot_kinds
  read(222,*) beta(i)
  end do
  close(222)

  val=0.0

  do i=1,this%object_lammps%nats
   do j=1,this%num_bisp_en

    val=val+beta((this%object_lammps%kind(i)-1)*this%num_bisp_en+j)*this%object_lammps%at_desc(i)%desc(j)
   
   end do
  end do

  deallocate(beta)

if (this%coul_flag) then

  call SNAP_coulomb_energy(this%set,this%R_screen,coul_energy)

  val=val+coul_energy(1)

end if

if (this%VdW_flag) then

        val=val+Har_to_kc*this%object_lammps%edisp

end if

if (this%shift_flag) then

        val=val+Har_to_kc*this%object_lammps%e_shift

end if


end subroutine get_energy

subroutine calc_force(this,vec,val,grad)
  implicit none
  class(force_field_min)         :: this
  double precision               :: val
  integer                        :: i,j
  double precision, allocatable  :: vec(:),grad(:)
  type(lammps_obj)               :: set_tmp
  double precision,allocatable   :: forces(:,:),forces_1(:,:),forces_model(:,:)
  character(len=120)             :: atom_string,set_newpos_string,cutoff_string
  real                           :: start,finish

do i=1,this%object_lammps%nats
 do j=1,3
    this%object_lammps%x(i,j)=vec((i-1)*3+j)
 end do
end do

  call lammps_open_no_mpi("lmp -screen none -log log.simple",this%object_lammps%lmp) 
  call this%object_lammps%setup_lammps_lattice(this%object_lammps%nkinds)     
  if (associated(this%object_lammps%at_desc)) then
  
  do i=1,this%object_lammps%nats
     if (allocated(this%object_lammps%at_desc(i)%desc)) deallocate(this%object_lammps%at_desc(i)%desc)
  end do
  
  deallocate(this%object_lammps%at_desc)
  end if
  !call this%object_lammps%get_bis(this%cutoff_dip,this%twojmax_dip)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!crucial point
 this%set(1)=this%object_lammps
 !call fit_dipoles(this%set,this%num_bisp_dip,trim(this%dipoles_file),this%lambda_dip,this%coeff_mask_dip,'VALID',&
  !     this%tot_charge,len_trim(this%dipoles_file))
   
!do i=1,this%object_lammps%nats
! if (allocated (this%object_lammps%at_desc(i)%desc)) deallocate(this%object_lammps%at_desc(i)%desc)
!end do

! deallocate(this%object_lammps%at_desc)

call this%object_lammps%get_bis(this%cutoff_en,this%twojmax_en)
!this%object_lammps%charge=this%set(1)%charge

! if (.NOT.this%coul_flag) then
!        this%object_lammps%charge=0.0d0
!end if

if (this%VdW_flag) then

 allocate(this%object_lammps%grads(3,this%object_lammps%nats))
 call this%object_lammps%grimme_d3(this%periodic_flag)

end if 

!if (this%shift_flag) then

!call this%object_lammps%shift_mol(this%k_AL,this%C_M0)

!end if

 call this%get_fval(vec,val)

set_tmp=this%object_lammps

  call lammps_open_no_mpi("lmp -screen none -log log.tmp",set_tmp%lmp)
  call set_tmp%setup_lammps_lattice(set_tmp%nkinds)
 
  !do i=1,this%object_lammps%nats

   !write(atom_string,*) 'set atom',i,'type',i
   !call lammps_command(set_tmp%lmp,trim(atom_string))

  !end do
  
  call set_tmp%get_der_bis(this%cutoff_en,this%twojmax_en,this%object_lammps%nkinds)

  if (.not. allocated(this%object_lammps%bisp_der)) allocate &
  (this%object_lammps%bisp_der(size(set_tmp%bisp_der,1),size(set_tmp%bisp_der,2)))
                this%object_lammps%bisp_der=set_tmp%bisp_der/A_to_B

!if (.not.this%coul_flag) then
!        this%object_lammps%charge=0.0d0
!end if



   call this%object_lammps%get_forces(this%R_screen,forces,'ENERGY',this%num_bisp_en,this%num_bisp_dip)

  !call set_tmp%get_der_bis(this%cutoff_dip,this%twojmax_dip,this%object_lammps%nats)
  !      this%object_lammps%bisp_der=set_tmp%bisp_der/A_to_B

!if (.NOT.this%coul_flag) then
!        this%object_lammps%charge=0.0d0
!end if
  

!call this%object_lammps%get_forces(this%R_screen,forces_1,'DIPOLE',this%num_bisp_en,this%num_bisp_dip)


deallocate(set_tmp%bisp_der,this%object_lammps%bisp_der)


if (.not.allocated(forces_model)) allocate(forces_model(this%object_lammps%nats,3))

 !forces_model=(forces_1+forces)*F_conv
 forces_model=forces*F_conv

if (this%VdW_flag) then
 
 forces_model=forces_model-transpose(this%object_lammps%grads)*F_conv
 
 deallocate(this%object_lammps%grads)

end if

!if (this%shift_flag) then

!       forces_model=forces_model + this%object_lammps%force_shift*F_conv

!end if

        call lammps_close(set_tmp%lmp)
        call lammps_close(this%object_lammps%lmp)

do i=1,this%object_lammps%nats 
if (allocated (this%object_lammps%at_desc(i)%desc)) deallocate(this%object_lammps%at_desc(i)%desc)
end do 
deallocate(this%object_lammps%at_desc)

 this%object_lammps%at_desc=> NULL()


if (.not.(allocated(grad))) allocate(grad(this%object_lammps%nats*3))

      
  do i=1,this%object_lammps%nats
   do j=1,3

    grad(((i-1)*3)+j)=-forces_model(i,j)

   end do
  end do

 if (this%minim_flag) then

   open(111, file="traj_min_molforge.xyz", action="write",position='append')

  write(111,*)this%object_lammps%nats
  write(111,*)'XXX'

 do i=1,this%object_lammps%nats
  write(111,*) this%object_lammps%label(this%object_lammps%kind(i)),vec(((i-1)*3)+1:(i*3))
  end do
 close(111)

end if


end subroutine calc_force

subroutine init_vel(this,temperature)
  implicit none
  class(force_field_min),intent(inout)         :: this
  double precision, intent(in)                 :: temperature
  integer                                      :: i
  integer                                      :: idist
  integer                                      :: N
  double precision, allocatable                :: X(:)
  double precision, dimension(3)               :: lin_mom_tot
  double precision                             :: Y
  double precision                             :: E_kin,T_system

  if (.not. allocated (this%object_lammps%v)) then
   allocate(this%object_lammps%v(this%object_lammps%nats,3))
  end if
  
  idist=3
  N=3
  allocate(X(N))

    do i=1,this%object_lammps%nats
  call dlarnv(idist,this%iseed,N,X)

  this%object_lammps%v(i,:)=dsqrt(boltz*temperature/(1822.89d0*this%object_lammps%mass(this%object_lammps%kind(i))))*X
  end do

  deallocate(X)

E_kin=0.0d0
lin_mom_tot=0.0d0

do i=1,this%object_lammps%nats

  lin_mom_tot=lin_mom_tot + this%object_lammps%v(i,:)*(amu_to_emass*this%object_lammps%mass(this%object_lammps%kind(i)))

end do

do i=1,this%object_lammps%nats

 this%object_lammps%v(i,:)=this%object_lammps%v(i,:)-(lin_mom_tot/(dble(this%object_lammps%nats)* &
        amu_to_emass* this%object_lammps%mass(this%object_lammps%kind(i))))

end do


do i=1,this%object_lammps%nats

  E_kin=E_kin+0.5d0*sum((this%object_lammps%v(i,:)**2))*(amu_to_emass*this%object_lammps%mass(this%object_lammps%kind(i)))

end do

  

  T_system=(2.0d0*E_kin)/((3.0d0*this%object_lammps%nats-3.0d0)*boltz)
  this%object_lammps%v(:,:)=dsqrt(temperature/T_system)*this%object_lammps%v(:,:)
 

end subroutine init_vel

 subroutine minimize_adam(this,max_iter,start_iter)
         implicit none
         class(force_field_min)        :: this
         integer                       :: iter,i,iter0,j
         integer, optional             :: max_iter,start_iter
         double precision              :: gradnorm,val
         double precision              :: E1,E2
         double precision, allocatable :: gradres(:),gradres2(:),vec(:)
         integer(kind=4)               ::  t1,rate,t2
         logical                       :: flag_new_struct

          if (.not.allocated (vec)) allocate(vec(this%object_lammps%nats*3))
          do i=1,this%object_lammps%nats
          do j=1,3
                vec((i-1)*3+j)=this%object_lammps%x(i,j)
          end do
          end do

          iter0=0

          if(present(start_iter)) iter0=start_iter
          if(present(max_iter)) this%max_iter=max_iter
          if(this%print_grad) open(this%print_grad_io,file='grad.dat')
          if(this%print_val) open(this%print_val_io,file='param.dat')

          if(allocated(gradres)) deallocate(gradres)
          if(allocated(gradres2)) deallocate(gradres2)
          allocate(gradres(size(vec)))
          allocate(gradres2(size(vec)))
          gradres=0.0d0
          gradres2=0.0d0
          iter=1
           
          E2=0.0
          call system_clock(t1,rate)
          do while (iter.le.this%max_iter)
           
           E1=E2
           call this%get_fgrad(vec,val,this%grad)
           E2=val

           open(222,file='grad_en_molforge.txt',action='write',position='append')
           write(222,*) iter+iter0,sqrt(gradnorm/size(this%grad)),val,&
                        maxval(this%grad),abs(E2-E1)
           close(222)
           
           
           if(this%print_grad) write(this%print_grad_io,*) this%grad
           gradnorm=0.0d0
           do i=1,size(this%grad)
           gradres(i)=this%beta1*gradres(i)+(1-this%beta1)*this%grad(i)
           gradres2(i)=this%beta2*gradres2(i)+(1-this%beta2)*this%grad(i)**2
           gradnorm=gradnorm+this%grad(i)**2
           enddo

           if(this%print_val) write(this%print_val_io,*) vec
           if(this%print_grad) write(this%print_grad_io,*) this%grad

           if(allocated(this%loc_lr))then
           vec=vec-this%lr*(gradres/(1-this%beta1**(iter+1)))&
                 /(sqrt(gradres2/(1-this%beta2**(iter+1)))+this%eps)*this%loc_lr
           else
            vec=vec-this%lr*(gradres/(1-this%beta1**(iter+1)))&
                 /(sqrt(gradres2/(1-this%beta2**(iter+1)))+this%eps)
           endif

           if (this%active_learn) then

           call this%get_prediction_err(vec,this%thresh_AL,iter,flag_new_struct,this%flag_energy,&
           this%flag_forces)
           if (flag_new_struct) then
           !call this%get_fgrad(vec,val,force)
           call system_clock(t2)
           open(222,file="AL_stats.txt",action="write",position="append")
           write(222,*) iter,(real(t2-t1)/real(rate))
           write(222,*)"GO"
           close(222)
           stop
           end if

           end if
         
          !!!!!!!!!!!!!!!!CHECK CONVERGENCE CRITERIA
          if ((sqrt(gradnorm/size(this%grad))<0.01).and.(maxval(this%grad)<0.1)&
           .and.(abs(E2-E1)<0.0001)) then
           call system_clock(t2)
           open(222,file="AL_stats.txt",action="write",position="append")
           write(222,*) "It converged after",iter,"iterations in",(real(t2-t1)/real(rate)),"seconds"
           close(222)

          stop
          end if

          iter=iter+1
          enddo

          if(this%print_grad) close(this%print_grad_io)
          if(this%print_val) close(this%print_val_io)
          if(allocated(gradres)) deallocate(gradres)
          if(allocated(gradres2)) deallocate(gradres2)

         return
         end subroutine minimize_adam



subroutine propagate_md(this,nsteps,dt)
implicit none
class(force_field_min),intent(inout)         :: this
double precision, allocatable                :: acc(:)
double precision, allocatable                :: force(:)
integer,intent(in)                           :: nsteps
double precision, intent(in)                 :: dt
double precision,allocatable                 :: vec(:),vel(:)
double precision                             :: val,E_kin,E_kin_new
double precision                             :: temp
double precision                             :: charge_1,charge_2
double precision                             :: start,finish,start_loop,finish_loop
double precision                             :: T_system
double precision                             :: resamplekin
integer                                      :: i,j
integer                                      :: iter
integer                                      :: call_active
logical                                      :: flag_new_struct
integer(kind=4)                              :: t1,rate,t2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (.not.allocated (vec)) allocate(vec(this%object_lammps%nats*3))
if (.not.allocated (vel)) allocate(vel(this%object_lammps%nats*3))

iter=0
this%counter=0

do i=1,this%object_lammps%nats
 do j=1,3

  vec((i-1)*3+j)=this%object_lammps%x(i,j)
  vel((i-1)*3+j)=this%object_lammps%v(i,j)

 end do
end do


open(111, file="traj_MD_molforge.xyz", action="write",position='append')

  write(111,*)this%object_lammps%nats
  write(111,*)'XXX'

 do i=1,this%object_lammps%nats
 write(111,*) this%object_lammps%label(this%object_lammps%kind(i)),vec(((i-1)*3)+1:(i*3)) 
 end do
  close(111)

if (.not.allocated(force)) allocate(force(this%object_lammps%nats*3))

force=0.0d0

if (this%shift_flag) then
this%C_M0=this%shift_CM(dble(iter))
end if

call this%get_fgrad(vec,val,force)

E_kin=0.0d0

do i=1,this%object_lammps%nats
 do j=1,3

  E_kin=E_kin+0.5d0*(vel((i-1)*3+j)**2)*amu_to_emass*this%object_lammps%mass(this%object_lammps%kind(i))

 end do
end do

temp=(2.0d0*E_kin)/((3.0d0*this%object_lammps%nats-3.0d0)*boltz)

open(111, file="etotal_kin_pot_temp_molforge.txt", action="write",position='append')

write(111,*) '# Total energy      Kinetic Energy   Potential Energy   Temperature'
write(111,*) E_kin*Har_to_kc+val, E_kin*Har_to_kc,val,temp

close(111)

if (.not.allocated(acc)) allocate(acc(this%object_lammps%nats*3))

acc=0.0d0

do i=1,this%object_lammps%nats
 do j=1,3
 
  acc((i-1)*3+j)=-(force((i-1)*3+j)/F_conv)/(1822.89d0*this%object_lammps%mass(this%object_lammps%kind(i)))

  end do
 end do

do while (iter < nsteps)
iter=iter +1

force=0.0

this%counter=this%counter+1

!E_kin=0.0

!do i=1,this%object_lammps%nats
! do j=1,3

!  E_kin=E_kin+0.5d0*(vel((i-1)*3+j)**2)*1822.89d0*this%object_lammps%mass(this%object_lammps%kind(i))

! end do
!end do

!temp=(2.0d0*E_kin)/(3.0d0*this%object_lammps%nats*boltz)

!if ((temp < this%temperature-50.0d0).or.(temp > this%temperature + 50.0d0)) then
!  vel=dsqrt(this%temperature/temp)*vel
!end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!calculate half-step velocities

vel=vel+0.5d0*acc*dt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!calculate full step positions

vec=(vec*A_to_B)+vel*dt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ACTIVE LEARNING CALL

vec=vec/A_to_B

if (this%active_learn) then
call this%get_prediction_err(vec,this%thresh_AL,iter,flag_new_struct,this%flag_energy,&
        this%flag_forces)

if (iter==this%nsteps) then

open(222,file="AL_stats.txt", action="write", position="append")

 write(222,*) iter,val,real(t2-t1)/real(rate)

close(222)
        
open(222,file=trim(this%record_file), action="write", position="append")

 write(222,*) "1"

close(222)

end if 

if (flag_new_struct) then
 call this%get_fgrad(vec,val,force)

open(222,file="AL_stats.txt",action="write",position="append")

 write(222,*) iter,val,real(t2-t1)/real(rate)
 write(222,*) "GO"

close(222)

 stop

 end if
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (this%shift_flag) then
this%C_M0=this%shift_CM(dble(iter))
end if

call this%get_fgrad(vec,val,force)


do i=1,this%object_lammps%nats
 do j=1,3

  acc((i-1)*3+j)=-(force((i-1)*3+j)/F_conv)/(1822.89d0*this%object_lammps%mass(this%object_lammps%kind(i)))

 end do
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!calculate full step velocities

vel=vel+0.5d0*acc*dt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
E_kin=0.0

do i=1,this%object_lammps%nats
 do j=1,3

  E_kin=E_kin+0.5d0*(vel((i-1)*3+j)**2)*1822.89d0*this%object_lammps%mass(this%object_lammps%kind(i))

 end do
end do

temp=(2.0d0*E_kin)/(3.0d0*this%object_lammps%nats*boltz)

if (this%rampa_flag) then
call this%update_temperature(iter)
end if

E_kin_new= resamplekin(E_kin,((3.0d0*this%object_lammps%nats-3.d0)/2.0d0)*this%temperature_final &
        *boltz,(3*this%object_lammps%nats)-3,100.0d0)


vel=dsqrt(E_kin_new/E_kin)*vel

E_kin=0.0

do i=1,this%object_lammps%nats
 do j=1,3

  E_kin=E_kin+0.5d0*(vel((i-1)*3+j)**2)*1822.89d0*this%object_lammps%mass(this%object_lammps%kind(i))

 end do
end do

temp=(2.0d0*E_kin)/((3.0d0*this%object_lammps%nats-3.0d0)*boltz)


if (mod(iter,5)==0) then

if (this%active_learn) then

open(111, file="etotal_kin_pot_temp_molforge.txt", action="write",position="append")
  write(111,*) E_kin*Har_to_kc+val, E_kin*Har_to_kc,val,temp,this%error
close(111)

else

open(111, file="etotal_kin_pot_temp_molforge.txt", action="write",position="append")
  write(111,*) E_kin*Har_to_kc+val, E_kin*Har_to_kc,val,temp
close(111)

end if

open(111, file="traj_MD_molforge.xyz", action="write",position='append')

 write(111,*)this%object_lammps%nats
  write(111,*)'Iter',iter

 do i=1,this%object_lammps%nats
 write(111,*) this%object_lammps%label(this%object_lammps%kind(i)),vec(((i-1)*3)+1:(i*3)) 
 end do
  close(111)


end if

end do

if (allocated(this%object_lammps%grads)) deallocate(this%object_lammps%grads)

if (associated(this%object_lammps%at_desc)) then
 do i=1,this%object_lammps%nats
  if (allocated (this%object_lammps%at_desc(i)%desc)) deallocate(this%object_lammps%at_desc(i)%desc)
 end do
 deallocate(this%object_lammps%at_desc)
 this%object_lammps%at_desc=> NULL()
end if

end subroutine propagate_md

subroutine control_structure(this,vec,sigma,thresh,iter)
class(force_field_min),intent(inout)                    :: this
double precision,allocatable,intent(in)                 :: vec(:)
double precision,intent(in)                             :: sigma
double precision,intent(in)                             :: thresh
double precision                                        :: dist_max,dist_tmp
integer                                                 :: i,j,k
integer,intent(in)                                      :: iter

if (iter==1) then 

do i=1,size(this%set_AL)

 call lammps_open_no_mpi("lmp -screen none -log log.ebete",this%set_AL(i)%lmp)

 call this%set_AL(i)%setup_lammps_lattice(this%set_AL(i)%nkinds)

 call this%set_AL(i)%get_bis(this%cutoff_en,this%twojmax_en)

 call lammps_close(this%set_AL(i)%lmp)

end do

end if

this%object_lammps%x=transpose(reshape(vec,(/3,this%object_lammps%nats/)))

call lammps_open_no_mpi("lmp -screen none -log log.simple",this%object_lammps%lmp)

call this%object_lammps%setup_lammps_lattice(this%object_lammps%nkinds)

call this%object_lammps%get_bis(this%cutoff_en,this%twojmax_en)

call lammps_close(this%object_lammps%lmp)


do i=1,this%object_lammps%nats

dist_tmp=0.0d0
dist_max=0.0d0
 

 do j=1,size(this%set_AL)
  do k=1,this%set_AL(j)%nats
   if (this%set_AL(j)%kind(k).eq.this%object_lammps%kind(i)) then

    dist_tmp=dexp(-sum((this%set_AL(j)%at_desc(k)%desc(:)-this%object_lammps%at_desc(i)%desc(:))**2)/((sigma)**2))
    
     if (dist_tmp >= dist_max) then
          dist_max=dist_tmp
     end if

   end if
 
 end do

end do

if (dist_max < thresh) then

open(111, file="new_geo_AL.xyz", action="write")

  write(111,*)this%object_lammps%nats
  write(111,*)'XXX'

  do j=1,this%object_lammps%nats
   write(111,*) this%object_lammps%label(this%object_lammps%kind(j)),vec(((j-1)*3)+1:(j*3))
  end do
 
close(111)
   
!open(111, file=trim(this%geometry_file), action="write",position='append')
!     write(111,*)this%object_lammps%nats
!  write(111,*)this%object_lammps%cell(1,:),this%object_lammps%cell(2,:),this%object_lammps%cell(3,:),&
!                     this%object_lammps%nkinds

!  do j=1,this%object_lammps%nats
!   write(111,*) this%object_lammps%label(this%object_lammps%kind(j)),vec(((j-1)*3)+1:(j*3)), &
!           this%object_lammps%kind(j),this%object_lammps%mass(this%object_lammps%kind(j))
!  end do


!close(111)

write(*,*) 'New structure found after',iter,'steps of MD'

stop
 
end if

end do

do i=1,this%object_lammps%nats
if (allocated (this%object_lammps%at_desc(i)%desc)) deallocate(this%object_lammps%at_desc(i)%desc)
end do
deallocate(this%object_lammps%at_desc)

end subroutine control_structure

subroutine get_prediction_err(this,vec,thresh,iter,flag_new_struct,flag_energy,flag_forces)
class(force_field_min),intent(inout)                    :: this
double precision,intent(inout)                          :: thresh
double precision,allocatable,intent(in)                 :: vec(:)
integer,intent(in)                                      :: iter
logical, intent(out)                                    :: flag_new_struct
logical, intent(in)                                     :: flag_energy,flag_forces
double precision                                        :: tmp,error,ave_atom
double precision,allocatable                            :: tmp_forces(:)
integer                                                 :: tot_kinds,N,tot_atom
integer                                                 :: size_ref,start_snap_force
integer                                                 :: counter,ez_cons_rows,end_loop
integer                                                 :: start_loop,help_counter
double precision,dimension(:,:),allocatable             :: X,X_t,A,D,D_t,test,test_copy
double precision,dimension(:),allocatable               :: theta,y,y_pred,AP
double precision,allocatable                           :: x_new(:,:),K_mat(:,:)
double precision,dimension(:),allocatable               :: mu,sigma
integer                                                 :: i,j,k,l
integer                                                 :: comp
double precision                                        :: NUMERATOR,DENOMINATOR, s_z
character(len=30)                                        :: filename

!!!!!!!!!!!!!!!!!!!!!!!!!!!

flag_new_struct=.false.
call get_tot_kinds(this%set_AL,tot_kinds)
call get_ave_atoms(this%set_AL,ave_atom,tot_atom)

N=this%num_bisp_en*tot_kinds

if (iter==1) then

do i=1,size(this%set_AL)

 call lammps_open_no_mpi("lmp -screen none -log log.ebete",this%set_AL(i)%lmp)

 call this%set_AL(i)%setup_lammps_lattice(this%set_AL(i)%nkinds)

if (this%flag_energy) then

        call this%set_AL(i)%get_bis(this%cutoff_en,this%twojmax_en)

end if

if (this%flag_forces) then

       call this%set_AL(i)%get_der_bis(this%cutoff_en,this%twojmax_en,this%set_AL(i)%nkinds)

end if

 call lammps_close(this%set_AL(i)%lmp)

end do


ez_cons_rows = count(this%coeff_mask_en)

if ((flag_energy).and.(.not.flag_forces)) then

  size_ref=size(this%set_AL) + ez_cons_rows
  start_loop=size(this%set_AL)+1
  end_loop=size(this%set_AL)
  help_counter=size(this%set_AL)

else if ((flag_forces).and.(.not.flag_energy)) then

  size_ref=3*tot_atom + ez_cons_rows
  start_loop=3*tot_atom+1
  end_loop=3*tot_atom
  help_counter=3*tot_atom

else if ((flag_forces).and.(flag_energy)) then

  size_ref=size(this%set_AL)+3*tot_atom+ez_cons_rows
  start_loop=size(this%set_AL)+3*tot_atom+1
  end_loop=size(this%set_AL)+3*tot_atom 
  help_counter=size(this%set_AL)+3*tot_atom

end if

if (this%lambda_en.ne.0.0) then

  size_ref=size_ref+(this%num_bisp_en-1)*tot_kinds

end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(A(size_ref,N))
allocate(this%SNAP_matrix_A(size_ref,N))
A=0.0
allocate(y(size_ref),y_pred(size_ref))
y=0.0

open(100,file="fit_energy_rms.dat")

do j=1,end_loop

  read(100,*) y(j)

 end do

close(100)

y(1:size(this%set_AL))=y(1:size(this%set_AL))*this%weight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (this%lambda_en.ne.0.0) then

counter=0

 do j=start_loop,size_ref-ez_cons_rows

  if (mod(j-help_counter,this%num_bisp_en-1)==1) then
    counter=counter+1
  end if

   A(j,j-help_counter+counter) = sqrt(this%lambda_en)

 end do

end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!BUILDING THE SNAP MATRIX

if ((this%flag_energy)) then

do j=1,size(this%set_AL)

 do i=1,this%set_AL(j)%nats

  do k=1,this%num_bisp_en

A(j,(this%set_AL(j)%kind(i)-1)*this%num_bisp_en+k) = A(j,(this%set_AL(j)%kind(i)-1)*this%num_bisp_en+k)&
        + this%set_AL(j)%at_desc(i)%desc(k)

  end do

 end do

end do

end if

if ((this%flag_forces).and.(.not.this%flag_energy)) then
        start_snap_force=0
else if ((this%flag_forces).and.(this%flag_energy)) then
        start_snap_force=size(this%set_AL)
end if

if (this%flag_forces) then
        do j=1,size(this%set_AL)

  do i=1,this%set_AL(j)%nkinds

   do l=1,this%set_AL(j)%nats

    do comp=1,3

     do k=1,this%num_bisp_en-1

        A(start_snap_force+(j-1)*this%set_AL(j)%nats*3+(l-1)*3+comp,this%num_bisp_en*(i-1)+1+k)=&
                this%set_AL(j)%bisp_der((i-1)*((this%num_bisp_en-1)*3)+(this%num_bisp_en-1)*(comp-1)+k,l)

     end do

    end do

   end do

  end do

 end do

end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FLAG E0 TO REGULARIZE THE FIRST COEFFICIENTS OF BISPECTRUM

counter=0

do j=1,tot_kinds

 if (this%coeff_mask_en(j)) then
  
  counter=counter+1
  A(size_ref - ez_cons_rows + counter,(j-1)*this%num_bisp_en + 1)= 1

 end if

end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(D(size_ref,N))
allocate(D_t(N,size_ref))
allocate(test(N,N))
allocate(this%SNAP_prediction_matrix(N,N))

D=A
this%SNAP_matrix_A=A

D(1:size(this%set_AL),:)=D(1:size(this%set_AL),:)*this%weight



D_t=transpose(D)
test=matmul(D_t,D)

call mat_inv(test,N)
this%SNAP_prediction_matrix=test

deallocate(test)

if (.not.allocated(theta)) allocate(theta(N))
theta=0.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        open(100,file="snapcoeff_energy")

        do j=1,N

         read(100,*) theta(j)

        end do

        close(100)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

y_pred=matmul(D,theta)

NUMERATOR=0.0d0

do j=1,size_ref

 NUMERATOR=NUMERATOR+(y_pred(j)-y(j))**2

end do

DENOMINATOR=dble(size_ref-N-1)

this%s_z=dsqrt(NUMERATOR/DENOMINATOR)

open(113,file='fattore_di_errore',action='write',position='append')
write(113,*)this%s_z 
close(113)
end if


this%object_lammps%x=transpose(reshape(vec,(/3,this%object_lammps%nats/)))

call lammps_open_no_mpi("lmp -screen none -log log.simple",this%object_lammps%lmp)

call this%object_lammps%setup_lammps_lattice(this%object_lammps%nkinds)

if (this%flag_energy) then

 call this%object_lammps%get_bis(this%cutoff_en,this%twojmax_en)

end if

if (this%flag_forces) then

  call this%object_lammps%get_der_bis(this%cutoff_en,this%twojmax_en,this%object_lammps%nkinds)

end if

call lammps_close(this%object_lammps%lmp)

if ((this%flag_energy).and.(.not.this%flag_forces)) then

if (.not.allocated(x_new)) allocate(x_new(N,1))
x_new=0.0d0

do i=1,this%object_lammps%nats
 do k=1,this%num_bisp_en

   x_new((this%object_lammps%kind(i)-1)*this%num_bisp_en+k,1) = x_new((this%object_lammps%kind(i)-1)*this%num_bisp_en+k,1)&
           + this%object_lammps%at_desc(i)%desc(k)

 end do
end do

else if ((this%flag_forces).and.(.not.this%flag_energy)) then

if (.not.allocated(x_new)) allocate(x_new(N,this%object_lammps%nats*3))
x_new=0.0d0

do i=1,this%object_lammps%nkinds
 do l=1,this%object_lammps%nats
  do comp=1,3
   do k=1,this%num_bisp_en-1

   x_new((i-1)*(this%num_bisp_en)+k+1,(l-1)*3+comp) = this%object_lammps&
           %bisp_der((i-1)*((this%num_bisp_en-1)*3)+(this%num_bisp_en-1)*(comp-1)+k,l)

   end do
  end do
 end do
end do

else if ((this%flag_forces).and.(this%flag_energy)) then

if (.not.allocated(x_new)) allocate(x_new(N,1+this%object_lammps%nats*3))
x_new=0.0d0
do i=1,this%object_lammps%nats
 do k=1,this%num_bisp_en

   x_new((this%object_lammps%kind(i)-1)*this%num_bisp_en+k,1) = x_new((this%object_lammps%kind(i)-1)*this%num_bisp_en+k,1)&
           + this%object_lammps%at_desc(i)%desc(k)

 end do
end do

do i=1,this%object_lammps%nkinds
 do l=1,this%object_lammps%nats
  do comp=1,3
   do k=1,this%num_bisp_en-1

   x_new((i-1)*(this%num_bisp_en)+k+1,(l-1)*3+comp+1) = this%object_lammps&
           %bisp_der((i-1)*((this%num_bisp_en-1)*3)+(this%num_bisp_en-1)*(comp-1)+k,l)

   end do
  end do
 end do
end do

end if
 
if ((this%flag_energy).and.(.not.this%flag_forces)) then

        if (.not.allocated(K_mat)) allocate(K_mat(N,1))

else if ((this%flag_forces).and.(.not.this%flag_energy)) then

 if (.not.allocated(K_mat)) allocate (K_mat(N,this%object_lammps%nats*3))
 if (.not.allocated(tmp_forces)) allocate (tmp_forces(this%object_lammps%nats*3))

else if ((this%flag_forces).and.(this%flag_energy)) then

 if (.not.allocated(K_mat)) allocate (K_mat(N,this%object_lammps%nats*3+1))
 if (.not.allocated(tmp_forces)) allocate (tmp_forces(this%object_lammps%nats*3+1))

end if 

K_mat=0.0d0

K_mat=matmul(this%SNAP_prediction_matrix,x_new)

tmp=0.0d0

if ((this%flag_energy).and.(.not.this%flag_forces)) then

do i=1,this%num_bisp_en*tot_kinds

   tmp=tmp+K_mat(i,1)*x_new(i,1)

end do

else if ((this%flag_forces).and.(.not.this%flag_energy)) then
tmp_forces=0.0d0

do j=1,this%object_lammps%nats*3
 do i=1,N
         tmp_forces(j)=tmp_forces(j)+K_mat(i,j)*x_new(i,j)
 end do
end do

else if ((this%flag_forces).and.(this%flag_energy)) then

tmp_forces=0.0d0

do j=1,this%object_lammps%nats*3+1
 do i=1,N
         tmp_forces(j)=tmp_forces(j)+K_mat(i,j)*x_new(i,j)
 end do
end do


end if

if ((this%flag_energy).and.(.not.this%flag_forces)) then

this%error=this%s_z*dsqrt(1+tmp)

else if ((this%flag_forces).and.(.not.this%flag_energy)) then

if (.not.allocated(this%error_forces)) allocate(this%error_forces(3*this%object_lammps%nats))

 this%error_forces=this%s_z*dsqrt(1+tmp_forces)
 this%error=maxval(this%error_forces)

else if ((this%flag_forces).and.(this%flag_energy)) then

if (.not.allocated(this%error_forces)) allocate(this%error_forces(3*this%object_lammps%nats+1))

 this%error_forces=this%s_z*dsqrt(1+tmp_forces)
 this%error=maxval(this%error_forces)
end if

        thresh=this%factor_thresh*this%s_z

if (this%error>thresh) then

        write(*,*) "Configuration found after",iter,"steps"

flag_new_struct=.true.

 open(221,file="errors_predicted",action="write",position="append")
  write(221,*) this%error
 close(221)

open(111, file="new_geo_AL.xyz", action="write")

  write(111,*)this%object_lammps%nats
  write(111,*)'XXX'

  do j=1,this%object_lammps%nats
   write(111,*) this%object_lammps%label(this%object_lammps%kind(j)),vec(((j-1)*3)+1:(j*3))
  end do

close(111)

!open(111, file=trim(this%geometry_file), action="write",position='append')
!     write(111,*)this%object_lammps%nats
!  write(111,*)this%object_lammps%cell(1,:),this%object_lammps%cell(2,:),this%object_lammps%cell(3,:),&
!                     this%object_lammps%nkinds

!  do j=1,this%object_lammps%nats
!   write(111,*) this%object_lammps%label(this%object_lammps%kind(j)),vec(((j-1)*3)+1:(j*3)), &
!           this%object_lammps%kind(j),this%object_lammps%mass(this%object_lammps%kind(j))
!  end do


!close(111)


end if

deallocate(K_mat,x_new)

end subroutine get_prediction_err

function shift_CM (this,t)
implicit none
double precision, intent(in)          :: t
double precision                      :: shift_CM
class(force_field_min),intent(inout)  :: this

        if (t<10000) then
        shift_CM = (2.9d0*A_to_B) + ((1.1d0*A_to_B)/(1000.0d0*41.49d0))*(t*this%timestep)
else
        shift_CM = (4.0d0*A_to_B) + (1.1d0*A_to_B)*dsin((t-10000)*this%timestep*0.000015136d0) ! 10 ps period
end if
return

end function shift_CM

subroutine update_temperature (this,iter)
class(force_field_min),intent(inout)  :: this
integer, intent(in)                   :: iter

if (iter<20000) then
        this%temperature=this%temperature_in+((this%temperature_final-this%temperature_in)/20000.0d0)*dble(iter)
else
        this%temperature=this%temperature_final
end if

end subroutine update_temperature

end module force_field_min_class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Stochastic velocity rescale, as described in
! Bussi, Donadio and Parrinello, J. Chem. Phys. 126, 014101 (2007)
!
! This subroutine implements Eq.(A7) and returns the new value for the kinetic energy,
! which can be used to rescale the velocities.
! The procedure can be applied to all atoms or to smaller groups.
! If it is applied to intersecting groups in sequence, the kinetic energy
! that is given as an input (kk) has to be up-to-date with respect to the previous rescalings.
!
! When applied to the entire system, and when performing standard molecular dynamics (fixed c.o.m. (center of mass))
! the degrees of freedom of the c.o.m. have to be discarded in the calculation of ndeg,
! and the c.o.m. momentum HAS TO BE SET TO ZERO.
! When applied to subgroups, one can chose to:
! (a) calculate the subgroup kinetic energy in the usual reference frame, and count the c.o.m. in ndeg
! (b) calculate the subgroup kinetic energy with respect to its c.o.m. motion, discard the c.o.m. in ndeg
!     and apply the rescale factor with respect to the subgroup c.o.m. velocity.
! They should be almost equivalent.
! If the subgroups are expected to move one respect to the other, the choice (b) should be better.
!
! If a null relaxation time is required (taut=0.0), the procedure reduces to an istantaneous
! randomization of the kinetic energy, as described in paragraph IIA.
!
! HOW TO CALCULATE THE EFFECTIVE-ENERGY DRIFT
! The effective-energy (htilde) drift can be used to check the integrator against discretization errors.
! The easiest recipe is:
! htilde = h + conint
! where h is the total energy (kinetic + potential)
! and conint is a quantity accumulated along the trajectory as minus the sum of all the increments of kinetic
! energy due to the thermostat.
!
!module rescale
!implicit none
!contains

function resamplekin(kk,sigma,ndeg,taut)
  implicit none
  double precision              :: resamplekin
  double precision,  intent(in)  :: kk    ! present value of the kinetic energy of the atoms to be thermalized (in arbitrary units)
  double precision,  intent(in)  :: sigma ! target average value of the kinetic energy (ndeg k_b T/2)  (in the same units as kk)
  integer, intent(in)  :: ndeg  ! number of degrees of freedom of the atoms to be thermalized
  double precision,  intent(in)  :: taut  ! relaxation time of the thermostat, in units of 'how often this routine is called'
  double precision :: factor,rr
  double precision, external :: gasdev

  if(taut>0.1) then
    factor=exp(-1.0/taut)
  else
    factor=0.0
  end if
  rr = gasdev()
  resamplekin = kk + (1.0-factor)* (sigma*(sumnoises(ndeg-1)+rr**2)/ndeg-kk) &
               + 2.0*rr*sqrt(kk*sigma/ndeg*(1.0-factor)*factor)

contains 

double precision function sumnoises(nn)
  implicit none
  integer, intent(in) :: nn
! returns the sum of n independent gaussian noises squared
! (i.e. equivalent to summing the square of the return values of nn calls to gasdev)
double precision, external :: gamdev,gasdev
  if(nn==0) then
    sumnoises=0.0
  else if(nn==1) then
    sumnoises=gasdev()**2
  else if(modulo(nn,2)==0) then
    sumnoises=2.0*gamdev(nn/2)
  else
    sumnoises=2.0*gamdev((nn-1)/2) + gasdev()**2
  end if
end function sumnoises

end function resamplekin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THE FOLLOWING ROUTINES ARE TRANSCRIBED FROM NUMERICAL RECIPES

double precision function gamdev(ia)
! gamma-distributed random number, implemented as described in numerical recipes

implicit none
integer, intent(in) :: ia
integer j
double precision am,e,s,v1,v2,x,y
double precision, external :: ran1
if(ia.lt.1)stop 'bad argument in gamdev'
if(ia.lt.6)then
  x=1.
  do 11 j=1,ia
    x=x*ran1()
11  continue
  x=-log(x)
else
1 v1=2.*ran1()-1.
    v2=2.*ran1()-1.
  if(v1**2+v2**2.gt.1.)goto 1
    y=v2/v1
    am=ia-1
    s=sqrt(2.*am+1.)
    x=s*y+am
  if(x.le.0.)goto 1
    e=(1.+y**2)*exp(am*log(x/am)-s*y)
  if(ran1().gt.e)goto 1
endif
gamdev=x
end function gamdev

double precision function gasdev()
! gaussian-distributed random number, implemented as described in numerical recipes

implicit none
integer, save :: iset = 0
double precision, save :: gset
double precision, external :: ran1
double precision fac,rsq,v1,v2
if(iset==0) then
1       v1=2.*ran1()-1.0d0
  v2=2.*ran1()-1.0d0
  rsq=v1**2+v2**2
  if(rsq.ge.1..or.rsq.eq.0.)goto 1
  fac=sqrt(-2.*log(rsq)/rsq)
  gset=v1*fac
  gasdev=v2*fac
  iset=1
else
  gasdev=gset
  iset=0
end if
end function gasdev

FUNCTION ran1()
! random number generator
INTEGER IA,IM,IQ,IR,NTAB,NDIV
double precision ran1,AM,EPS,RNMX
PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
  NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
INTEGER j,k,iv(NTAB),iy
SAVE iv,iy
DATA iv /NTAB*0/, iy /0/
INTEGER, SAVE :: idum=0 !! ATTENTION: THE SEED IS HARDCODED
if (idum.le.0.or.iy.eq.0) then
  idum=max(-idum,1)
  do 11 j=NTAB+8,1,-1
    k=idum/IQ
    idum=IA*(idum-k*IQ)-IR*k
    if (idum.lt.0) idum=idum+IM
    if (j.le.NTAB) iv(j)=idum
11  continue
  iy=iv(1)
endif
k=idum/IQ
idum=IA*(idum-k*IQ)-IR*k
if (idum.lt.0) idum=idum+IM
j=1+iy/NDIV
iy=iv(j)
iv(j)=idum
ran1=min(AM*iy,RNMX)
return
      END function ran1

      !end module rescale
