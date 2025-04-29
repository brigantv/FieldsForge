  module lammps_class
  
  use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int,C_char
  use LAMMPS
  !use MPI
  use atoms_class
  use parameters
  implicit none
  
  type,extends(atoms_group)        ::  lammps_obj
   type(C_ptr)                      :: lmp
   double precision, allocatable    :: bisp_der(:,:)
   double precision, allocatable    :: grads(:,:),stress(:,:)
   double precision                 :: edisp
   double precision                 :: e_shift
   double precision, allocatable    :: force_shift(:,:)
   integer,allocatable              :: topology(:)
   integer                          :: n_mol
   integer                          :: counter
   contains

   procedure     :: setup_lammps_lattice
   procedure     :: get_bis
   procedure     :: get_der_bis
   procedure     :: get_forces
   procedure     :: grimme_d3
   procedure     :: shift_mol

  end type lammps_obj

  contains
 
  subroutine number_bispec(twojmax,components)
  implicit none
  integer,intent(in)            :: twojmax
  integer,intent(out)           :: components
  double precision              :: order_coeff

  if (modulo(twojmax,2)==0) then
   order_coeff=(twojmax/2.0)+1
   components=order_coeff*(order_coeff+1)*(2*order_coeff+1)/6.0
  else
   order_coeff=(twojmax+1)/2.0
   components=order_coeff*(order_coeff+1)*(order_coeff+2)/3.0
  end if

  components=components + 1

  end subroutine number_bispec

  subroutine setup_lammps_lattice(this,nkinds)
  use atoms_class
  implicit none
  class(lammps_obj),intent(in)          :: this
  integer,intent(in)                    :: nkinds
  character(len=100)                    :: a1,a2,a3
  character(len=200)                    :: region_string
  character(len=250)                    :: create_atoms_string,mass_string,create_box_string
  integer                               :: j

  write(a1,*) 'a1',this%cell(1,1),this%cell(1,2),this%cell(1,3)
  write(a2,*) 'a2',this%cell(2,1),this%cell(2,2),this%cell(2,3)
  write(a3,*) 'a3',this%cell(3,1),this%cell(3,2),this%cell(3,3)

  write(region_string,*) '-50',this%cell(1,1)-50,'-50',this%cell(2,2)-50,'-50',this%cell(3,3)-50,this%cell(2,1),&
          this%cell(3,1),this%cell(3,2)
  
  call lammps_command(this%lmp,"units real") 
  call lammps_command(this%lmp,"dimension 3")
  call lammps_command(this%lmp,"boundary p p p")
  call lammps_command(this%lmp,"atom_style charge")
  call lammps_command(this%lmp,"neighbor 0.3 bin")
  call lammps_command(this%lmp,"neigh_modify one 10000")
  call lammps_command(this%lmp,'lattice custom 1.0'&
          //trim(a1)//trim(a2)//trim(a3)//' basis 0 0 0')
  call lammps_command(this%lmp,'region id_1 prism'//trim(region_string)//' units box')
  call lammps_command(this%lmp,"atom_modify map yes")
  call lammps_command(this%lmp,'box tilt large')
  write(create_box_string,*) 'create_box',nkinds,'id_1'
  call lammps_command(this%lmp,trim(create_box_string))
  
  do j=1,this%nats

   create_atoms_string=""

   write(create_atoms_string,*)'create_atoms',this%kind(j),'single',(this%x(j,1)),(this%x(j,2)),(this%x(j,3)),&
           'remap yes units box'
  call lammps_command(this%lmp,trim(create_atoms_string))

  end do

  do j=1,this%nkinds

    write(mass_string,*)"mass",j,this%mass(j)
    call lammps_command(this%lmp,mass_string)

  end do

  if (nkinds.ne.this%nkinds) then
   do j=this%nkinds+1,nkinds
    write(mass_string,*)"mass",j,'1'
    call lammps_command(this%lmp,mass_string)
   end do
  end if


  end subroutine setup_lammps_lattice
   
  subroutine get_der_bis(this,cutoff,twojmax,types)
  use atoms_class
  implicit none
  class(lammps_obj),intent(inout)                       :: this
  integer                                               :: i,j,k,pos,m
  integer (C_int), dimension(:),pointer                 :: id
  integer                                               :: nlocal
  integer,intent(in)                                    :: types
  integer,intent(out)                                   :: twojmax
  integer                                               :: components
  integer,allocatable                                   :: store_kind(:)
  real (C_double), dimension(:,:), pointer              :: x,bis_comp,bispec,der_bis
  double precision,intent(in)                           :: cutoff
  character(len=12000)                                    :: der_bis_string
                                                               
  
  call lammps_command(this%lmp,"pair_style zero 30")
  call lammps_command(this%lmp,"pair_coeff * *")
  write(der_bis_string,*)'compute der_bis all snad/atom',cutoff,'1',twojmax

  do i=1,types
   der_bis_string=trim(der_bis_string)//' 0.5'
  end do

  do i=1,types
   der_bis_string=trim(der_bis_string)//' 1'
  end do

  call lammps_command(this%lmp,trim(der_bis_string))
  
  call lammps_command(this%lmp,"run 0")

  call lammps_extract_compute(der_bis,this%lmp,'der_bis',LMP_STYLE_ATOM,LMP_TYPE_ARRAY)
  
  call lammps_extract_atom(id,this%lmp,"id")
  
  if (allocated(this%bisp_der)) deallocate(this%bisp_der)   
  if (.not.allocated(this%bisp_der)) allocate(this%bisp_der(size(der_bis,1),size(der_bis,2)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!columns refer to atoms, rows to order and component  

  do k=1,this%nats

   !pos=FINDLOC(id,k,1)
   do i=1,this%nats
    if (id(i)==k) then
   this%bisp_der(:,k)=der_bis(:,i)
   end if
   end do
  end do

  call lammps_command(this%lmp,"uncompute der_bis")

  !open(11,file='debug_snad',position="append")
  
    !     do k=1,size(this%bisp_der,1)
    !     write(11,*) this%bisp_der(k,:)
    !     end do

  !close(11)

  end subroutine get_der_bis

  subroutine get_bis(this,cutoff,twojmax)
  use atoms_class
  implicit none
  class(lammps_obj),intent(inout)                       :: this
  integer                                               :: i,j,k,pos,m
  integer (C_int), dimension(:),pointer                 :: id
  integer                                               :: nlocal
  integer,intent(out)                                   :: twojmax
  integer                                               :: components
  integer,allocatable                                   :: store_kind(:)
  logical                                               :: exist
  real (C_double), dimension(:,:), pointer              :: bispec=>NULL(),der_bis=>NULL()
  double precision,intent(in)                           :: cutoff
  character(len=150)                                    :: cutoff_string,der_bis_string
  real                                                  :: start,finish


  call lammps_command(this%lmp,"pair_style zero 30")
  call lammps_command(this%lmp,"pair_coeff * *")

  write(cutoff_string,*)'compute bispec all sna/atom',cutoff,'1',twojmax

  do i=1,this%nkinds
   
   cutoff_string=trim(cutoff_string)//' 0.5'
  
  end do

  do i=1,this%nkinds

   cutoff_string=trim(cutoff_string)//' 1'

  end do

  call lammps_command(this%lmp,trim(cutoff_string))
  call lammps_command(this%lmp,"run 0")

  call lammps_extract_atom(id,this%lmp,"id")
  call number_bispec(twojmax,components)
  
  !if (.not.associated(this%at_desc)) allocate(this%at_desc(this%nats))
  !allocate(this%at_desc(this%nats)) 

  call lammps_extract_compute(bispec,this%lmp,'bispec',LMP_STYLE_ATOM,LMP_TYPE_ARRAY)


  do k=1,this%nats

   !pos=FINDLOC(id,k,1)
  
   if (.not.allocated(this%at_desc(k)%desc)) allocate(this%at_desc(k)%desc(components))
   this%at_desc(k)%desc(1)=1.0

  do i=1,this%nats
    if (id(i)==k) then

    do m=1,components-1
     
     this%at_desc(k)%desc(m+1)=bispec(m,i)

    end do
   end if
  end do
   end do  

  call lammps_command(this%lmp,"uncompute bispec")

  end subroutine get_bis

  subroutine get_forces(this,screen_rad,F,bispec_flag,num_bisp_en,num_bisp_dip)
  implicit none
  class(lammps_obj),intent(in)                      :: this
  integer,intent(in)                                :: num_bisp_en,num_bisp_dip
  double precision, allocatable,intent(out)         :: F(:,:)
  double precision,intent(in)                       :: screen_rad
  double precision                                  :: dump
  integer                                           :: i,j,k,l,m
  double precision,allocatable                      :: beta_en(:),beta_dip(:)
  character(len=150)                                :: atom_string
  double precision,allocatable                      :: rel_dist(:,:),coord(:,:)
  character(len=6),intent(in)                       :: bispec_flag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Forces due to SNAP

allocate(beta_en((num_bisp_en-1)*this%nkinds))
allocate(beta_dip((num_bisp_dip-1)*this%nkinds))

open(222,file='snapcoeff_energy',action='read')

do i=1,this%nkinds

 read(222,*)
 do k=1,num_bisp_en-1

 read(222,*) beta_en((i-1)*(num_bisp_en-1)+k)

 end do

end do

close(222)


open(222,file='snapcoeff_dipoles',action='read')

do i=1,this%nkinds
 
 read(222,*)
 
  do k=1,num_bisp_dip-1

   read(222,*) beta_dip((i-1)*(num_bisp_dip-1)+k)

  end do

 end do

close(222)

beta_en=beta_en/Har_to_kc

allocate(rel_dist(this%nats,this%nats),coord(this%nats,3))

allocate(F(this%nats,3))

F=0.0

call this%dist_ij

coord=A_to_B*this%x

flag1: if (bispec_flag == "ENERGY") then
   
atom1:  do i=1,this%nats
   
nkinds:   do j=1,this%nkinds

cartesian_comp: do m=1,3

bis_comp:   do k=1,num_bisp_en-1

    F(i,m)=F(i,m)+this%bisp_der((j-1)*3*(num_bisp_en-1)+(m-1)*(num_bisp_en-1)+k,i)* & 
            beta_en((j-1)*(num_bisp_en-1)+k)

   end do bis_comp
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Forces due to electrostatic
 
!different_atoms:   if (i.ne.j) then
   
!    rel_dist(i,j)=A_to_B*this%dist(i,1,j,1)
   
!dumping:    if (rel_dist(i,j) <= screen_rad) then

 !    dump=0.5*(1-dcos(PI*rel_dist(i,j)/screen_rad))

  !  else

   !  dump=1.0
  
 !   end if dumping

  !  F(i,m)=F(i,m)+dump*(this%charge(i)*this%charge(j)*(coord(i,m)-coord(j,m)))/(rel_dist(i,j)**3)

   ! if (rel_dist(i,j) <= screen_rad) then
  
    ! F(i,m)=F(i,m)-0.5*this%charge(i)*this%charge(j)*(dsin(PI*rel_dist(i,j)/screen_rad))*PI* &
 !         (coord(i,m)-coord(j,m))/(screen_rad*(rel_dist(i,j)**2))
!
  !  end if

   !end if different_atoms
   
   end do cartesian_comp 

  end do nkinds

 end do atom1

end if flag1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! if (bispec_flag == "DIPOLE") then

 !do i=1,this%nats

  !do l=1,this%nats

   !do j=1,this%nats

    !if (l.ne.j) then

     !rel_dist(l,j)=A_to_B*this%dist(l,1,j,1)

      !if (rel_dist(l,j) <= screen_rad) then

       !dump=0.5*(1-dcos(PI*rel_dist(l,j)/screen_rad))

      !else

      ! dump=1.0

     ! end if

    !do m=1,3

     !do k=1,num_bisp_dip-1

      !F(i,m)=F(i,m)+dump*(this%charge(l)/rel_dist(l,j))*this%bisp_der((j-1)*3*(num_bisp_dip-1)+(m-1)*(num_bisp_dip-1)+k,i)* &
      !  beta_dip((this%kind(j)-1)*(num_bisp_dip-1)+k)
     
     !end do

    !end do

   ! end if

  ! end do

 ! end do

! end do

!end if
end subroutine get_forces

  subroutine grimme_d3(this,periodic_flag)
  use dftd3_api
  implicit none

  class(lammps_obj),intent(inout)        :: this
  double precision, allocatable          :: coords(:,:)
  double precision, allocatable          :: grads(:,:)
  double precision                       :: latVecs(3,3)
  logical                                :: periodic_flag

 ! integer, parameter :: species(nAtoms) = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
   !  & 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, &
   !  & 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]

  ! Lattice vectors in Angstrom as found in dna.xyz/dna.poscar
  ! They must be converted to Bohr before passed to dftd3
  !real(wp), parameter :: latVecs(3, 3) = reshape([&
   !    &  8.0000000000E+00,   0.0000000000E+00,   0.0000000000E+00, &
   !    &  0.0000000000E+00,   8.0000000000E+00,   0.0000000000E+00, &
   !    &  0.0000000000E+00,   0.0000000000E+00,   1.5000000000E+01  &
   !    & ] * AA__Bohr, [3, 3])

  !integer, parameter :: nSpecies = 4
  !character(2),  :: speciesNames(nSpecies) = [ 'N ', 'C ', 'O ', 'H ']

  integer,allocatable           :: atnum(:)
  type(dftd3_input)             :: input
  type(dftd3_calc)              :: dftd3
  double precision              :: edisp

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialize input
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! You can set input variables if you like, or just leave them on their
  ! defaults, which are the same as the dftd3 program uses.

  !! Threebody interactions (default: .false.)
  !input%threebody = .true.
  !
  !! Numerical gradients (default: .false.)
  !input%numgrad = .false.
  !
  !! Cutoffs (below you find the defaults)
  !input%cutoff = sqrt(9000.0_wp)
  !input%cutoff_cn = sqrt(1600.0_wp)

  ! Initialize dftd3
  call dftd3_init(dftd3, input)

  ! Choose functional. Alternatively you could set the parameters manually
  ! by the dftd3_set_params() function.
  call dftd3_set_functional(dftd3, func='pbe', version=4, tz=.false.)

  allocate(atnum(this%nats))
  ! Convert species name to atomic number for each atom
  atnum(:) = get_atomic_number(this%label(this%kind))

  allocate(coords(3,this%nats))
  coords=transpose(this%x)*A_to_B
  !allocate(this%grads(3,this%nats))

  if (.not.periodic_flag) then
  ! Calculate dispersion and gradients for non-periodic case
  call dftd3_dispersion(dftd3, coords, atnum, this%edisp, this%grads)
  !write(*, "(A)") "*** Dispersion for non-periodic case"
  !open(111, file='VdW_ener_MolForge_20220225.txt',action='write',position='append')
  !write(111,*)  this%edisp*Har_to_kc
  !close(111)
  !write(*, "(A)") "Gradients [au]:"
  !write(*, "(3ES20.12)") this%grads
   ! write(*, *) this%grads(1,:)
  !!write(*, *)
  else
  latVecs(1,:)=this%cell(:,1)*A_to_B
  latVecs(2,:)=this%cell(:,2)*A_to_B
  latVecs(3,:)=this%cell(:,3)*A_to_B
  ! Calculate dispersion and gradients for periodic case
  if (.not.allocated(this%stress)) allocate(this%stress(3,3))
  call dftd3_pbc_dispersion(dftd3, coords, atnum, latVecs, this%edisp, this%grads, this%stress)
  !write(*, "(A)") "*** Dispersion for periodic case"
  !write(*, "(A,ES20.12)") "Energy [au]:", edisp
  !write(*, "(A)") "Gradients [au]:"
  !write(*, "(3ES20.12)") grads
  !write(*, "(A)") "Stress [au]:"
  !write(*, "(3ES20.12)") stress
  end if

  end subroutine grimme_d3

subroutine shift_mol(this,k,C_M0)
implicit none
class(lammps_obj),intent(inout)                   :: this
integer                                           :: i,j
double precision,allocatable                      :: C_M(:,:),mass_mol(:)
double precision,intent(in)                       :: k,C_M0

this%n_mol=2

allocate(C_M(this%n_mol,3),mass_mol(this%n_mol))
if (.not. allocated(this%force_shift)) allocate(this%force_shift(this%nats,3))

C_M=0.0d0
mass_mol=0.0d0
this%force_shift=0.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!calculation of masses and centers of mass of molecules
do i=1,this%nats

   mass_mol(this%topology(i))=mass_mol(this%topology(i))+this%mass(this%kind(i))*amu_to_emass
   C_M(this%topology(i),:)=C_M(this%topology(i),:)+this%mass(this%kind(i))*amu_to_emass*this%x(i,:)*A_to_B

end do

do j=1,this%n_mol

C_M(j,:)=C_M(j,:)/mass_mol(j)

end do

open(111, file='dist_CM.txt',action='write',position='append')
write(111,*) 'Distance between centers of mass (angstrom):',norm2(C_M(1,:)-C_M(2,:))/A_to_B
close(111)

this%e_shift=0.5d0*k*(norm2(C_M(1,:)-C_M(2,:))-C_M0)**2

do i=1,this%nats
  
 if (this%topology(i)==1) then

 this%force_shift(i,:)=k*(C_M0-norm2(C_M(1,:)-C_M(2,:)))*(C_M(1,:)-C_M(2,:)) &
         *amu_to_emass*this%mass(this%kind(i))/(norm2(C_M(1,:)-C_M(2,:))*mass_mol(1))
 
else
 
         this%force_shift(i,:)=k*(C_M0-norm2(C_M(1,:)-C_M(2,:)))*(C_M(2,:)-C_M(1,:)) &
         *amu_to_emass*this%mass(this%kind(i))/(norm2(C_M(1,:)-C_M(2,:))*mass_mol(2))
 
 end if


end do



end subroutine shift_mol

end module lammps_class 

