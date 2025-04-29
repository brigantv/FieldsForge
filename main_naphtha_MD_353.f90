        program main
        use lammps_class
        use geom_class
        use SNAP_class
        use max_class
        use dist_class
        use gradmin_class
        use force_field_min_class
        use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
        use LAMMPS
        !use MPI

        implicit none

        type(lammps_obj),allocatable              :: set(:)
        type(lammps_obj),allocatable              :: set_tmp(:)
        integer                                   :: nconfig,tot_kinds,num_bisp_en,num_bisp_dip,twojmax_dip,twojmax_en
        integer                                   :: i,j
        character(len=120)                        :: geometry_file,energy_file,dipoles_file,shift_file,atom_string, &
                                                        shift_geo_file,forces_file
        logical,dimension(:),allocatable          :: coeff_mask_en,coeff_mask_dip
        double precision                          :: lambda_en,cutoff_en
        double precision                          :: lambda_dip,cutoff_dip
        logical                                   :: dipole_flag,energy_flag,md_flag,VdW_flag,coul_flag,minim_flag, &
                train_ff,rampa_flag,periodic_flag
        logical                                   :: flag_forces,flag_energy
        double precision,dimension(:),allocatable :: energies,coul_energy,snap_energy
        double precision,dimension(:),allocatable :: DFT_forces
        double precision                          :: R_screen
        double precision,dimension(:),allocatable :: tot_charge
        character(len=5)                          :: set_type_en,set_type_dip
        double precision,allocatable              :: forces(:,:),forces_1(:,:),forces_model(:,:)
        double precision,allocatable              :: force(:,:)
        double precision                          :: E_plus,E_minus,ave_atom
        integer                                   :: tot_atom
        integer                                   :: atom,direction,k
        double precision,dimension(:),allocatable :: coul_energy_plus,coul_energy_minus
        double precision                          :: temperature,timestep
        integer                                   :: tot_steps_md
        type(force_field_min),pointer             :: min_force_field
        double precision                          :: val,step
        double precision,allocatable              :: vec(:),grad(:)
        type(adam)                                :: minimizer
        double precision,allocatable              :: VdW_en(:),test_1(:,:),test_2(:,:)
        double precision,allocatable              :: grads(:,:)
        double precision                          :: edisp
        double precision                          :: test(3,3)
        integer                                       :: INFO
        integer                                     :: len_shift_file
        integer                                      :: idist
        integer,dimension(4)                         :: iseed
        integer                                      :: N
        double precision, allocatable                ::X_M(:), Y(:),X(:),mean(:),sigma(:)
        type(dist1D)                                 :: list
        double precision                            :: y_tmp,y_tmp_2
        character (len=255)                             ::      cwd
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!What do you want to do?
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Achtung: you cannot turn only the VdW flag but it must be always turned together with coul flag

allocate(min_force_field)

        train_ff=.false.
        VdW_flag=.true. !refers only to FFs
        periodic_flag=.true.
        coul_flag=.false. !refers only to FFs
        min_force_field%md_flag=.false.
        min_force_field%temperature_in=10.0d0
        min_force_field%temperature_final=50.0d0
        min_force_field%timestep=1.0d0*41.49d0
        min_force_field%nsteps=5000
        min_force_field%iseed=[1216,364,1903,4059]
        min_force_field%periodic_flag=.true.
        min_force_field%minim_flag=.false.
        min_force_field%active_learn=.true.
        min_force_field%sigma_AL=dsqrt(2.0d0)*10.0d0
        min_force_field%thresh_AL=0.5d0
        min_force_field%rampa_flag=.false.
        min_force_field%factor_thresh=1.5d0
        min_force_field%shift_flag=.false.
        min_force_field%k_AL=0.006d0
        min_force_field%C_M0=3.5d0*A_to_B
        min_force_field%weight=dsqrt(3.0d0*1620.0d0)

        !this two flags refer both to minimization and molecular dynamics
        min_force_field%VdW_flag=.true. 
        min_force_field%coul_flag=.false.
        call getcwd(cwd)
!!!!!!!!/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!INPUT PARAMETERS
        nconfig=22
        min_force_field%nconfig_AL=22

        geometry_file="/home/users/brigantv/active_learning/naphtha_MD_353/geo_tr_AL_++"
        energy_file="/home/users/brigantv/active_learning/naphtha_MD_353/ener_tr_AL_++"
        dipoles_file="/home/users/brigantv/active_learning/naphtha_MD_353/dipoles_tr_AL_++"
        forces_file="/home/users/brigantv/active_learning/naphtha_MD_353/forces_tr_AL_++"
        
        !shift_file="/home/valeriobriganti/Desktop/MolForge_develop/FitSnap/SCO/centers_mass_AL_++"
        
        set_type_en="TRAIN"
        set_type_dip="TRAIN"
        flag_energy=.true.
        flag_forces=.true.

        len_shift_file=len_trim(shift_file)
        call import_geom(set,nconfig,trim(geometry_file),len_trim(geometry_file))
        !call import_geom(set_tmp,nconfig,trim(geometry_file),len_trim(geometry_file))
        call import_geom(min_force_field%set_AL,min_force_field%nconfig_AL,trim(geometry_file),len_trim(geometry_file))
        
        lambda_en=0.1
        twojmax_en=8
        cutoff_en=3.5
        R_screen=4.0*A_to_B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CHARGES
        
        lambda_dip=0.0
        cutoff_dip=4.0
        twojmax_dip=8

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
        call get_tot_kinds(set,tot_kinds)
        allocate(coeff_mask_en(tot_kinds))
        allocate(coeff_mask_dip(tot_kinds))
        allocate(tot_charge(size(set)))

!the coeffiecients that you set to be true are the ones that you put equal to zero
        
        coeff_mask_en=.true.
        coeff_mask_en(1)=.false.

        coeff_mask_dip=.true.
        coeff_mask_dip(1)=.false.

        tot_charge=0.0
        min_force_field%object_lammps%n_mol=1
        call number_bispec(twojmax_en,num_bisp_en)
        call number_bispec(twojmax_dip,num_bisp_dip)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (train_ff) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!READING ENERGIES (this should be in fit energy subroutine)
                
        allocate(energies(size(set)),VdW_en(size(set)))
        VdW_en=0.0

        open(100,file=trim(energy_file))
        
        do j=1,nconfig
        
         read(100,*) energies(j)

        end do

        close(100)

if (flag_forces) then
call get_ave_atoms(set,ave_atom,tot_atom)
 allocate(DFT_forces(3*tot_atom))

 open(100,file=trim(forces_file))

  do j=1,3*tot_atom

         read(100,*) DFT_forces(j)

  end do

        close(100)
         
         if (VdW_flag) then

          do i=1,nconfig

           allocate(set(i)%grads(3,set(i)%nats))
           call set(i)%grimme_d3(periodic_flag)

            do j=1,set(i)%nats
             do k=1,3

              DFT_forces(((i-1)*set(i)%nats*3)&
              +(j-1)*3+k)=DFT_forces(((i-1)*set(i)%nats*3) +(j-1)*3+k)-set(i)%grads(k,j)

             end do
            end do

           deallocate(set(i)%grads)
          end do

         end if


        DFT_forces=-DFT_forces*F_conv

        end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        
       do i=1,nconfig
         
         call lammps_open_no_mpi("lmp -screen none -log log.simple",set(i)%lmp)
         
         call set(i)%setup_lammps_lattice(set(i)%nkinds)

         call set(i)%get_bis(cutoff_dip,twojmax_dip)
         
         if (.not.allocated(set(i)%grads)) allocate(set(i)%grads(3,set(i)%nats))

         call set(i)%grimme_d3(periodic_flag)
          
         VdW_en(i)=set(i)%edisp*Har_to_kc
         
         deallocate(set(i)%grads)
         call lammps_close(set(i)%lmp)

       end do
         call fit_dipoles(set,num_bisp_dip,trim(dipoles_file),lambda_dip,coeff_mask_dip,set_type_dip,tot_charge,&
                 len_trim(dipoles_file))
         call SNAP_coulomb_energy(set,R_screen,coul_energy)

          if ((coul_flag).and.(VdW_flag)) then

           energies=energies-coul_energy-VdW_en

          else if (coul_flag) then

           energies=energies-coul_energy

          else if (VdW_flag) then

           energies=energies-VdW_en


          end if

               open(111, file='fit_energy_rms.dat', action='write')

if (flag_energy) then

         do i=1,nconfig

          write(111,*) energies(i)

         end do

end if
         
if (flag_forces) then

         do i=1,3*tot_atom

          write(111,*) DFT_forces(i)

         end do
 end if
close(111)
        
        do i=1,nconfig

         call lammps_open_no_mpi("lmp -screen none -log log.bis_en",set(i)%lmp)

         call set(i)%setup_lammps_lattice(set(i)%nkinds)

         call set(i)%get_bis(cutoff_en,twojmax_en)
 
        call set(i)%get_der_bis(cutoff_en,twojmax_en,set(i)%nkinds)
         
         call lammps_close(set(i)%lmp)

        end do
         
        call fit_energy_forces(set,num_bisp_en,energies,DFT_forces,lambda_en,coeff_mask_en,set_type_en,snap_energy,&
              flag_energy,flag_forces,min_force_field%weight)
        
       if ((coul_flag).and.(VdW_flag)) then

         energies=energies+coul_energy+VdW_en

        else if (coul_flag) then

         energies=energies+coul_energy

         else if (VdW_flag) then

           energies=energies+VdW_en

        end if

 open(111,file="energie_fitting.dat",action="write")        
 do i=1,min_force_field%nconfig_AL

                if ((coul_flag).and.(VdW_flag)) then

        write(111,*)  snap_energy+coul_energy+VdW_en

        else if (coul_flag) then

write(111,*)  snap_energy+coul_energy
         else if (VdW_flag) then
write(111,*)  snap_energy+VdW_en

        end if

 end do

end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       do i=1,nconfig

        call lammps_open_no_mpi("lmp -screen none -log log.simple",set(i)%lmp)

        call set(i)%setup_lammps_lattice(set(i)%nkinds)

        call set(i)%get_bis(cutoff_dip,twojmax_dip)

        call lammps_close(set(i)%lmp)

       end do
 if (.not. allocated(min_force_field%object_lammps%topology)) &
         allocate(min_force_field%object_lammps%topology(min_force_field%object_lammps%nats))
        min_force_field%R_screen=R_screen
        min_force_field%set=set
        min_force_field%tot_charge=tot_charge
        min_force_field%cutoff_dip=cutoff_dip
        min_force_field%twojmax_dip=twojmax_dip
        min_force_field%cutoff_en=cutoff_en
        min_force_field%twojmax_en=twojmax_en
        min_force_field%num_bisp_en=num_bisp_en
        min_force_field%lambda_en=lambda_en
        min_force_field%num_bisp_dip=num_bisp_dip
        min_force_field%coeff_mask_dip=coeff_mask_dip
        min_force_field%coeff_mask_en=coeff_mask_en
        min_force_field%object_lammps=set(1)
        min_force_field%object_lammps%topology=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
        min_force_field%dipoles_file=dipoles_file
        min_force_field%geometry_file=geometry_file
        min_force_field%energy_file=energy_file
        min_force_field%flag_energy=flag_energy
        min_force_field%flag_forces=flag_forces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Minimization routine

if (min_force_field%minim_flag) then

        min_force_field%lr=0.0001d0
        call min_force_field%minimize_adam(max_iter=500000,start_iter=0)

end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!LAMMPS molecular dynamics

if (min_force_field%md_flag) then

  call min_force_field%init_vel(min_force_field%temperature_in)
  call min_force_field%propagate_md(min_force_field%nsteps,min_force_field%timestep)

   
end if

end program
