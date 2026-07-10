        program main
        use parameters
        use lammps_class
        use geom_class
        use SNAP_class
        use max_class
        use dist_class
        use gradmin_class
        use force_field_min_class
        use input_class
        use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
        use LAMMPS
        !use MPI

        implicit none

        type(lammps_obj),allocatable              :: set(:)
        type(lammps_obj),allocatable              :: set_tmp(:)
        integer                                   :: nconfig,tot_kinds,num_bisp_en,num_bisp_dip,twojmax_dip,twojmax_en
        integer                                   :: i,j
        character(len=200)                        :: geometry_file,energy_file,dipoles_file,forces_file,record_file
        character(len=120)                        :: shift_file,atom_string,shift_geo_file
        logical,dimension(:),allocatable          :: coeff_mask_en,coeff_mask_dip
        double precision                          :: lambda_en,cutoff_en
        double precision                          :: lambda_dip,cutoff_dip
        logical                                   :: dipole_flag,energy_flag,md_flag,VdW_flag,coul_flag,minim_flag, &
                train_ff,rampa_flag,periodic_flag
        logical                                   :: flag_forces,flag_energy,flag_stress
        double precision,dimension(:),allocatable :: energies,coul_energy,snap_energy
        double precision,dimension(:),allocatable   :: DFT_forces
        double precision,dimension(:,:),allocatable :: DFT_stress
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
        double precision,allocatable              :: VdW_en(:),test_1(:,:),test_2(:,:)
        double precision,allocatable              :: grads(:,:)
        double precision                          :: edisp
        double precision                          :: test(3,3)
        integer                                   :: INFO
        integer                                   :: len_shift_file
        integer                                   :: idist
        integer,dimension(4)                      :: iseed
        integer                                   :: N
        double precision, allocatable             :: X_M(:), Y(:),X(:),mean(:),sigma(:)
        type(dist1D)                              :: list
        double precision                          :: y_tmp,y_tmp_2
        character (len=255)                       :: cwd
        type(snap_input)                          :: inp
        character(len=200)                        :: input_filename

        allocate(min_force_field)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Read input file
        !input_filename = "snap.input"
        call get_command_argument(1, input_filename)
        call read_snap_input(trim(input_filename), inp)

        train_ff                          = inp%train_ff
        VdW_flag                          = inp%VdW_flag
        periodic_flag                     = inp%periodic_flag
        coul_flag                         = inp%coul_flag
        min_force_field%md_flag           = inp%md_flag
        min_force_field%temperature_in    = inp%temperature_in
        min_force_field%temperature_final = inp%temperature_final
        min_force_field%timestep          = inp%timestep
        min_force_field%nsteps            = inp%nsteps
        min_force_field%iseed             = inp%iseed
        min_force_field%periodic_flag     = inp%periodic_flag
        min_force_field%minim_flag        = inp%minim_flag
        min_force_field%active_learn      = inp%active_learn
        min_force_field%sigma_AL          = inp%sigma_AL
        min_force_field%thresh_AL         = inp%thresh_AL
        min_force_field%rampa_flag        = inp%rampa_flag
        min_force_field%factor_thresh     = inp%factor_thresh
        min_force_field%shift_flag        = inp%shift_flag
        min_force_field%k_AL              = inp%k_AL
        min_force_field%C_M0              = inp%C_M0
        min_force_field%weight            = inp%weight
        min_force_field%VdW_flag          = inp%VdW_flag
        min_force_field%coul_flag         = inp%coul_flag
        min_force_field%lr                = inp%lr

        call getcwd(cwd)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!INPUT PARAMETERS
        nconfig                    = inp%nconfig
        min_force_field%nconfig_AL = inp%nconfig_AL

        geometry_file = trim(inp%geometry_file)
        energy_file   = trim(inp%energy_file)
        dipoles_file  = trim(inp%dipoles_file)
        forces_file   = trim(inp%forces_file)
        record_file   = trim(inp%record_file)

        set_type_en  = inp%set_type_en
        set_type_dip = inp%set_type_dip
        flag_energy  = inp%flag_energy
        flag_forces  = inp%flag_forces
        flag_stress  = inp%flag_stress

        len_shift_file=len_trim(shift_file)
        call import_geom(set,nconfig,trim(geometry_file),len_trim(geometry_file))
        call import_geom(set_tmp,nconfig,trim(geometry_file),len_trim(geometry_file))
        call import_geom(min_force_field%set_AL,min_force_field%nconfig_AL,trim(geometry_file),len_trim(geometry_file))

        lambda_en  = inp%lambda_en
        twojmax_en = inp%twojmax_en
        cutoff_en  = inp%cutoff_en
        R_screen   = inp%R_screen    ! already in Bohr (converted in read_snap_input)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CHARGES

        lambda_dip  = inp%lambda_dip
        cutoff_dip  = inp%cutoff_dip
        twojmax_dip = inp%twojmax_dip

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call get_tot_kinds(set,tot_kinds)
        allocate(coeff_mask_en(tot_kinds))
        allocate(coeff_mask_dip(tot_kinds))
        allocate(tot_charge(size(set)))

        coeff_mask_en  = inp%ezero_en(1:tot_kinds)
        coeff_mask_dip = inp%ezero_dip(1:tot_kinds)

        tot_charge=0.0
        min_force_field%object_lammps%n_mol = inp%n_mol
        call number_bispec(twojmax_en,num_bisp_en)
        call number_bispec(twojmax_dip,num_bisp_dip)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!2!!!!!!!!!!!!!!!!!!!!!!

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

if (flag_stress) then
allocate(DFT_stress(3*nconfig,3))
!open(100,file=trim(stress_file))
!do j=0,nconfig-1
!        read(100,*) DFT_stress(3*j+1,:)
!        read(100,*) DFT_stress(3*j+2,:)
!        read(100,*) DFT_stress(3*j+3,:)
!end do
!close(100)
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i=1,nconfig

         call lammps_open_no_mpi("lmp -screen none -log log.simple",set(i)%lmp)
         call set(i)%setup_lammps_lattice(set(i)%nkinds)
        ! call set(i)%calc_vol                   !calculate the volume of the cell for every structure in the set
        ! call set(i)%calc_Js                    !calculate the jacobians
        ! call set(i)%calc_Gs                    !calculate the metric tensors for the direct and reciprocal lattice
        ! set(i)%Gcon=(4.0*PI)**2*set(i)%Gcon    !rescaling the metric tensor of the reciprocal lattice to
         call set(i)%get_bis(cutoff_dip,twojmax_dip)

        ! if (.not.allocated(set(i)%grads)) allocate(set(i)%grads(3,set(i)%nats))
        ! call set(i)%grimme_d3(periodic_flag)
        ! VdW_en(i)=set(i)%edisp*Har_to_kc
        ! deallocate(set(i)%grads)
         call lammps_close(set(i)%lmp)
       end do
        call fit_dipoles(set,num_bisp_dip,trim(dipoles_file),lambda_dip,coeff_mask_dip,set_type_dip,tot_charge,&
                 len_trim(dipoles_file))
        call SNAP_coulomb_energy(set,R_screen,coul_energy)

        ! testing the routine for the Ewald summation
        ! insert here the call to the subroutine to test
        !set(1)%comp_max=3.0
        !set(1)%ew_alpha=2.0
        !call set(1)%get_ewald

        

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

         call set(i)%get_stress(cutoff_en,twojmax_en)

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
        allocate(min_force_field%fixed_atoms(min_force_field%object_lammps%nats))
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
        min_force_field%object_lammps%topology=inp%topology(1:min_force_field%object_lammps%nats)
        min_force_field%fixed_atoms=inp%fixed_atoms(1:min_force_field%object_lammps%nats)
        min_force_field%dipoles_file=dipoles_file
        min_force_field%geometry_file=geometry_file
        min_force_field%energy_file=energy_file
        min_force_field%record_file=record_file
        min_force_field%flag_energy=flag_energy
        min_force_field%flag_forces=flag_forces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Minimization routine

if (min_force_field%minim_flag) then

        call min_force_field%minimize_adam(max_iter=inp%max_iter_adam,start_iter=0)

end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!LAMMPS molecular dynamics

if (min_force_field%md_flag) then

  call min_force_field%init_vel(min_force_field%temperature_in)
  call min_force_field%propagate_md(min_force_field%nsteps,min_force_field%timestep)

end if



end program
