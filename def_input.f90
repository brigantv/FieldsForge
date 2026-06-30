module input_class

  use parameters
  implicit none

  integer, parameter :: MAX_KINDS = 20
  integer, parameter :: MAX_ATOMS = 1000

  type snap_input

    ! Workflow
    logical            :: train_ff
    logical            :: VdW_flag
    logical            :: coul_flag
    logical            :: periodic_flag

    ! Training data
    integer            :: nconfig
    character(len=200) :: geometry_file
    character(len=200) :: energy_file
    character(len=200) :: forces_file
    character(len=200) :: dipoles_file

    ! Fit targets
    logical            :: flag_energy
    logical            :: flag_forces
    logical            :: flag_stress
    character(len=5)   :: set_type_en
    character(len=5)   :: set_type_dip

    ! SNAP energy descriptor
    integer            :: twojmax_en
    double precision   :: cutoff_en
    double precision   :: lambda_en
    double precision   :: weight
    logical            :: ezero_en(MAX_KINDS)   ! .true. = constrain E0 to zero

    ! SNAP dipole descriptor
    integer            :: twojmax_dip
    double precision   :: cutoff_dip
    double precision   :: lambda_dip
    logical            :: ezero_dip(MAX_KINDS)

    ! Coulomb screening radius (user gives Angstroms; stored in Bohr)
    double precision   :: R_screen

    ! MD
    logical            :: md_flag
    integer            :: nsteps
    double precision   :: timestep           ! user: fs; stored: atomic time units (*41.49)
    double precision   :: temperature_in
    double precision   :: temperature_final
    integer            :: iseed(4)

    ! Geometry minimization
    logical            :: minim_flag
    double precision   :: lr
    integer            :: max_iter_adam

    ! Active learning
    logical            :: active_learn
    integer            :: nconfig_AL          ! defaults to nconfig if left 0
    double precision   :: sigma_AL
    double precision   :: thresh_AL
    double precision   :: factor_thresh

    ! Molecule topology (only first nats elements are used)
    integer            :: n_mol
    integer            :: topology(MAX_ATOMS)
    integer            :: fixed_atoms(MAX_ATOMS)  ! 1 = fixed, 0 = free

    ! Misc
    logical            :: rampa_flag
    logical            :: shift_flag
    double precision   :: k_AL
    double precision   :: C_M0               ! user: Angstroms; stored: Bohr

  end type snap_input

contains

  subroutine read_snap_input(filename, inp)
    character(len=*),  intent(in)  :: filename
    type(snap_input),  intent(out) :: inp

    ! Namelist mirrors (can't use derived-type fields directly in namelist)
    logical            :: train_ff, VdW_flag, coul_flag, periodic_flag
    integer            :: nconfig
    character(len=200) :: geometry_file, energy_file, forces_file, dipoles_file
    logical            :: flag_energy, flag_forces, flag_stress
    character(len=5)   :: set_type_en, set_type_dip
    integer            :: twojmax_en, twojmax_dip
    double precision   :: cutoff_en, cutoff_dip, lambda_en, lambda_dip, weight
    logical            :: ezero_en(MAX_KINDS), ezero_dip(MAX_KINDS)
    double precision   :: R_screen
    logical            :: md_flag, minim_flag, active_learn
    integer            :: nsteps, max_iter_adam, nconfig_AL
    double precision   :: timestep, temperature_in, temperature_final
    integer            :: iseed(4)
    double precision   :: lr, sigma_AL, thresh_AL, factor_thresh
    integer            :: n_mol, topology(MAX_ATOMS), fixed_atoms(MAX_ATOMS)
    logical            :: rampa_flag, shift_flag
    double precision   :: k_AL, C_M0

    namelist /SNAP/ &
      train_ff, VdW_flag, coul_flag, periodic_flag,                     &
      nconfig, geometry_file, energy_file, forces_file, dipoles_file,   &
      flag_energy, flag_forces, flag_stress, set_type_en, set_type_dip, &
      twojmax_en, cutoff_en, lambda_en, weight, ezero_en,               &
      twojmax_dip, cutoff_dip, lambda_dip, ezero_dip,                   &
      R_screen, md_flag, nsteps, timestep, temperature_in,              &
      temperature_final, iseed, minim_flag, lr, max_iter_adam,          &
      active_learn, nconfig_AL, sigma_AL, thresh_AL, factor_thresh,     &
      n_mol, topology, fixed_atoms, rampa_flag, shift_flag, k_AL, C_M0

    ! Defaults (match the original hardcoded values in main.f90)
    train_ff          = .true.
    VdW_flag          = .false.
    coul_flag         = .false.
    periodic_flag     = .true.
    nconfig           = 0
    geometry_file     = ""
    energy_file       = ""
    forces_file       = ""
    dipoles_file      = ""
    flag_energy       = .true.
    flag_forces       = .true.
    flag_stress       = .false.
    set_type_en       = "TRAIN"
    set_type_dip      = "TRAIN"
    twojmax_en        = 8
    cutoff_en         = 4.0d0
    lambda_en         = 0.0d0
    weight            = 1.0d0
    ezero_en          = .false.
    twojmax_dip       = 8
    cutoff_dip        = 4.0d0
    lambda_dip        = 0.0d0
    ezero_dip         = .false.
    R_screen          = 4.0d0          ! Angstroms
    md_flag           = .false.
    nsteps            = 1000
    timestep          = 1.0d0          ! femtoseconds
    temperature_in    = 300.0d0
    temperature_final = 300.0d0
    iseed             = [1216, 364, 1903, 4059]
    minim_flag        = .false.
    lr                = 0.001d0
    max_iter_adam     = 500000
    active_learn      = .false.
    nconfig_AL        = 0
    sigma_AL          = dsqrt(2.0d0) * 10.0d0
    thresh_AL         = 0.5d0
    factor_thresh     = 1.5d0
    n_mol             = 1
    topology          = 1              ! all atoms in molecule 1 by default
    fixed_atoms       = 0              ! all atoms free by default
    rampa_flag        = .false.
    shift_flag        = .false.
    k_AL              = 0.006d0
    C_M0              = 3.5d0          ! Angstroms

    open(99, file=trim(filename), status='old', action='read')
    read(99, nml=SNAP)
    close(99)

    ! Unit conversions
    R_screen = R_screen * A_to_B      ! Angstroms -> Bohr
    C_M0     = C_M0     * A_to_B      ! Angstroms -> Bohr
    timestep = timestep * 41.49d0     ! fs -> atomic time units

    if (nconfig_AL == 0) nconfig_AL = nconfig

    ! Copy to output type
    inp%train_ff          = train_ff
    inp%VdW_flag          = VdW_flag
    inp%coul_flag         = coul_flag
    inp%periodic_flag     = periodic_flag
    inp%nconfig           = nconfig
    inp%geometry_file     = geometry_file
    inp%energy_file       = energy_file
    inp%forces_file       = forces_file
    inp%dipoles_file      = dipoles_file
    inp%flag_energy       = flag_energy
    inp%flag_forces       = flag_forces
    inp%flag_stress       = flag_stress
    inp%set_type_en       = set_type_en
    inp%set_type_dip      = set_type_dip
    inp%twojmax_en        = twojmax_en
    inp%cutoff_en         = cutoff_en
    inp%lambda_en         = lambda_en
    inp%weight            = weight
    inp%ezero_en          = ezero_en
    inp%twojmax_dip       = twojmax_dip
    inp%cutoff_dip        = cutoff_dip
    inp%lambda_dip        = lambda_dip
    inp%ezero_dip         = ezero_dip
    inp%R_screen          = R_screen
    inp%md_flag           = md_flag
    inp%nsteps            = nsteps
    inp%timestep          = timestep
    inp%temperature_in    = temperature_in
    inp%temperature_final = temperature_final
    inp%iseed             = iseed
    inp%minim_flag        = minim_flag
    inp%lr                = lr
    inp%max_iter_adam     = max_iter_adam
    inp%active_learn      = active_learn
    inp%nconfig_AL        = nconfig_AL
    inp%sigma_AL          = sigma_AL
    inp%thresh_AL         = thresh_AL
    inp%factor_thresh     = factor_thresh
    inp%n_mol             = n_mol
    inp%topology          = topology
    inp%fixed_atoms       = fixed_atoms
    inp%rampa_flag        = rampa_flag
    inp%shift_flag        = shift_flag
    inp%k_AL              = k_AL
    inp%C_M0              = C_M0

  end subroutine read_snap_input

end module input_class
