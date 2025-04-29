  module geom_class
  
  use atoms_class
  use lammps_class

  implicit none
  
  contains

  subroutine get_centroid(set)

  type(lammps_obj), allocatable      :: set(:)
  integer                            :: i,j
  double precision,dimension(3)                   :: x_centr
  
  open(222,file='centroid_for_dipole',action='write')
  
  do i=1,size(set)
   
   x_centr=0.0
   
   do j=1,set(i)%nats

   x_centr=x_centr+set(i)%x(j,:)
   
   end do

   x_centr=x_centr/set(i)%nats

  write(222,*)'OriXYZ', x_centr(1),',',x_centr(2),',',x_centr(3)

  end do

  close(222) 

  end subroutine get_centroid

  subroutine import_geom(set,nconfig,file_input,len_file_inp)
    
  type(lammps_obj), allocatable    ::  set(:)
  integer                          :: nconfig, i,j,nats,ntypes
  character(len=100),allocatable   :: tmp(:,:)
  character(len=100),dimension(10) :: tmp_cell_nkinds
  character(len=100)                :: filename
  character(len=*)               :: file_input
  integer,intent(in)               :: len_file_inp
  character(len=150)               ::file_inp
  !declare dimension allocatable above together with the fixed one, why everything should be allocatable?
  
  file_inp=file_input(1:len_file_inp)

  allocate(set(nconfig))
  !allocate(tmp_cell_nkinds(10))
  
  open(1,file=trim(file_inp))

  do i=1,nconfig

   read(1,*)nats
   read(1,*)tmp_cell_nkinds(:)

   allocate(tmp(nats,6))

   do j=1,nats

    read(1,*) tmp(j,:)

   end do

   filename="file_temporaneo"

   open(unit=2,status="replace",file=filename)

   write(2,*)nats
   write(2,*)trim(tmp_cell_nkinds(1)),' ',trim(tmp_cell_nkinds(2)),' ',trim(tmp_cell_nkinds(3)),' '&
           ,trim(tmp_cell_nkinds(4)),' ',trim(tmp_cell_nkinds(5)),' ',trim(tmp_cell_nkinds(6)),' ',trim(tmp_cell_nkinds(7))&
           ,' ',trim(tmp_cell_nkinds(8)),' ',trim(tmp_cell_nkinds(9)),' ',trim(tmp_cell_nkinds(10))

   do j=1,nats

    write(2,*)trim(tmp(j,1)),' ',trim(tmp(j,2)),' ',trim(tmp(j,3)),' ',trim(tmp(j,4)),' ',trim(tmp(j,5)),&
            ' ',trim(tmp(j,6))
   
   end do
   
   deallocate(tmp)
   !deallocate(tmp_cell_nkinds)

   close(2)

   call set(i)%read_extended_xyz(3,filename)

  end do

  close(1)

  end subroutine import_geom

  subroutine shift_geom(set,filename,len_shift_file)
  type(lammps_obj), allocatable,intent(inout)      :: set(:)
  integer                                          :: i,j
  character(len=150),intent(in)                     :: filename
  integer,intent(in)                               :: len_shift_file
  double precision,allocatable                     :: shifter(:,:)
 character(len=150)                                :: file_shift

  allocate(shifter(size(set),3))

  file_shift=filename(1:len_shift_file)
   open(unit=222,file=file_shift,action='read')

  do j=1,size(set)

   read(222,*) shifter(j,1),shifter(j,2),shifter(j,3)
  
  end do

  close(222)

  do i=1,size(set)
   
   
   do j=1,set(i)%nats
      

      set(i)%x(j,:)=set(i)%x(j,:)-shifter(i,:)
     
   end do

  end do
  
  end subroutine shift_geom

 end module geom_class
