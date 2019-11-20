! this module provides parameters that choose the correct kind of REAL,COMPLEX and INTEGER variables 


module precision
  integer, parameter   :: sint=selected_int_kind(9)    ! integers up to 10^9 usually 4 byte (MPI_INTEGER4)
  integer, parameter   :: lint=selected_int_kind(19)   ! integers up to 10^19 usually 8 byte (MPI_INTEGER8) 

  integer, parameter   :: spreal=selected_real_kind(6)   ! REAL/COMPLEX with 7 digits usually 4 byte (MPI_REAL4,MPI_COMPLEX8)
  integer, parameter   :: dpreal=selected_real_kind(12)  ! REAL/COMPLEX with 12 digits usually 8 byte (MPI_REAL8,MPI_COMPLEX16) 

  integer, parameter   :: maxinteger=2147483647        ! maximum single precision integer 2**(31)-1
  
  integer, parameter   :: nprint=0   ! controls the details printed out, production run = 0
  logical,public       :: print_infor = .true.
  real(dpreal),public :: eps_sp,eps_dp
  real(dpreal),parameter :: pi = acos(-1.0_dpreal)           ! pi and imaginary unit 
  complex(dpreal),parameter :: iu = (0.0_dpreal,1.0_dpreal)


 ! this subroutine initializes the machine precision for sp and dp



end module precision
