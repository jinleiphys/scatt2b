module constants
 use precision
 implicit none



 real(dpreal),parameter :: hbarc=197.3269718_dpreal         ! NIST Ref 02.12.2014   ! MeV.fm           
 real(dpreal),parameter :: finec=137.03599_dpreal
 real(dpreal),parameter :: amu=931.49432_dpreal  !MeV
 real(dpreal),parameter :: e2=1.43997_dpreal         ! MeV.fm
 real(dpreal),parameter :: zero=0.0_dpreal 
 real(dpreal),parameter :: convert=1.0_dpreal    ! convert integer to real
 
CONTAINS

SUBROUTINE printconstants
 IMPLICIT NONE
 
 write(*,*)'***************** CONSTANTS **********************'
 write(*,'(" * hbarc=",1f9.5," MeV.fm",5x, "e^2=",1f8.5," MeV.fm *")') hbarc,e2
 write(*,'(" * so, alpha= 1/",1f9.5,$)')hbarc/e2
 write(*,'(7x,"amu=",1f8.4, " MeV  *")') amu
 write(*,*)'**************************************************'

 
END SUBROUTINE printconstants
end module constants

