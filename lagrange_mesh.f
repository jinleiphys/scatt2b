      module  lagrange_mesh_single_channel
      use mesh     
      implicit none 
      ! note the delta function structure fi(xj) = \lambda_i^{-1/2} \delta_{ij}
      real*8,dimension(:),allocatable :: lagrange_func 
      real*8,dimension(:,:), allocatable :: kinetic_matrix
      real*8,dimension(:,:),allocatable :: l_barrier_matrix 
      complex*16,dimension(:,:),allocatable :: V_matrix
      real*8,dimension(:),allocatable ::  mesh_rr, mesh_rw
      real*8 :: rmax_rmatrix
      contains
      
c***********************************************************************      
       subroutine initial_lagrange_func (rmatch_rmatrix)
       !
       ! before calling this subroutine please initial rmax_rmatrix and nr
       !
       ! this subroutine is used to initial the Lagrange function 
       ! for the purpose to reconstruct the scattering wave function 
       
       !!!!!! please gives the matching radius 
c***********************************************************************  
       use gauss
       use mesh 
       use precision
       implicit none 
       integer :: ir 
       real*8 :: rmatch_rmatrix
       rmax_rmatrix=rmatch_rmatrix
       
       if(allocated(lagrange_func)) deallocate(lagrange_func)
       if(allocated(mesh_rr)) deallocate(mesh_rr)
       if(allocated(mesh_rw)) deallocate(mesh_rw)
       
       allocate(lagrange_func(1:nr))
       allocate(mesh_rr(1:nr))
       allocate(mesh_rw(1:nr))
       
       call gauleg(nr,0.0_dpreal,1.0_dpreal,mesh_rr,mesh_rw)
       
       do ir=1, nr 
         lagrange_func(ir) = 1/ sqrt(rmax_rmatrix * mesh_rw(ir))
       end do 
       
       end subroutine 
c-----------------------------------------------------------------------
c***********************************************************************      
       subroutine T_and_Bloch(mu)
       ! this subroutine is used to compute the kinetic and Bloch operator
       ! in the matrix form
c***********************************************************************      
       use precision
       use constants
       implicit none 
       integer :: ir, irp
       real*8 :: f1,f2,f3
       real*8 ::  xi, xj
       real*8 :: mu
       
       
       if(allocated(kinetic_matrix)) deallocate(kinetic_matrix)
       allocate(kinetic_matrix(1:nr, 1:nr))
       
       do ir=1 , nr 
         do irp=1, nr
         
            f1=0.0_dpreal
            f2=0.0_dpreal
            f3=0.0_dpreal
            xi= mesh_rr(ir)
            xj= mesh_rr(irp)
            
            if(ir==irp) then 
               f1= hbarc**2 / (6.0_dpreal * rmax_rmatrix**2 * xi * (1-xi) * mu) 
               f2= 4.0_dpreal * nr * (nr+1) + 3.0
               f3= (1.0_dpreal - 6.0_dpreal *xi )/( xi * (1-xi) )
               kinetic_matrix(ir,irp)= f1 * ( f2 + f3 )
            else 
               f1= (-1.0_dpreal)**(ir+irp) * hbarc**2 /2.0_dpreal/rmax_rmatrix/rmax_rmatrix/
     +             sqrt( xi * xj * (1-xi) * (1-xj) ) / mu
               f2= nr * (nr+1)  + 1.0 + ( xi + xj - 2.0_dpreal * xi * xj )/( (xi-xj)**2 )
               f3= 1.0_dpreal/(1.0_dpreal-xi)  + 1.0_dpreal/(1.0_dpreal -xj)
               

               kinetic_matrix(ir,irp)= f1* (f2-f3)
            end if          
         end do
       end do     
       
 
       
       end subroutine
c-----------------------------------------------------------------------
c***********************************************************************      
       subroutine centrifugal_barrier (l,mu)
       ! this subroutine is used to compute centrifugal_barrier in 
       ! lagrange mesh basis 
       ! V_l = \frac{\hbarc^2 l (l+1)}{2 * mu * r}
       ! input : mu  in MeV
       !         l 
c***********************************************************************     
       use constants 
       use precision
       implicit none 
       integer :: l , ir 
       real*8 :: mu ,xi 
       if(allocated(l_barrier_matrix)) deallocate(l_barrier_matrix)
       allocate(l_barrier_matrix(1:nr, 1:nr))
       l_barrier_matrix=0.0_dpreal
       
       do ir=1 , nr 
         xi = mesh_rr(ir)
         l_barrier_matrix(ir,ir) = ( hbarc**2 * l * (l+1) )/ ( 2.0_dpreal * mu * rmax_rmatrix * rmax_rmatrix * xi * xi  )
       end do 
       
       end subroutine
c-----------------------------------------------------------------------
c***********************************************************************      
       subroutine lagrange_V(Vpot)
       ! this subroutine is used to compute the V-matrix in the Lagrange 
       ! mesh basis 
       ! input V is given in the mesh point in the step of hcm 
       ! this subroutine interpolate the Vpot to give the correct mesh sets
c***********************************************************************      
       use interpolation
       use precision
       implicit none 
       complex*16,dimension(0:irmatch),intent(in) :: vpot
       integer :: ir 
       real*8 ::  r 
       
       if(allocated(V_matrix)) deallocate(V_matrix)
       allocate(V_matrix(1:nr,1:nr))
       V_matrix=0.0_dpreal
       
       do ir = 1, nr 
         r = rmax_rmatrix*  mesh_rr(ir)
         V_matrix(ir,ir) = FFC(r/hcm, vpot, 1+irmatch)
       end do 

       
       end subroutine
c-----------------------------------------------------------------------
c***********************************************************************      
       subroutine R_matrix(l,mu,ecm,vpot,cph,gc,gcp,fc,fcp)
       ! this subroutine is used to compute the single channel R-matrix 
c***********************************************************************   
       use precision  
       use constants 
       implicit none
       integer :: l 
       real*8 :: mu, ecm ,k 
       real*8 :: gc,gcp,fc,fcp
       complex*16,dimension(0:irmatch) :: vpot
       complex*16,dimension(1:nr,1:nr) :: C_minus_E
       complex*16,dimension(1:nr) ::  B_vector , X_vector
       integer :: ir,INFO,NRHS
       integer,dimension(1:nr) :: IPIV
       real*8 :: N
       complex*16 ::  Rmatrix , Zmatrix , Smatrix
       complex*16 :: hc ,hc1 !H(+), H(-)
       complex*16 :: hcp ! derivatives of H(+)
       complex*16,dimension(1:nr) :: scattwf 
       complex*16 :: wf_bound
       real*8 ::  cph
       
       
       k=sqrt(2.*mu*ecm/(hbarc**2))
       hc=cmplx(gc,fc,kind=8)   
       hc1=cmplx(gc,-fc,kind=8)
       hcp=cmplx(gcp,fcp,kind=8)

       
       call lagrange_V(Vpot)
       call centrifugal_barrier (l,mu)

       C_minus_E=0.0_dpreal 
       B_vector=0.0_dpreal
       C_minus_E = kinetic_matrix + l_barrier_matrix + V_matrix
       
       do ir=1, nr 
          C_minus_E(ir,ir) =  C_minus_E(ir,ir) - ecm 
          B_vector(ir) = (-1.0_dpreal)**ir / sqrt( rmax_rmatrix* mesh_rr(ir) * ( 1.0-mesh_rr(ir) ) )
       end do 
       
       NRHS=1
       X_vector= B_vector
       
       call ZGESV( nr, NRHS, C_minus_E, nr, IPIV, X_vector, nr, INFO ) 
       If(INFO/=0) stop "error in calling ZGESV" 
       
       
       ! normalize factor 
       N = hbarc**2 / (2.0_dpreal * mu * rmax_rmatrix)
       
       Rmatrix=0.0_dpreal 
       do ir=1,nr 

            Rmatrix = Rmatrix + B_vector(ir) * X_vector(ir)
       end do 

       Rmatrix=Rmatrix * N
       Zmatrix= (hc-k*rmax_rmatrix*Rmatrix*hcp) / (k*rmax_rmatrix)
       
       Smatrix= CONJG(Zmatrix) / Zmatrix
       
       write(*,*) "Smatrix=", Smatrix, 0.5_dpreal*log(Smatrix)/iu
       write(*,*) "phase=",real(0.5_dpreal*log(Smatrix)/iu) * 180.0_dpreal /pi
       
       
       ! compute the scattering wave function  ! still testing 
       
       wf_bound = 0.5* iu * (hc1-hc*Smatrix) * exp(iu*cph)
       do ir=1,nr
       
        scattwf(ir)=lagrange_func(ir)*X_vector(ir)*wf_bound *N / Rmatrix  
       
       write(99,*) mesh_rr(ir)* rmax_rmatrix , real(scattwf(ir)), aimag(scattwf(ir))
       end do 
       
       
       
C      stop
       
       end subroutine
c-----------------------------------------------------------------------      
      end module