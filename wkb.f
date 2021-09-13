      module wkb 
!     calculating the wkb wave function by using the formula suggested 
!     by Nucl Phys A 259(1976) 99-121, Eq(3.1)
      implicit none
      
      
      
      contains 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      subroutine wkbwf()
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use mesh,only:irmatch,hcm,rmatch,nr
       !  matching radius irmatch*hcm 
       !  irmatch: index
       !  hcm : step (unite in fm)
       
       use channels
       ! channel index 
       
       use systems,only:zp,massp,zt,masst,elab
       !  projectile mass number and charge number: massp, zp 
       !  the target mass number and charge number: masst, zt
       !  incoming energy in lab frame : unit in MeV 
       
       use constants,only:amu,hbarc,e2,zero
       ! constants used in the calculation 
       ! amu:  atomic mass unit (MeV)
       ! hbarc :  hbar * c 
       ! e2 : e^2 charge unit square
       ! zero : numerical 0 in double precision
       
       use precision
       !  double precision for different compiler 
       
       use pot
       ! potentials use in the scattering calculation 
       
       use interpolation, only: ffc
       
       use gauss, only: gauleg
       implicit none
       integer :: l,s,j,nch ! channel index 
       integer :: ifail ! for subroutine coul90
       real*8 :: mu ! reduce mass in MeV 
       real*8 :: k ! wavenumber in fm^{-1}
       real*8 :: ecm ! C.M. frame bombing energy in MeV
       real*8 :: ls ! ls=0.5_dpreal*(j*(j+1)-l*(l+1)-s*(s+1))
       real*8 :: r0 ! classical turning point 
       complex*16,dimension(0:irmatch) :: wf ! wave  function
       integer :: ir, i
       real*8 :: eta, b ,r 
       real*8, dimension(1:nr) :: rr, rrw
       complex*16 :: f1, kr, krr, Veff, Veffp
       
      
       
       ecm=elab*masst/(massp+masst) ! compute the energy in C.M. frame (MeV)
       mu=amu*(masst*massp)/(massp+masst) ! reduced mass  (MeV)
       k=sqrt(2.*mu*ecm/(hbarc**2)) ! wavenumber (fm^{-1})
       eta=zp*zt*e2*mu/hbarc/hbarc/k
    
       wf=0.0_dpreal
       
       do nch=1,alpha2b%nchmax

          l=alpha2b%l(nch)
          s=alpha2b%s(nch)
          j=alpha2b%j(nch)
          ls=0.5_dpreal*(j*(j+1)-l*(l+1)-s*(s+1))
          
          ! obtain the potential
          call potr(zp*zt,ls)
C         v=0.0_dpreal
       
          ! find the classical turning point 
          b=sqrt(l*(l+1.0_dpreal))/k
          r0=eta/k+sqrt((eta/k)**2+b**2)
          
C         r0=0
          write(*,*) "Classical turning point r0=", r0
          write(*,*) "mu=",mu
          
          write(14,*) "&l=",l
          write(31,*)  "&l=",l

          do ir=1, irmatch
          
            r=ir*hcm
            Veff=hbarc**2 * (l+0.5_dpreal)**2 / 2.0_dpreal / mu / r**2 + v(ir)
            
             write(31,*) r, real(veff)
            
            if (r<r0) then 
               call gauleg(nr,r,r0,rr,rrw)
               kr=k*sqrt(Veff/Ecm-1.0_dpreal)          
               f1=0.0_dpreal  
               do i=1,nr
                 Veffp=hbarc**2 * (l+0.5_dpreal)**2 / 2.0_dpreal / mu / rr(i)**2 
     &                 + FFC(rr(i)/hcm,v(0:irmatch),irmatch+1) 
                 krr=k*sqrt(Veffp/Ecm-1.0_dpreal)
                 f1=f1+krr*rrw(i)         
               end do 
               wf(ir)=0.5*sqrt(k/kr)*exp(-f1)
            
            
            else 
               call gauleg(nr,r0,r,rr,rrw)
               kr=k*sqrt(1.0_dpreal-Veff/Ecm)
               f1=0.0_dpreal 
               do i=1, nr
                 Veffp=hbarc**2 * (l+0.5_dpreal)**2 / 2.0_dpreal / mu / rr(i)**2 
     &                 + FFC(rr(i)/hcm,v(0:irmatch),irmatch+1) 
                 krr=k*sqrt(1.0_dpreal-Veffp/Ecm)
                 f1=f1+krr*rrw(i)
               end do             
               wf(ir)=sqrt(k/kr) * sin(f1+0.25*pi)
            end if 
            
            write(14,*) r, real(wf(ir)), aimag(wf(ir))
          end do 
          
            
                    
          

       end do
             
      
      end subroutine 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      
      
      
      
      
      
      
      
      
      
      end module 