      module scattwf
      implicit none
      complex*16,dimension(:,:),allocatable :: wf ! wave  function
      real*8,dimension(:),allocatable :: nfc,ngc,nfcp,ngcp ! used for coul90
      real*8,dimension(:,:),allocatable ::  pl


      contains
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine scatt2b()
c     calculate the scattering wave function for a given two body system
c     the coupling coffcient for this two body problem
c     should be | l (jp jt)s  ; J M>
c     sp : spin of projectile
c     st : spin of target
c     l  : angular momentum between projectile and target
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       use mesh,only:irmatch,hcm
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

       use coulfunc
       ! Coulomb/Bessel function

       use lagrange_mesh_single_channel
       implicit none

       integer :: l,s,j,nch ! channel index
       integer :: ifail ! for subroutine coul90
       real*8 :: mu ! reduce mass in MeV
       real*8 :: k ! wavenumber in fm^{-1}
       real*8 :: rho ! dimensionless factor for r * k used in coul90
       real*8 :: ecm ! C.M. frame bombing energy in MeV
       real*8 :: eta ! Sommerfeld parameter
       complex*16 :: sl,nl !  s-matrix and normalization parameter
       real*8,dimension(0:lmax) :: cph  !Coulomb phase-shift
       integer :: r0 ! starting index for solving the differential equations r0= 2*l
       real*8 :: ls ! ls=0.5_dpreal*(j*(j+1)-l*(l+1)-s*(s+1))
       complex*16,dimension(1:5) :: wfmatch ! 5 points to match the Coulomb/Bessel function
       complex*16,dimension(:,:),allocatable :: Upot ! potential used for calculations
       complex*16, dimension(0:lmax) ::  smat


       if(allocated(wf)) deallocate(wf) ! deallocate the wave function
       if(allocated(Upot)) deallocate(Upot) ! deallocate the potential used in the calculation


       allocate(nfc(0:lmax),ngc(0:lmax),nfcp(0:lmax),ngcp(0:lmax))
       allocate(wf(0:irmatch,1:alpha2b%nchmax))
       allocate(Upot(0:irmatch,1:alpha2b%nchmax))


       ecm=elab*masst/(massp+masst) ! compute the energy in C.M. frame (MeV)
       mu=amu*(masst*massp)/(massp+masst) ! reduced mass  (MeV)
       k=sqrt(2.*mu*ecm/(hbarc**2)) ! wavenumber (fm^{-1})
       rho=(irmatch-2)*hcm*k        ! used for coul90
       eta=zp*zt*e2*mu/hbarc/hbarc/k ! Sommerfeld parameter
       write(*,*) "eta=",eta


       ! test lagrange mesh
       call initial_lagrange_func(rmax-hcm*2.)
       call T_and_Bloch(mu)
       !

       ! check the step size
       if (k*hcm>0.2)  then
         write(*,*) 'warning!please decrease the value of hcm,',
     +     'it should be smaller than ', 0.2/k
C        stop
       end if

       ! compute the Coulomb/Bessel function used for matching
       ! if eta = 0 returning Bessel function
       ! if eta /= 0 returning Coulomb function
       call coul90(rho,eta,zero,lmax,nfc,ngc,nfcp,ngcp,0,ifail)
       if (ifail/=0) then
       write(*,*) 'coul90: ifail=',ifail; stop
       endif

       ! compute the Coulomb phase-shift
       call coulph(eta,cph,lmax)
       write(*,*)"elastic S-matrix elements for the incoming channel"

       ! solve the differential equations for each channel
       do nch=1,alpha2b%nchmax

          l=alpha2b%l(nch)
          s=alpha2b%s(nch)
          j=alpha2b%j(nch)
          ls=0.5_dpreal*(j*(j+1)-l*(l+1)-s*(s+1))

          ! obtain the potential
          call potr(zp*zt,ls)
          Upot(0:irmatch,nch)=v

          ! solve the differential equation
          r0=2*l
            call sch_enhanced_numerov(r0,mu,ecm,Upot(0:irmatch,nch),l,wf(:,nch))
C         call sch_numerov(r0,mu,ecm,Upot(0:irmatch,nch),l,wf(:,nch))

          ! matching the boundary conditions
          wfmatch(1:5)=wf(irmatch-4:irmatch,nch)
          call matching(l,k,wfmatch,sl,nl)
          smat(l)=sl
          ! compute the phase-shift
          write(2,*) l, 0.5_dpreal*log(sl)/iu, real(0.5_dpreal*log(sl)/iu) * 180.0_dpreal /pi
          ! renormalize the wave function
          wf(:,nch)=wf(:,nch)*nl
          wf(:,nch)=wf(:,nch)*exp(iu*cph(l))
          call chan_out(nch,wf(0:irmatch,nch),sl)


          ! test lagrange mesh
          call R_matrix(l,mu,ecm,Upot,cph(l),ngc(l),ngcp(l),nfc(l),nfcp(l))
          !

       end do

       call plcos()
       call angular_distribution(cph,smat,k,eta)
C      call plotwfxz()
       write(*,150)
150   format('********************************************************')

      end subroutine
c-----------------------------------------------------------------------
      subroutine plotwfxz()
c     plot the |\chi(\vec{r})| in xz plane
       use mesh
       use channels
       ! channel index
       use precision
       use constants
       use interpolation
       use systems,only:zp,massp,zt,masst,elab
       !  projectile mass number and charge number: massp, zp
       !  the target mass number and charge number: masst, zt
       !  incoming energy in lab frame : unit in MeV
       use pyplot_module, only : pyplot
       implicit none
       complex*16,allocatable,dimension(:,:) :: sumwf
       real*8,allocatable,dimension(:) :: xx
       integer :: ith,ir,l,nch, ix , iz , nxmax
       real*8 :: x,z,theta,r,k,mu,ecm
       integer :: istat
       type(pyplot)             :: plt   !! pytplot handler

       ecm=elab*masst/(massp+masst) ! compute the energy in C.M. frame (MeV)
       mu=amu*(masst*massp)/(massp+masst) ! reduced mass  (MeV)
       k=sqrt(2.*mu*ecm/(hbarc**2)) ! wavenumber (fm^{-1})
       nxmax=nint(irmatch/2.)
       allocate(sumwf(-nxmax/5:nxmax/5,-nxmax/5:nxmax/5))
       allocate(xx(-nxmax/5:nxmax/5))
       sumwf=0.0_dpreal
      do ix=-nxmax/5, nxmax/5
         x=ix*hcm*5
         xx(ix)=x
         do iz=-nxmax/5, nxmax/5
            z=iz*hcm*5
            r=sqrt(x**2 + z**2)
            theta=acos(z/r)
            do nch=1,alpha2b%nchmax
               l=alpha2b%l(nch)
c               sumwf=sumwf+ 2. * sqrt(pi*(2.*l+1.))* iu**l * wf(ir,nch) * pl(l,ith) / k / r
               sumwf(iz,ix)=sumwf(iz,ix)+ 2. * sqrt(pi*(2.*l+1.))* iu**l * FFC(r/hcm,wf(:,nch),irmatch+1) *
     &          FFR4((theta*180./pi-thmin)/thinc,pl(l,:),nth) / k / r

            end do
         end do
      end do

C      do iz=-nxmax,nxmax,5
C           z=iz*hcm
C           write(765,*)( abs(sumwf(iz,ix)),ix=-nxmax,nxmax)
C           write(766,*) iz*hcm
C           do ix=-nxmax, nxmax,5
C             write(769,*) iz*hcm, ix*hcm ,abs(sumwf(iz,ix))
C             write(770,*) iz*hcm, ix*hcm ,abs( sumwf(iz,ix)-exp(iu*k*z) )
C             write(771,*) iz*hcm, ix*hcm ,abs( exp(iu*k*z) )
C           end do
C      end do


      !contour plot:
      call plt%initialize(grid=.false.,xlabel='Z (fm)',
     &                      ylabel='X (fm)',figsize=[10,10],
     &                      title='2dplot', real_fmt='*',
     &                      axisbelow=.false.)
      call plt%add_contour(xx, xx, abs(sumwf), linestyle='-',
     &                      filled=.true., cmap='bone', colorbar=.true.,
     &                     istat=istat)
      call plt%savefig('contour.png',pyfile='contour.py',istat=istat)





      end subroutine


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sch_enhanced_numerov(r0,mu,ecm,vpot,l,rwfl)
c     mu! reduce mass
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       use mesh, only: irmatch,hcm
       use constants, only:hbarc,e2
       use precision
       implicit none
       integer,intent(in) :: r0
       real*8,intent(in) :: ecm ! energy
       complex*16,dimension(0:irmatch),intent(in) :: vpot
       integer,intent(in) :: l
       real*8,intent(in) :: mu  ! reduced mass
       integer :: ir
       real*8 :: r
       complex*16,dimension(0:irmatch),intent(out) :: rwfl   ! partial radial wave function
       complex*16,dimension(0:irmatch) :: kl
       complex*16,dimension(0:irmatch) :: Tx
       complex*16,dimension(0:irmatch) :: Wx


        kl=0.0d0
        rwfl=0.0d0
        Tx=0.0d0


c      Numerov method to solve the differential equation
       rwfl(r0)=0      ! boundary condition


       ir=r0+1; r=ir*hcm
       rwfl(ir)=hcm**(l+1)  ! arbitrary value
       if(l>140) rwfl(r0+1)=0.000000000000000000001


       kl(ir)=2.*mu*ecm/hbarc**2-l*(l+1)/r**2-2.*mu*Vpot(ir)/hbarc**2
       Tx(ir)=-hcm**2/12.0d0*kl(ir)
       Wx(ir)=(1-Tx(ir))*rwfl(ir)

       rwfl(r0+2)=2.*rwfl(r0+1)-hcm**2*kl(r0+1)*rwfl(r0+1)
       ir=r0+2; r=ir*hcm
       kl(ir)=2.*mu*ecm/hbarc**2-l*(l+1.)/r**2-2.*mu*Vpot(ir)/hbarc**2
       Tx(ir)=-hcm**2/12.0d0*kl(ir)
       Wx(ir)=(1-Tx(ir))*rwfl(ir)


       do ir=r0+2 ,irmatch-1

        kl(ir+1)=2.*mu*ecm/hbarc**2-l*(l+1.)/((ir+1.)*hcm)**2-2.*mu*Vpot(ir+1)/hbarc**2
        Tx(ir+1)=-hcm**2/12.0d0*kl(ir+1)
        Wx(ir+1)=(2+12.*Tx(ir)+12.*Tx(ir)**2)*Wx(ir)-Wx(ir-1)
        rwfl(ir+1)=Wx(ir+1)/(1.-Tx(ir+1))

       end do


      end subroutine

c----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sch_numerov(r0,mu,ecm,vpot,l,rwfl)
c     mu! reduce mass
c     Numerov method to solve the differential equation
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       use mesh, only: irmatch,hcm
       use constants, only:hbarc,e2
       use precision
       implicit none
       integer,intent(in) :: r0
       real*8,intent(in) :: ecm ! energy
       real*8,intent(in) :: mu  ! reduced mass
       complex*16,dimension(0:irmatch),intent(in) :: vpot
       integer,intent(in) :: l
       complex*16,dimension(0:irmatch) :: rwfl   ! partial radial wave function
       real*8 :: r
       real*8 :: const
       integer :: ir
       complex*16,dimension(1:irmatch) :: kl

       rwfl=0.0_dpreal
c      Numerov method to solve the differential radial equation
       rwfl(r0)=0      ! boundary condition

       rwfl(r0+1)=hcm**(l+1)  ! arbitrary value
       !if(l>140) rwfl(r0+1)=0.000000000000000000001

       ir=r0+1; r=ir*hcm
       kl(ir)=2.*mu*ecm/hbarc**2-l*(l+1)/r**2-2.*mu*Vpot(ir)/hbarc**2

       rwfl(r0+2)=2.*rwfl(r0+1)-hcm**2*kl(r0+1)*rwfl(r0+1)
       ir=r0+2; r=ir*hcm
       const=hcm**2/12.
       kl(ir)=2.*mu*ecm/hbarc**2-l*(l+1.)/r**2-2.*mu*Vpot(ir)/hbarc**2


       do ir=r0+2 ,irmatch-1
        kl(ir+1)=2.*mu*ecm/hbarc**2-l*(l+1.)/((ir+1.)*hcm)**2-
     &                  2.*mu*Vpot(ir+1)/hbarc**2

        rwfl(ir+1)=((2.-10.*const*kl(ir))*rwfl(ir)-(1.+const*kl(ir-1))
     &             *rwfl(ir-1))/(1.+const*kl(ir+1))

       end do


c****


      end subroutine

c----------------------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine matching(l,k,wf,sl,nl) !wf has dimension (5)
c     nl*wf=0.5*i*(H(-)-sl*H(+))
c     nl*wfp=0.5*i*k*(H'(-)-sl*H'(+))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       use precision, only:pi,iu
       use mesh,only:hcm
       implicit none

       integer,intent(in) :: l
       real*8,intent(in) :: k
       complex*16,intent(in),dimension(1:5) :: wf ! wavefunction at rmatch
       complex*16 :: wfp !derivative of wf
       complex*16,intent(out) :: nl !Normalization parameter
       complex*16 :: hc,hc1  !H(+),H(-)
       complex*16 ::hcp,hcp1 ! derivatives of H(+),H(-)
       complex*16,intent(out) :: sl ! S-matrix

       hc=cmplx(ngc(l),nfc(l),kind=8)
       hc1=cmplx(ngc(l),-nfc(l),kind=8)
       hcp=cmplx(ngcp(l),nfcp(l),kind=8)
       hcp1=cmplx(ngcp(l),-nfcp(l),kind=8)

       wfp=(-wf(5)+8.*wf(4)-8.*wf(2)+wf(1))/12./hcm
       nl=(hc1*hcp*iu*k-hc*hcp1*iu*k)/(2.*(hcp*wf(3)*k-hc*wfp))
       sl=(hc1*wfp-hcp1*wf(3)*k)/(hc*wfp-hcp*wf(3)*k)


      end subroutine matching
c----------------------------------------------------------------------



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine chan_out(nch,rwfl,sl)
c     output subroutine
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       use mesh,only:irmatch,hcm
       implicit none
       integer,intent(in) :: nch
       integer :: ir
       complex*16,dimension(0:irmatch) :: rwfl
       real*8,dimension(0:irmatch) :: rwflr,rwfli ! real part and imaginary part
       complex*16,intent(in) :: sl

       rwflr=real(rwfl)
       rwfli=aimag(rwfl)
c-----------------------------------------------------------------------
c*** print elastic S-matrix.
       write(3,*)"elastic S-matrix"
       write(3,100)real(sl),aimag(sl),nch
100    format(F13.9,2X,F13.9,2X,I3,2X)
c-----------------------------------------------------------------------

       write(*,101) nch,real(sl),aimag(sl)
101    format('nch=',I3,2X,'S-matrix = ('
     &           ,F13.9,',',F13.9,')')
c-----------------------------------------------------------------------
       write(4,102)nch
102    format('@nch=',I3)

       do ir=0,irmatch
          write (4,*) hcm*ir, rwflr(ir),rwfli(ir)
       end do
c----------------------------------------------------------------------
       write (4,*) "& "
      end subroutine

c----------------------------------------------------------------------
      subroutine plcos()
      use mesh
      use precision
      use channels
      implicit none
      integer ::  ith, l
      real*8 :: theta



      if (.not. allocated(pl)) allocate( pl(0:lmax,1:nth) )
      pl=0.0_dpreal

      do ith=1,nth
        theta=thmin+ thinc*(ith-1)
        pl(0,ith) = 1.0_dpreal
        pl(1,ith) = cos(theta*pi/180.)
        do l=2,lmax
          pl(l,ith)=dble(2.*l-1.)/dble(l)*cos(theta*pi/180.)*pl(l-1,ith)-dble(l-1.)/dble(l)*pl(l-2,ith)
        end do
      end do

      end subroutine



c-----------------------------------------------------------------------
      subroutine angular_distribution(cph,smat,k,eta)
      use mesh
      use precision
      use channels
      use systems
      implicit none
      complex*16, dimension(0:lmax) ::  smat
      real*8,dimension(0:lmax) :: cph
      integer :: l
      real*8 :: dsdw, dsdw_Ruth
      integer :: ith
      complex*16 ::  fc, fn , f
      real*8 :: theta,k,eta,theta_radian, sin2theta



      do ith=1, nth
         theta=thmin+ thinc*(ith-1)
         theta_radian= theta* pi / 180.
         sin2theta = sin( 0.5*theta_radian )**2
         fc=-eta * exp( -iu*eta*log(sin2theta) + 2*iu*cph(0) ) / ( 2.0_dpreal*k*sin2theta )

         fn=0.0_dpreal
         do l=0, lmax
           fn=fn + (2.*l+1.)*exp(2*iu*cph(l))*(smat(l)-1.)*pl(l,ith) / (2.*iu*k)
         end do
         f= fc+fn
         dsdw= abs(f) ** 2 * 10
         dsdw_Ruth = eta**2 * 10 / (4 * k**2 * sin2theta**2 )
          write(199,*) theta, dsdw, dsdw_Ruth
         if(abs(eta) > 1e-6) dsdw = dsdw/dsdw_Ruth

         write(99,*) theta, dsdw
      end do






      end subroutine
c--------------------------------------------------------------------------

      end module
