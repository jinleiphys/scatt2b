c----------------------------------------------------------------------
c   The full CENTRAL (local) potential will be assumed of the form
c
c    V(r)=Vcou(r)+ Vc(r)+Vsur(r)
c
c    Vcou(r) = coulomb central   !! real
c    Vc(r) =  woods-saxon potential   !! complex
c    Vsur(r) = surface potential   !complex
c----------------------------------------------------------------------
      module pot
C     implicit none
       complex*16,allocatable,dimension(:) :: V ! total potential
       real*8  :: uv,av,rv ! parameters of real part of W-S
       real*8  :: uw,aw,rw ! parameters of imaginary part of W-S
       real*8  :: vsov,rsov,asov ! real part spin-orbit potential for projectile
       real*8  :: vsow,rsow,asow ! imaginary part spin-orbit potential for projectile
       real*8  :: vd,avd,rvd  ! real surface part
       real*8  :: wd,awd,rwd  ! imaginary surface
       real*8  :: a1,a2    ! mass for radius conversion
       real*8  :: rc ! for coulomb potential
       integer :: ptype ! select potential type
       real*8, allocatable,dimension(:) :: rvec

      contains
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine potr(z12,ls)
c     subroutine to calculate the potential
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       use mesh,only:irmatch,hcm
       use precision
       implicit none
       real*8 :: a13,r !a1^0.333333+a2^0.333333
       character(len=1) :: kpot1
       integer :: kpot2
       complex*16 :: vc,vls,vsur
       real*8 :: wsv,wsw,vcou,vlsv,vlsw,vsurv,vsurw   !real part and imaginary part
       real*8 :: z12  ! For Coulomb potential
       real*8 ::ls !l*s
       integer :: ir
       real*8, dimension(0:irmatch) :: vr,vi
       character*40 filename

       vc=0.0d0;vcou=0.0d0;vls=0.0d0;vsur=0.0d0
       wsv=0.0d0;wsw=0.0d0;vcou=0.0d0;vlsv=0.0d0;vlsw=0.0d0
       vsurv=0.0d0;vsurw=0.0d0


      if (allocated(v)) deallocate(v)
      allocate(v(0:irmatch))
      v=0.0d0

c-----------------------------------------------------------------------
       a13=a1**(1./3.)+a2**(1./3.)
       if (a13<1e-4) then
       write(*,*) 'a1=',a1,'a2=',a2
       write(*,*) "please check input mass of potential"
c       stop
       end if
c-----------------------------------------------------------------------



c-----------------------------------------------------------------------
       if (ptype==41) then !read external potentials
c-----------------------------------------------------------------------
         if (allocated(rvec)) deallocate(rvec)
         allocate(rvec(0:irmatch))    !!! 0 or 1  ??
         do ir=0, irmatch
          rvec(ir)=ir*hcm
         end do
         filename='fort.41'
         call extpot(filename,vr,vi,irmatch)
         do ir=0,irmatch
           v(ir)=cmplx(uv*vr(ir),uw*vi(ir),kind=8)
         end do
       end if 
c-----------------------------------------------------------------------




c      calculate the potential
       do ir=1,irmatch
       r=ir*hcm
       vc=0.0d0
       vls=0.0d0
       vsur=0.0d0
       vcou=0.0d0
       select case(ptype)
c-----------------------------------------------------------------------
       case(1) ! Woods-Saxon
c-----------------------------------------------------------------------
       !central potential
       wsv=-ws(r,uv,rv*a13,av)
       wsw=-ws(r,uw,rw*a13,aw)
       vc=cmplx(wsv,wsw,kind=8)
       
       ! spin-orbit potential
       vlsv=wsso(r,vsov,rsov*a13,asov)*ls
       vlsw=wsso(r,vsow,rsow*a13,asow)*ls
       vls=cmplx(vlsv,vlsw,kind=8)
       
       ! surface potential
       vsurv=4*avd*dws(r,vd,rvd*a13,avd)
       vsurw=4*awd*dws(r,wd,rwd*a13,awd)
       vsur=cmplx(vsurv,vsurw,kind=8)
       
       ! Coulomb potential
       vcou=vcoul(r,z12,rc*a13)
c-----------------------------------------------------------------------
       case(2) !Gaussian
c-----------------------------------------------------------------------
       ! central pot
       wsv=-gausspot(r,uv,rv*a13,av)
       wsw=-gausspot(r,uw,rw*a13,aw)
       vc=cmplx(wsv,wsw,kind=8)
       
       ! spin - orbit pot
       vlsv=gausder(r,vsov,rsov*a13,asov)*ls
       vlsw=gausder(r,vsow,rsow*a13,asow)*ls
       vls=cmplx(vlsv,vlsw,kind=8)
        
       ! Coulomb potential
       vcou=vcoul(r,z12,rc*a13)
c-----------------------------------------------------------------------
       case(3)
c-----------------------------------------------------------------------
       vc=MT(r)
       
c-----------------------------------------------------------------------
       case(41) !read potential
c-----------------------------------------------------------------------
       vsurv=4*avd*dws(r,vd,rvd*a13,avd)
       vsurw=4*awd*dws(r,wd,rwd*a13,awd)
       vsur=cmplx(vsurv,vsurw,kind=8)
       vcou=vcoul(r,z12,rc*a13)
       end select

       v(ir)=vc+vcou+vls+vsur+v(ir)

       end do


       
       do ir=1,irmatch
          write(30,*) ir*hcm,real(v(ir)),aimag(v(ir))
       end do
       write(30,*) "&"


      end subroutine potr


c-----------------------------------------------------------------------


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c             functions of different types of potential
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c *** Woods-Saxon (volume)
       function ws(r,v0,r0,a)
       implicit none
       real*8 r,v0,r0,a,ws,aux
       ws=0d0
        if (abs(v0).lt.1e-6) return
        if (a>1e-6) then
           aux=exp(-(r-r0)/a)
         ws=v0/(1d0+1d0/aux)
        else
         write(0,*)'WS: a too small!'
        end if
        return
        end function
c-----------------------------------------------------------------------
c *** Gaussian
      function gausspot(r,v0,r0,a)
      implicit none
       real*8 r,v0,r0,gausspot,a
         if (a.gt.1e-6) then
           gausspot=V0*exp(-(r-r0)**2/a**2)
             else
               write(*,*)'a too small in gausspot!'
               stop
         endif
         return
      end function

c-----------------------------------------------------------------------
c Coulomb potential
      FUNCTION VCOUL(R,z12,Rc)
          use constants
          implicit none
          real*8 r,rc,rc2,aux,vcoul,z12

          RC2=RC*2d0
          aux=e2*Z12
          vcoul=0
          if (z12.lt.1e-4) return
          if (rc.lt.1e-6) rc=1e-6

          IF(R.GT.RC)GO TO 1
          VCOUL=AUX*(3.-(R/RC)**2)/RC2
          RETURN
1         VCOUL=AUX/R
          RETURN
        END

c-----------------------------------------------------------------------
c *** Spin-orbit with WS derivative
c     This form factor will be then multiplied by l.s
c     We use fresco definiton for spin-orbit
      function wsso(r,v0,r0,a)
        implicit none
        real*8 r,v0,r0,a,wsso,aux,conls
        parameter(conls=2.0)
        wsso=0.0d0
         if (r<1e-6) r=1e-6
         if (a>1e-6) then
          aux=exp(-(r-r0)/a)
          wsso=-2d0*conls*v0*aux/(1d0+aux)**2/a/r
         else
          write(0,*)'WS spin-orbit : a too small!'
         endif
         return
      end function

c-----------------------------------------------------------------------
c *** Spin-orbit with Gauss derivative
c     This form factor will be then multiplied by l.s
c     We use fresco definiton for spin-orbit
        function gausder(r,v0,r0,a)
      implicit none
      real*8 r,v0,r0,a,gausder,conls,rh
      parameter(conls=2.0)
      gausder=0.0d0
       if (r<1e-6) r=1e-6
         if (a>1e-6) then
            rh=(r-r0)/a
            gausder=-exp(-rh**2)**rh*conls*v0/(r*a)
         else
           write(0,*)'Gauss spin-orbit : a too small!'
       endif
       return
      end function

c-----------------------------------------------------------------------
c *** WS derivative
      function dws(r,v0,r0,a)
      implicit none
      real*8 r,v0,r0,a,dws,aux
        if (r<1e-6) r=1e-6
           if (a>1e-6) then
             aux=exp(-(r-r0)/a)
             dws=-v0*aux/(1d0+aux)**2/a
           else
             write(0,*)'derivative WS: a too small!';stop
         endif
         return
      end function
c-----------------------------------------------------------------------
c** Malfliet-Tjon potential 

       function MT(r)
       implicit none 
       real*8 :: r 
       real*8 :: VA, VR, muA, muR
       complex*16 :: MT
       
       VA= 626.8932
       VR= 1438.7228
       muA=1.550
       muR=3.11
       
       
       MT= dcmplx(VR*exp(-muR*r)/r - VA*exp(-muA*r)/r,2.0)
       
       
       end function

c-----------------------------------------------------------------------
c Read external potential
c Format:
c 1st line: header
c 2nd line: npoints, rstep, rfirst
c Next lines: real, imag
      subroutine extpot(filename,vr,vi,nr)
        use interpolation
        implicit none
        logical uu
        integer npt,n,nr
        real*8 vr(0:nr),vi(0:nr)
        integer,parameter:: kpot=50
        character*40 filename,header
        real*8:: r,rstep,rfirst,x,y1,y2
        real*8,allocatable::xv(:),faux(:),haux(:)
        real*8,parameter:: alpha=0

        uu = .false.
        inquire(file=filename,exist=uu)
        if (.not.uu) then
          write(0,*) 'Potential file:', filename,' not found!'
          stop
        endif
        write(*,'(8x, "Potential file:",a20)') filename
        open(kpot,file=filename,status='old')
        rewind(kpot)
        read(kpot,*) header
        npt=0
20      read(kpot,*,end=200) x,y1,y2
        npt=npt+1
        goto 20

200     if (npt.lt.1) goto 300
        rewind(kpot)

        read(kpot,*) header
        read(kpot,*) npt,rstep,rfirst
        allocate(faux(npt),haux(npt),xv(npt))
        faux(:)=0; haux(:)=0; xv(:)=0
        do n=1,npt
           read(kpot,*,err=300)  y1,y2
           xv(n)=rfirst+rstep*(n-1)
           faux(n)=y1
           haux(n)=y2
        enddo

      write(*,'(3x,"=> read:",i3," points")') npt
      write(*,250) xv(1),xv(npt),xv(2)-xv(1)
250   format(/,5x,"[Radial grid: Rmin=",1f7.3," fm,",
     &      " Rmax=",1f7.1," fm,"
     &      " Step=",1f7.3," fm]",/)

        do n=0,nr
          r=rvec(n)
          vr(n)=fival(r,xv,faux,npt,alpha)
          vi(n)=fival(r,xv,haux,npt,alpha)
        enddo

        deallocate(xv,faux,haux)
        return
300     write(*,*)'error reading ext. potential !'
        stop
        return
      end subroutine



      end module
