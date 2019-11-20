      module input
      implicit none
      logical,dimension(1:9999) ::  written
      contains
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine initialize()
c     parameter initialization
c       namelist /global/ hcm,rmatch,lmax,elab,lmin
c
c       namelist /system/ namep,massp,zp,jp,namet,masst,zt,jt
c                 					
c       namelist /potential/ ptype,a1,a2,rc,uv,av,
c                           rv,uw,aw,rw,vsov,rsov,asov,
c                           vsow,rsow,asow,vd,avd,rvd,wd,awd,rwd

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use systems
      use channels
      use mesh
      use pot
      use precision
      implicit none


      namelist /global/ hcm,rmatch,lmax,elab,lmin,nr

      namelist /system/ namep,massp,zp,jp,namet,masst,zt,jt


      namelist /potential/ ptype,a1,a2,rc,uv,av,
     &                      rv,uw,aw,rw,vsov,rsov,asov,
     &                      vsow,rsow,asow,vd,avd,rvd,wd,awd,rwd

       written=.false.
       written(1)=.true.;

C       open (unit=5,file='test.in')
c-----------------------------------------------------------------------
c /global/
       hcm=0.05_dpreal;rmatch=50.0_dpreal
       lmax=30;elab=10.0_dpreal
       lmin=0
       nr=80
       read(5,nml=global)
       irmatch=nint(rmatch/hcm)
       rmax=rmatch
c-----------------------------------------------------------------------
c/system/
       namep='null';massp=0.0d0;zp=0.0d0;jp=0.0d0
       namet='null';masst=0.0d0;zt=0.0d0;jt=0.0d0
       read(5,nml=system)
c-----------------------------------------------------------------------
c /potential/ parameter
c      if (ecmbh<0.000001) then

       ptype=1;a1=0.0d0;a2=0.0d0;rc=1.3d0
       uv=0.0d0;av=0.1d0;rv=0.0d0
       uw=0.0d0;aw=0.1d0;rw=0.0d0
       vsov=0.0d0;rsov=0.0d0;asov=0.1d0
       vsow=0.0d0;rsow=0.0d0;asow=0.1d0
       vd=0.0d0;avd=0.1d0;rvd=0.0d0
       wd=0.0d0;awd=0.1d0;rwd=0.0d0
       read(5,nml=potential)
c-----------------------------------------------------------------------
      end subroutine


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine check()
c     check the input file and give the local copy of input
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use systems
      use channels
      use mesh
      use pot
      use constants
      implicit none
      
     
      namelist /global/ hcm,rmatch,lmax,elab,lmin

      namelist /system/ namep,massp,zp,jp,namet,masst,zt,jt


      namelist /potential/ ptype,a1,a2,rc,
     &                      uv,av,rv,uw,aw,rw,vsov,rsov,asov,
     &                      vsow,rsow,asow,vd,avd,rvd,wd,awd,rwd
     
       write(1,nml=global)

       write(1,nml=system)

       write(1,nml=potential)


c-----------------------------------------------------------------------

c ***print parameters

      write(*,70) nint(rmatch/hcm),hcm,rmatch
70    format('Centre-of-mass Range is ',I5,'*',f5.3,' fm.',
     &       ' Maximum at ',f8.3,' fm.' )

      write(*,80)lmin,lmax
80    format('Range of angular momentum l is ',I1,' <= l <=',I3)

c***print reaction systems
      write(*,90)
90    format('*******************Reaction systems*********************')
      write(*,100) namep,massp,zp,jp
100   format('Project=',A5,' MASS = ',F7.4,' Z = ',F5.1, ' Jp = ',f4.1)

      write(*,110),namet,masst,zt,jt
110   format('Target =',A5,' MASS = ',F7.4,' Z = ',F5.1, ' Jt = ',f4.1)

      write(*,150)
150   format('********************************************************')

c***print potential
      write(*,160)
160   format('The following POTENTIALS are defined :')
      write(*,170)
170   format(1X,'KP1#',2x,'KP2#',7X,'TYPE',12X,'SHAPE')

180   format(3X,I1,2x,I2,5x,'Coulomb',10X,'CHARGE',10x,'a1=',F6.3,
     &   '   a2=',
     &       F6.3,'   rc=',F6.3)

190   format(3x,I1,2x,I2,5x,'Volume',10x,I2,'=',A15,"v=",F6.3,'    av=',
     &       F6.3,'   rv=',F6.3,' w=',F6.3, ' aw=',F6.3, ' rw=',
     &        F6.3)

200   format(3x,I1,2x,I2,5x,'Projtl S.O.',5x,I2,'=',A15,'vsov=',F6.3,
     &       '  rsov=',F6.3,'  asov=',F6.3,'  vsow=',F6.3,
     &         '  rsow=',F6.3,'  asow=',F6.3)

210   format(3x,I1,2x,I2,5x,'Surface',9x,I2,'=',A15,'vd=',F6.3,
     &       '  rvd=',F6.3,'  avd=',F6.3,'  wd=',F6.3,
     &         '  rwd=',F6.3,'  awd=',F6.3)


      write(*,180) 1,1,a1,a2,rc
      write(*,190)1,1,ptype,potype(ptype),uv,av,rv,
     &           uw,aw,rw
      if (vsov/=0 .or. vsow/=0 ) then
         write(*,200) 1,1,ptype,potype(ptype),vsov,rsov,asov,
     &                vsow,rsow,asow

      end if
      if (vd/=0 .or. wd/=0 ) then
         write(*,210) 1,1,ptype,potype(ptype),vd,rvd,avd,wd,rwd,awd
      end if

      write(*,150)
      end subroutine


c-----------------------------------------------------------------------
      function potype(a)
      integer :: a
      character(len=15) :: potype
      select case(a)
      case(1)
         potype='Woods-Saxon'
      case(2)
         potype='Gaussian'
      end select
      end function

c-----------------------------------------------------------------------
cWrite output files units and names
       subroutine fkind()
       character*40 flkind(9999)
       integer writf(9999),nwrit,i
       flkind(:) = ' '
       flkind(1)='local copy of input'
       written(1)=.TRUE.
       flkind(2)='phase-shifts '
       written(2)=.TRUE.
       flkind(3)='S-matrix '
       written(3)=.TRUE.
       flkind(4)='radial part of Wf '
       written(4)=.TRUE.

       flkind(8)='channels coupling index'
       written(8)=.TRUE.

       flkind(30)='potential used in the calcuation'
       written(30)=.TRUE.


       nwrit = 0
       do i=1,9999
        if(written(i)) then
        flkind(i) = trim(flkind(i))//'.'
        nwrit = nwrit+1
        writf(nwrit) = i
        endif
       enddo
       write(*,990) (writf(i),flkind(writf(i)),i=1,nwrit)
990    format(/'  The following files have been created:',
     X  /(2x,2(i3,':',1x,a40)))
       return
       end subroutine





      end module
