      module channels
      use precision
      use systems

       integer :: lmax !maximum value of l
       integer :: lmin ! minimum value of l 

!      2-body channels
       Type nch_2b
!  |(ls)j>
         integer :: nchmax ! maximum number 
         integer,dimension(:),allocatable :: l
         real*8,dimension(:),allocatable  :: s, j
       End type

       type(nch_2b) :: alpha2b
       
       contains

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine alpha_2b()
c     index of two body projectile and target system
c     in the form of  (l (jp jt) s )   J 
c
c     for a given channel index number,
c     one can directly get the l s J quantum numbers 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer :: l
      real*8 :: s,J
      integer :: ns, nJ
      integer :: nch  ! channel index {(l (jb jx)s) J}

      alpha2b%nchmax=0
      
      do l=lmin,lmax
      
         do ns=nint(2.*abs(jp-jt)),nint(2.*(jp+jt)),2
            s=ns/2.0_dpreal
            
            do nJ=nint(2.*(l+s)),nint(2.*abs(l-s)),2
               J=nJ/2.0_dpreal
               alpha2b%nchmax=alpha2b%nchmax+1
            end do! nJ
         end do ! ns
      end do ! l
      
      write(8,10)alpha2b%nchmax
C     write(*,10)alpha2b%nchmax
10    format('there are',I3,1X,'reaction channels')

      allocate(alpha2b%l(1:alpha2b%nchmax))
      allocate(alpha2b%s(1:alpha2b%nchmax))
      allocate(alpha2b%j(1:alpha2b%nchmax))
     
      nch=1
      write(8,20)
      write(8,30)
C     write(*,20)
C     write(*,30)
      
      do l=lmin,lmax
      
         do ns=nint(2.*abs(jp-jt)),nint(2.*(jp+jt)),2
            s=ns/2.0_dpreal
            
            do nJ=nint(2.*(l+s)),nint(2.*abs(l-s)),2
               J=nJ/2.0_dpreal
               alpha2b%l(nch)=l
               alpha2b%s(nch)=s
               alpha2b%j(nch)=J
               write(8,40)nch,l,jp,jt,s,j
               nch=nch+1
            end do ! nJ
         end do ! ns
      end do !l 
      
      
20    format('---The coupling coefficients are')
30    format(' a2b','|','(',' l ','(',' jp ',' jt ',')',' s ',')',' J ')
40    format(I4,1x,I3,2x,f3.1,2x,f3.1,2x,f3.1,2x,f4.1)      
      
      end subroutine

c-----------------------------------------------------------------------


      end module channels
