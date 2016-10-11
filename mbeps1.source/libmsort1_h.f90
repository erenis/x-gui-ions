!-----------------------------------------------------------------------
! Interface file for libmsort1.f
      module libmsort1_h
      implicit none
!
      interface
         subroutine PPORDER1L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,nx&
     &,mx,mx1,npbmx,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, nx, mx, mx1, npbmx, ntmax
         integer, intent(inout) :: irc
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(idimp,npbmx,mx1), intent(inout) :: ppbuff
         integer, dimension(mx1), intent(inout) :: kpic
         integer, dimension(2,mx1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1), intent(inout) :: ihole
         end subroutine
      end interface
!
      interface
         subroutine PPORDERF1L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx, &
     &mx1,npbmx,ntmax,irc)
         implicit none
         integer, intent(in) :: idimp, nppmx, mx1, npbmx, ntmax
         integer, intent(inout) :: irc
         real, dimension(idimp,nppmx,mx1), intent(inout) :: ppart
         real, dimension(idimp,npbmx,mx1), intent(inout) :: ppbuff
         integer, dimension(mx1), intent(inout) :: kpic
         integer, dimension(2,mx1), intent(inout) :: ncl
         integer, dimension(2,ntmax+1,mx1), intent(in) :: ihole
         end subroutine
      end interface
!
      end module
