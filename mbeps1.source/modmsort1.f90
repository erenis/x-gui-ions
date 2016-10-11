!-----------------------------------------------------------------------
!
      module msort1
!
! Fortran90 wrappers to 1d OpenMP PIC library libmsort1.f
! morder1 creates list of particles which are leaving tile, buffers
!         outgoing particles, and copies buffers into particle array
!         calls PPORDER1L
! morderf1 buffers outgoing particles and copies buffers into particle
!          array
!          calls PPORDERF1L
! wmporder1 generic procedure to perform particle reordering into tiles
!           calls mporder1 or morderf1
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: august 4, 2016
!
      use libmsort1_h
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      subroutine mporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,mx,irc)
! performs particle reordering into tiles
      implicit none
      integer, intent(in) :: nx, mx
      real, intent(inout) :: tsort
      integer, intent(inout) :: irc
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(inout) :: ppbuff
      integer, dimension(:), intent(inout) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
! local data
      integer :: idimp, nppmx, npbmx, ntmax, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      npbmx = size(ppbuff,2)
      mx1 = size(kpic,1); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPORDER1L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,nx,mx,mx1, &
     &npbmx,ntmax,irc)
! record time
      call dtimer(dtime,itime,1)
      tsort = tsort + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine morderf1(ppart,ppbuff,kpic,ncl,ihole,tsort,irc)
! performs particle reordering into tiles,
! does not create list of particles which are leaving tile
      implicit none
      real, intent(inout) :: tsort
      integer, intent(inout) :: irc
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(inout) :: ppbuff
      integer, dimension(:), intent(inout) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(in) :: ihole
! local data
      integer :: idimp, nppmx, npbmx, ntmax, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      npbmx = size(ppbuff,2)
      mx1 = size(kpic,1); ntmax = size(ihole,2) - 1
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PPORDERF1L(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx,mx1,npbmx,&
     &ntmax,irc)
! record time
      call dtimer(dtime,itime,1)
      tsort = tsort + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine wmporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,mx,list,&
     &irc)
! generic procedure to perform particle reordering into tiles
! list = (true,false) = list of particles leaving tiles found in push
      implicit none
      integer, intent(in) :: nx, mx
      real, intent(inout) :: tsort
      integer, intent(inout) :: irc
      logical, intent(in) :: list
      real, dimension(:,:,:), intent(inout) :: ppart
      real, dimension(:,:,:), intent(inout) :: ppbuff
      integer, dimension(:), intent(inout) :: kpic
      integer, dimension(:,:), intent(inout) :: ncl
      integer, dimension(:,:,:), intent(inout) :: ihole
! do not calculate list of particles leaving tile
      if (list) then
         call morderf1(ppart,ppbuff,kpic,ncl,ihole,tsort,irc)
! calculate list of particles leaving tile
      else
         call mporder1(ppart,ppbuff,kpic,ncl,ihole,tsort,nx,mx,irc)
      endif
      if (irc /= 0) then
         write (*,*) 'wmporder1 error: irc=', irc
         stop
      endif
      end subroutine
!
      end module
