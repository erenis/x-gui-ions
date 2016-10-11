!-----------------------------------------------------------------------
!
      module mdiag1
!
! Fortran90 wrappers to 1d OpenMP PIC library libmdiag1.f
! get_funit returns an unconnected fortran unit number
! dafopennc1 opens new binary file for complex 1d scalar data.
! dafwritec1 writes record in direct access binary file
! mcspect1 performs frequency analysis of complex time series
!          calls CSPECT1
! micspect1 performs incremental frequency analysis of complex time
!           series for one time step
!           calls ICSPECT1
! mvpdist1 calculates 1 component velocity distribution, velocity
!          moments, and entropy with segmented particle array
!          calls VPDIST1
! mvdist1 calculates 1 component velocity distribution, velocity
!         moments, and entropy with standard particle array
!         calls VPDIST1
! settraj1 sets test charge distribution by setting a particle id in
!          particle location 3 or 5
!          calls STPTRAJ1 or STPTRAJ13
! mptraj1 copies tagged particles in ppart to array partt
!         calls PTRAJ1
! setmbeam1 marks beam particles by setting a particle id in particle
!           location 3 or 5
!           calls STPBEAM1 or STPBEAM13
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: october 5, 2016
!
      use libmdiag1_h
      implicit none
!
      contains
!
!-----------------------------------------------------------------------
      function get_funit(start) result(funit)
! this function returns an unconnected fortran unit number,
! starting with unit = start.  returns -1 if none found
      integer, intent(in) :: start
      integer :: funit
! local data
      integer :: i
      logical :: connected
      funit = -1
! check connection status
      do i = start, 99
         inquire(unit=i,opened=connected)
         if (.not.connected) then
            funit = i
            exit
         endif
      enddo
      end function
!
!-----------------------------------------------------------------------
      subroutine dafopennc1(fc,iunit,nrec,fname)
! this subroutine opens new direct access binary file
! fc = data array to be written in each record
! iunit = fortran unit number to be used 
! nrec = returns initial record number
! fname = file name
      implicit none
      integer, intent(inout) :: iunit, nrec
      complex, dimension(:), intent(in) :: fc
      character(len=*), intent(in) :: fname
! local data
      integer :: lrec
      inquire(iolength=lrec) fc(1); lrec = lrec*size(fc)
      iunit = get_funit(iunit)
      open(unit=iunit,file=fname,form='unformatted',access='direct',    &
     &recl=lrec,status='replace')
      nrec = 1
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine dafwritec1(fc,tdiag,iunit,nrec,nx)
! this subroutine writes record in direct access binary file
! fc = data array to be written
! iunit = fortran unit number to be used 
! nrec = record number for write (then updated to next record)
! nx = number of elements to be written in record
      implicit none
      integer, intent(in) :: iunit, nx
      integer, intent(inout) :: nrec
      real, intent(inout) :: tdiag
      complex, dimension(:), intent(in) :: fc
! local data
      integer :: j
      integer, dimension(4) :: itime
      double precision :: dtime
      call dtimer(dtime,itime,-1)
      write (unit=iunit,rec=nrec) (fc(j),j=1,nx)
      nrec = nrec + 1
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mcspect1(fc,wm,pkw,t0,dt,tdiag,nt,iw,modesx)
      integer, intent(in) :: nt, iw, modesx
      real, intent(in) :: t0, dt
      real, intent(inout) :: tdiag
      complex, dimension(:,:), intent(in) :: fc
      real, dimension(:), intent(in) :: wm
      real, dimension(:,:,:), intent(inout) :: pkw
! local data
      integer :: ntd, iwd, modesxd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ntd = size(fc,1); modesxd = size(fc,2)
      iwd = size(wm,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call CSPECT1(fc,wm,pkw,t0,dt,nt,iw,modesx,ntd,iwd,modesxd)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine micspect1(fc,wm,pkw,pks,time,t0,tdiag,nt,iw,modesx,nx, &
     &norm)
      integer, intent(in) :: nt, iw, modesx, nx, norm
      real, intent(in) :: time, t0
      real, intent(inout) :: tdiag
      complex, dimension(:), intent(in) :: fc
      real, dimension(:), intent(in) :: wm
      real, dimension(:,:,:), intent(inout) :: pkw
      double precision, dimension(:,:,:), intent(inout) :: pks
! local data
      integer :: iwd, modesxd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      modesxd = size(fc,1); iwd = size(wm,1)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call ICSPECT1(fc,wm,pkw,pks,time,t0,nt,iw,modesx,nx,norm,iwd,     &
     &modesxd)
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mvpdist1(ppart,kpic,sfv,fvm,tdiag,np,nmv)
! calculates 1d velocity distribution, velocity moments, and entropy
! with segmented particle array
      integer, intent(in) :: np, nmv
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      integer, dimension(:), intent(in) :: kpic
      real, dimension(:,:,:), intent(inout) :: sfv
      real, dimension(:,:), intent(inout) :: fvm
! local data
      integer :: idimp, nppmx, nmvf, idimv, mx1
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      nmvf = size(sfv,1); idimv = size(sfv,2)
      mx1 = size(kpic,1)
! initialize timer
      call dtimer(dtime,itime,-1)
      if (idimv==1) then
         call VPDIST1(ppart,kpic,sfv,fvm,idimp,nppmx,mx1,np,nmv,nmvf)
!     else if (idimv==3) then
!        call VDIST13(part,sfv,fvm,idimp,np,nmv,nmvf)
      endif
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mvdist1(part,fv,fvm,tdiag,np,nmv)
! calculates 1d velocity distribution, velocity moments, and entropy
! with standard particle array
      integer, intent(in) :: np, nmv
      real, intent(inout) :: tdiag
      real, dimension(:,:), intent(in) :: part
      real, dimension(:,:), intent(inout) :: fv, fvm
! local data
      integer :: idimp, nmvf, idimv
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(part,1)
      nmvf = size(fv,1); idimv = size(fv,2)
! initialize timer
      call dtimer(dtime,itime,-1)
      if (idimv==1) then
         call VDIST1(part,fv,fvm,idimp,np,nmv,nmvf)
!     else if (idimv==3) then
!        call VDIST13(part,fv,fvm,idimp,np,nmv,nmvf)
      endif
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine setptraj1(ppart,kpic,iprobt,nst,vtx,vtsx,dvtx,np,nprobt&
     &)
! sets test charge distribution by setting a particle id in particle
! location 3 or 5
! nst = type of test particle distribution
!   1 = uniformly distribution in real space
!   2 = uniform distribution in velocity space
!   3 = velocity slice at vtsx +- dvtx/2
      real, dimension(:,:,:), intent(inout) :: ppart
      integer, dimension(:), intent(in) :: kpic
      integer, dimension(:), intent(inout) :: iprobt
      integer, intent(in) :: nst, np
      integer, intent(inout) :: nprobt
      real, intent(in) :: vtx, vtsx, dvtx
! local data
      integer :: idimp, mx1, nppmx
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mx1 = size(kpic,1)
! call low level procedure
      if (idimp > 5) then
         call STPTRAJ13(ppart,kpic,iprobt,nst,vtx,vtsx,dvtx,idimp,nppmx,&
     &mx1,np,nprobt)
      else if (idimp > 2) then
         call STPTRAJ1(ppart,kpic,iprobt,nst,vtx,vtsx,dvtx,idimp,nppmx, &
     &mx1,np,nprobt)
      else
         write (*,*) 'setptraj1 error: idimp=', idimp
         stop
      endif
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mptraj1(ppart,kpic,partt,tdiag)
! this copies tagged particles in ppart to array partt
      implicit none
      real, intent(inout) :: tdiag
      real, dimension(:,:,:), intent(in) :: ppart
      integer, dimension(:), intent(in) :: kpic
      real, dimension(:,:), intent(inout) :: partt
! local data
      integer :: idimp, nppmx, mx1, nprobt
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      idimp = size(ppart,1); nppmx = size(ppart,2)
      mx1 = size(kpic,1); nprobt = size(partt,2)
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call PTRAJ1(ppart,kpic,partt,idimp,nppmx,mx1,nprobt)
! record time
      call dtimer(dtime,itime,1)
      tdiag = tdiag + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine setmbeam1(part,npx)
! marks beam particles by setting a particle id in particle location
! 3 or 5
      integer, intent(in) :: npx
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: idimp, nop
      idimp = size(part,1); nop = size(part,2)
! call low level procedure
      if (idimp > 5) then
         call STPBEAM13(part,npx,idimp,nop)
      else if (idimp > 2) then
         call STPBEAM1(part,npx,idimp,nop)
      else
         write (*,*) 'setbeamc1 error: idimp=', idimp
         stop
      endif
      end subroutine
!
      end module
