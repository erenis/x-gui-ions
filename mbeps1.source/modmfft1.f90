!-----------------------------------------------------------------------
!
      module mfft1
!
! Fortran90 wrappers to 1d OpenMP PIC library libmfft1.f
! mfft1_init calculates tables needed by 1d FFTs
!            calls WFFT1RINIT
! mfft1r wrapper function for in place scalar 1d real/complex FFT
!        calls FFT1RXX
! mfft1rn wrapper function for in place vector 1d real/complex FFT
!         calls FFT1R2X, FFT1R3X
! mfft1cr wrapper function for scalar 1d complex to real FFT
!         calls FFT1RXX
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: august 30, 2016
!
      use libmfft1_h
      implicit none
!
! t = scratch array for mfft1rn
      complex, dimension(:), allocatable :: t
      integer :: szt = 0
      save
!
      private :: t, szt
!
      contains
!
!-----------------------------------------------------------------------
      subroutine mfft1_init(mixup,sct,indx)
! calculates tables needed by 1d FFTs
      implicit none
      integer, intent(in) :: indx
      integer, dimension(:), intent(inout) :: mixup
      complex, dimension(:), intent(inout) :: sct
! local data
      integer :: nxhd
! extract dimensions
      nxhd = size(mixup,1)
! call low level procedure
      call WFFT1RINIT(mixup,sct,indx,nxhd)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mfft1r(f,isign,mixup,sct,tfft,indx)
! wrapper function for scalar 1d real/complex FFT
      implicit none
      integer, intent(in) :: isign, indx
      real, dimension(:), intent(inout) :: f
      integer, dimension(:), intent(in) :: mixup
      complex, dimension(:), intent(in) :: sct
      real, intent(inout) :: tfft
! local data
      integer :: nxd, nxhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      nxd = size(f,1)
      nxhd = size(mixup,1)
! check if required size of buffer has increased
      if (szt < nxhd) then
         if (szt /= 0) deallocate(t)
! allocate new buffers
         allocate(t(nxhd))
         szt = nxhd
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      call FFT1RXX(f,t,isign,mixup,sct,indx,nxd,nxhd)
! record time
      call dtimer(dtime,itime,1)
      tfft = tfft + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mfft1rn(f,isign,mixup,sct,tfft,indx)
! wrapper function for n component vector 1d real/complex FFT
      implicit none
      integer, intent(in) :: isign, indx
      real, dimension(:,:), intent(inout) :: f
      integer, dimension(:), intent(in) :: mixup
      complex, dimension(:), intent(in) :: sct
      real, intent(inout) :: tfft
! local data
      integer :: ndim, nxd, nxhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      ndim = size(f,1); nxd = size(f,2)
      nxhd = size(mixup,1)
! check if required size of buffer has increased
      if (szt < ndim*nxhd) then
         if (szt /= 0) deallocate(t)
! allocate new buffers
         allocate(t(ndim*nxhd))
         szt = ndim*nxhd
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
! call low level procedure
      select case(ndim)
      case (2)
         call FFT1R2X(f,t,isign,mixup,sct,indx,nxd,nxhd)
      case (3)
         call FFT1R3X(f,t,isign,mixup,sct,indx,nxd,nxhd)
      end select
! record time
      call dtimer(dtime,itime,1)
      tfft = tfft + real(dtime)
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine mfft1cr(fc,f,mixup,sct,tfft,indx)
! wrapper function for scalar 1d complex to real FFT
! input fc is complex, output f is real
      implicit none
      integer, intent(in) :: indx
      complex, dimension(:), intent(in) :: fc
      real, dimension(:), intent(inout) :: f
      integer, dimension(:), intent(in) :: mixup
      complex, dimension(:), intent(in) :: sct
      real, intent(inout) :: tfft
! local data
      integer :: j, n, isign, nxd, nxhd
      integer, dimension(4) :: itime
      double precision :: dtime
! extract dimensions
      isign = 1
      n = size(fc,1)
      nxd = size(f,1)
      nxhd = size(mixup,1)
! check if required size of buffer has increased
      if (szt < nxhd) then
         if (szt /= 0) deallocate(t)
! allocate new buffers
         allocate(t(nxhd))
         szt = nxhd
      endif
! initialize timer
      call dtimer(dtime,itime,-1)
      if (nxd.ge.(2*n)) then
         do j = 1, n
            f(2*j-1) = real(fc(j))
            f(2*j) = aimag(fc(j))
         enddo
      else
         write (*,*) 'mfft1cr f too small: nxd, 2*n = ', nxd, 2*n
         stop
      endif
! call low level procedure
      call FFT1RXX(f,t,isign,mixup,sct,indx,nxd,nxhd)
! record time
      call dtimer(dtime,itime,1)
      tfft = tfft + real(dtime)
      end subroutine
!
      end module
