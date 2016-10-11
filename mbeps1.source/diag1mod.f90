!-----------------------------------------------------------------------
!
      module diag1d
!
! Fortran90 interface to 1d PIC Fortran77 library diag1lib.f
! diag1mod.f contains diagnostic procedures:
!            defines module diag1d
! wtimer performs wall clock timing.
! get_funit returns an unconnected fortran unit number.
! bfopen => bfopen1 opens binary file for real 1d scalar data.
! bfopen => bfvopen1 opens binary file for real 1d vector data
! bfopen => bfcopen1 opens binary file for complex 1d scalar data.
! bfopen => bfvcopen1 opens binary file for complex 1d vector data
! writebf => ifwrite1 writes real 1d scalar data to a binary file.
!            calls FWRITE1
! writebf => ifvwrite1 writes real 1d vector data to a binary file.
!            call FWRITE1
! writebf => ifcwrite1 writes complex 1d scalar data to a binary file.
!            calls FCWRITE1
! writebf => ifvcwrite1 writes complex 1d vector data to a binary file.
!            call FCWRITE1
! readbf => ifread1 reads real 1d scalar data from a binary file.
!           calls FREAD1
! readbf => ifvread1 reads real 1d vector data from a binary file.
!           calls FREAD1
! readbf => ifcread1 reads complex 1d scalar data from a binary file.
!           calls FCREAD1
! readbf => ifvcread1 reads complex 1d vector data from a binary file.
!           calls FCREAD1
! vdist => ivdist1 calculates 1 or 3 component velocity distribution,
!          velocity moments, and entropy.
!          calls VDIST13 or VDIST1
! psdist =>  ipsdist1 calculates 1 or 3 component phase space
!            distribution, velocity moments, and entropy
!            calls PSDIST13 or PSDIST1
! settraj1 sets test charge distribution
! open_graphs open graphics device
!             calls GROPEN and SETNPLT
! close_graphs close graphics device
!              calls GRCLOSE
! reset_graphs reset graphics device
!              calls RSTSCRN
! displays => idscaler1 displays 1d scalar field in real space.
!             calls DISPS
! displayv => idvector1 displays amplitude of a 1d vector field in real
!             space.
!             calls DISPS
! displayfv => idisplayfv1 displays velocity distribution functions.
!              calls DISPR
! displayt => idisplayt1 displays time history of trajectories
! displayw => idisplayw1 displays time history of electric field,
!             kinetic and total energies.
!             calls DISPR
! grasp => igrasp13 displays particle in phase space.
!          calls GRASP13
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: july 4, 2013
!
      use globals, only: LINEAR, QUADRATIC, CUBIC
      use init1d, only: idrun, indx, ntp, nta, nte, psolve, tend, dt,   &
     &ndim, omx, omy, omz, ci, ntd, ntj, t0, ceng, pot1d, modesxp,      &
     &nprec, fpname, vpot1d, modesxa, narec, faname, em1d, modesxe,     &
     &nerec, fename, den1d, modesxd, ndrec, fdname, vcur1d, modesxj,    &
     &njrec, fjname
!
      implicit none
      private
      public :: LINEAR, QUADRATIC, CUBIC
      public :: idrun, indx, ntp, nta,  nte, psolve
      public :: tend, dt, ndim, omx, omy, omz, ci, ntd, ntj
      public :: t0, ceng
      public :: pot1d, modesxp, nprec, fpname
      public :: vpot1d, modesxa,  narec, faname
      public :: em1d, modesxe,  nerec, fename
      public :: den1d, modesxd, ndrec, fdname
      public :: vcur1d, modesxj, njrec, fjname
      public :: wtimer, get_funit
      public :: bfopen, writebf, readbf
!     public :: writef
      public :: vdist, psdist, settraj1
      public :: open_graphs, close_graphs, reset_graphs
      public :: displays, displayv, displayfv, displayt, displayw, grasp
!
! npl = (0,1) = display is (off,on)
      integer, save :: npl = 1
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine FCWRITE1(f,nx,nxv,iunit,nrec,lrec,name)
         implicit none
         integer :: nx, nxv, iunit, nrec, lrec
!        complex, dimension(*) :: f
         complex :: f
         character(len=*) :: name
         end subroutine
      end interface
      interface
         subroutine FCREAD1(f,nx,nxv,iunit,nrec,lrec,name,ierr)
         implicit none
         integer :: nx, nxv, iunit, nrec, lrec, ierr
!        complex, dimension(*) :: f
         complex :: f
         character(len=*) :: name
         end subroutine
      end interface
      interface
         subroutine FWRITE0(f,nxp,iunit,nrec,name)
         implicit none
         integer :: nxp, iunit, nrec
         real, dimension(nxp) :: f
         character(len=*) :: name
         end subroutine
      end interface
      interface
         subroutine VDIST1(part,fv,fvm,idimp,np,nmv,nmvf)
         integer :: idimp, np, nmv, nmvf
         real, dimension(idimp,np) :: part
         real, dimension(nmvf) :: fv
         real, dimension(3) :: fvm
         end subroutine
      end interface
      interface
         subroutine VDIST13(part,fv,fvm,idimp,np,nmv,nmvf)
         integer :: idimp, np, nmv, nmvf
         real, dimension(idimp,np) :: part
         real, dimension(nmvf,3) :: fv
         real, dimension(3,3) :: fvm
         end subroutine
      end interface
      interface
         subroutine PSDIST1(part,fps,fpsm,psm,nx,nxb,idimp,np,nmv,nmvf)
         integer :: nx, nxb, idimp, np, nmv, nmvf
         real, dimension(idimp,np) :: part
         real, dimension(nmvf,nxb) :: fps
         real, dimension(3,nxb) :: fpsm
         double precision, dimension(2,nxb) :: psm
         end subroutine
      end interface
      interface
         subroutine PSDIST13(part,fps,fpsm,psm,nx,nxb,idimp,np,nmv,nmvf)
         integer :: nx, nxb, idimp, np, nmv, nmvf
         real, dimension(idimp,np) :: part
         real, dimension(nmvf,3,nxb) :: fps
         real, dimension(3,3,nxb) :: fpsm
         double precision, dimension(2,3,nxb) :: psm
         end subroutine
      end interface
      interface
         subroutine GROPEN
         implicit none
         end subroutine
      end interface
      interface
         subroutine SETNPLT(nplt,irc)
         implicit none
         integer :: nplt, irc
         end subroutine
      end interface
      interface
         subroutine DISPR(f,label,xmin,xmax,isc,ist,mks,nx,nxv,ngs,chr,c&
     &hrs,irc)
         implicit none
         integer :: isc, ist, mks, nx, nxv, ngs, irc
         real :: xmin, xmax
         character(len=*) :: label, chr
         character(len=*), dimension(ngs) :: chrs
!        real, dimension(*) :: f
         real :: f
         end subroutine
      end interface
      interface
         subroutine GRASP13(part,label,itime,isc,nx,iyp,ixp,idimp,npx,np&
     &,irc)
         implicit none
         integer :: itime, isc, nx, iyp, ixp, idimp, npx, np, irc
         character(len=*) :: label
         real, dimension(idimp,np) :: part
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface bfopen
         module procedure bfopen1
         module procedure bfvopen1
         module procedure bfcopen1
         module procedure bfvcopen1
      end interface
!
      interface writebf
         module procedure ifwrite1
         module procedure ifvwrite1
         module procedure ifcwrite1
         module procedure ifvcwrite1
      end interface
!
      interface readbf
         module procedure ifread1
         module procedure ifvread1
         module procedure ifcread1
         module procedure ifvcread1
      end interface
!
      interface writef
         module procedure ifwrite0
      end interface
!
      interface vdist
         module procedure ivdist1
      end interface
!
      interface psdist
         module procedure ipsdist1
      end interface
!
      interface displays
         module procedure idscaler1
      end interface
!
      interface displayv
         module procedure idvector1
      end interface
!
      interface displayfv
         module procedure idisplayfv1
      end interface
!
      interface displayt
         module procedure idisplayt1
      end interface
!
      interface displayw
         module procedure idisplayw1
      end interface
!
      interface grasp
         module procedure igrasp13
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine wtimer(time,itime,icntrl)
! this subroutine performs wall clock timing
! input: itime, icntrl, output: itime, time
! itime = initial time, gets updated on each call
! icntrl = (-1,0,1) = (initialize,ignore,read) clock
!          clock should be initialized before it is read!
! time = elapsed time, in seconds,  since itime was set
         implicit none
         real, intent(out) :: time
         integer, intent(inout) :: itime
         integer, intent(in), optional :: icntrl
         integer :: ltime, COUNT_RATE, COUNT_MAX
         ltime = 1
         if (present(icntrl)) ltime = icntrl
         if (ltime.eq.0) return
! read clock and write time difference from last clock initialization
         if (ltime.eq.1) then
            ltime = itime
            call system_clock(itime,COUNT_RATE,COUNT_MAX)
            ltime = itime - ltime
            if (ltime.lt.0) ltime = ltime + COUNT_MAX
            time = dble(ltime)/dble(COUNT_RATE)
! initialize clock
         else
           call system_clock(itime,COUNT_RATE,COUNT_MAX)
         endif
         end subroutine wtimer
!
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
         end function get_funit
!
         subroutine bfopen1(f,nx,iunit,nrec,fname)
! this subroutine opens direct access binary file
! for real 1d scalar data
! nrec = (0,-1) open (old, new file), reset to 1 if successful
         implicit none
         integer :: nx, iunit, nrec
         real, dimension(:), pointer :: f
         character(len=*) :: fname
! local data
         integer :: lrec, ierr
         if (nrec > 0) return
         inquire(iolength=lrec) f(1)
         lrec = nx*lrec
         if (nrec==0) then
            open(unit=iunit,file=fname,form='unformatted',access='direct&
     &',recl=lrec,status='old',iostat=ierr)
         else
            open(unit=iunit,file=fname,form='unformatted',access='direct&
     &',recl=lrec,status='replace',iostat=ierr)
         endif
         if (ierr==0) nrec = 1
         end subroutine bfopen1
!
         subroutine ifwrite1(f,nx,iunit,nrec,name,order)
! writes a subset of real 1d scalar data to a direct access binary
! file
         implicit none
         integer :: nx, iunit, nrec
         real, dimension(:), pointer :: f
         character(len=*), optional :: name
         integer, optional :: order
! local data
         integer :: lrec, nxv, inorder
         character(len=1) :: noname = ' '
         nxv = size(f,1)
         if (nrec <= 0) then
            inquire(iolength=lrec) f(1)
            lrec = nx*lrec
            if (nrec < 0) iunit = get_funit(iunit)
         endif
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (present(name)) then
            if (inorder==LINEAR) then
               call FWRITE1(f(1),nx,nxv,iunit,nrec,lrec,name)
            else if (inorder==CUBIC) then
               call FWRITE1(f(3),nx,nxv,iunit,nrec,lrec,name)
            else
               call FWRITE1(f(2),nx,nxv,iunit,nrec,lrec,name)
            endif
         else
            if (inorder==LINEAR) then
               call FWRITE1(f(1),nx,nxv,iunit,nrec,lrec,noname)
            else if (inorder==CUBIC) then
               call FWRITE1(f(3),nx,nxv,iunit,nrec,lrec,noname)
            else
               call FWRITE1(f(2),nx,nxv,iunit,nrec,lrec,noname)
            endif
         endif
         end subroutine ifwrite1
!
         subroutine ifread1(f,nx,iunit,nrec,ierr,name,order)
! reads a subset of real 1d scalar data from a direct access binary
! file
         implicit none
         integer :: nx, iunit, nrec, ierr
         real, dimension(:), pointer :: f
         character(len=*), optional :: name
         integer, optional :: order
! local data
         integer :: lrec, nxv, inorder
         character(len=1) :: noname = ' '
         nxv = size(f,1)
         if (nrec <= 0) then
            inquire(iolength=lrec) f(1)
            lrec = nx*lrec
         endif
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (present(name)) then
            if (inorder==LINEAR) then
               call FREAD1(f(1),nx,nxv,iunit,nrec,lrec,name,ierr)
            else if (inorder==CUBIC) then
               call FREAD1(f(3),nx,nxv,iunit,nrec,lrec,name,ierr)
            else
               call FREAD1(f(2),nx,nxv,iunit,nrec,lrec,name,ierr)
            endif
         else
            if (inorder==LINEAR) then
               call FREAD1(f(1),nx,nxv,iunit,nrec,lrec,noname,ierr)
            else if (inorder==CUBIC) then
               call FREAD1(f(3),nx,nxv,iunit,nrec,lrec,noname,ierr)
            else
               call FREAD1(f(2),nx,nxv,iunit,nrec,lrec,noname,ierr)
            endif
         endif
         end subroutine ifread1
!
         subroutine bfvopen1(f,nx,iunit,nrec,fname)
! this subroutine opens direct access binary file
! for complex 1d vector data
! nrec = (0,-1) open (old, new file), reset to 1 if successful
         implicit none
         integer :: nx, iunit, nrec
         real, dimension(:,:), pointer :: f
         character(len=*) :: fname
! local data
         integer :: lrec, nnx, ierr
         if (nrec > 0) return
         nnx = size(f,1)*nx
         inquire(iolength=lrec) f(1,1)
         lrec = nnx*lrec
         if (nrec==0) then
            open(unit=iunit,file=fname,form='unformatted',access='direct&
     &',recl=lrec,status='old',iostat=ierr)
         else
            open(unit=iunit,file=fname,form='unformatted',access='direct&
     &',recl=lrec,status='replace',iostat=ierr)
         endif
         if (ierr==0) nrec = 1
         end subroutine bfvopen1
!
         subroutine ifvwrite1(f,nx,iunit,nrec,name,order)
! writes a subset of real 1d vector data to a direct access binary
! file
         implicit none
         integer :: nx, iunit, nrec
         real, dimension(:,:), pointer :: f
         character(len=*), optional :: name
         integer, optional :: order
! local data
         integer :: nnx, nxv, lrec, inorder
         character(len=1) :: noname = ' '
         nnx = size(f,1)*nx
         nxv = size(f,1)*size(f,2)
         if (nrec <= 0) then
            inquire(iolength=lrec) f(1,1)
            lrec = nnx*lrec
            if (nrec < 0) iunit = get_funit(iunit)
         endif
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (present(name)) then
            if (inorder==LINEAR) then
               call FWRITE1(f(1,1),nnx,nxv,iunit,nrec,lrec,name)
            else if (inorder==CUBIC) then
               call FWRITE1(f(1,3),nnx,nxv,iunit,nrec,lrec,name)
            else
               call FWRITE1(f(1,2),nnx,nxv,iunit,nrec,lrec,name)
            endif
         else
            if (inorder==LINEAR) then
               call FWRITE1(f(1,1),nnx,nxv,iunit,nrec,lrec,noname)
            else if (inorder==CUBIC) then
               call FWRITE1(f(1,3),nnx,nxv,iunit,nrec,lrec,noname)
            else
               call FWRITE1(f(1,2),nnx,nxv,iunit,nrec,lrec,noname)
            endif
         endif
         end subroutine ifvwrite1
!
         subroutine ifvread1(f,nx,iunit,nrec,ierr,name,order)
! reads a subset of real 1d vector data from a direct access binary
! file
         implicit none
         integer :: nx, iunit, nrec, ierr
         real, dimension(:,:), pointer :: f
         character(len=*), optional :: name
         integer, optional :: order
! local data
         integer :: nnx, nxv, lrec, inorder
         character(len=1) :: noname = ' '
         nnx = size(f,1)*nx
         nxv = size(f,1)*size(f,2)
         if (nrec <= 0) then
            inquire(iolength=lrec) f(1,1)
            lrec = nnx*lrec
         endif
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (present(name)) then
            if (inorder==LINEAR) then
               call FREAD1(f(1,1),nnx,nxv,iunit,nrec,lrec,name,ierr)
            else if (inorder==CUBIC) then
               call FREAD1(f(1,3),nnx,nxv,iunit,nrec,lrec,name,ierr)
            else
               call FREAD1(f(1,2),nnx,nxv,iunit,nrec,lrec,name,ierr)
            endif
         else
            if (inorder==LINEAR) then
               call FREAD1(f(1,1),nnx,nxv,iunit,nrec,lrec,noname,ierr)
            else if (inorder==CUBIC) then
               call FREAD1(f(1,3),nnx,nxv,iunit,nrec,lrec,noname,ierr)
            else
               call FREAD1(f(1,2),nnx,nxv,iunit,nrec,lrec,noname,ierr)
            endif
         endif
         end subroutine ifvread1
!
         subroutine bfcopen1(f,nx,iunit,nrec,fname)
! this subroutine opens direct access binary file
! for complex 1d scalar data
! nrec = (0,-1) open (old, new file), reset to 1 if successful
         implicit none
         integer :: nx, iunit, nrec
         complex, dimension(:), pointer :: f
         character(len=*) :: fname
! local data
         integer :: lrec, ierr
         if (nrec > 0) return
         inquire(iolength=lrec) f(1)
         lrec = nx*lrec
         if (nrec==0) then
            open(unit=iunit,file=fname,form='unformatted',access='direct&
     &',recl=lrec,status='old',iostat=ierr)
         else
            open(unit=iunit,file=fname,form='unformatted',access='direct&
     &',recl=lrec,status='replace',iostat=ierr)
         endif
         if (ierr==0) nrec = 1
         end subroutine bfcopen1
!
         subroutine ifcwrite1(f,nx,iunit,nrec,name,order)
! writes a subset of complex 1d scalar data to a direct access binary
! file
         implicit none
         integer :: nx, iunit, nrec
         complex, dimension(:), pointer :: f
         character(len=*), optional :: name
         integer, optional :: order
! local data
         integer :: lrec, nxv, inorder
         character(len=1) :: noname = ' '
         nxv = size(f,1)
         if (nrec <= 0) then
            inquire(iolength=lrec) f(1)
            lrec = nx*lrec
            if (nrec < 0) iunit = get_funit(iunit)
         endif
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (present(name)) then
            if (inorder==LINEAR) then
               call FCWRITE1(f(1),nx,nxv,iunit,nrec,lrec,name)
            else if (inorder==CUBIC) then
               call FCWRITE1(f(3),nx,nxv,iunit,nrec,lrec,name)
            else
               call FCWRITE1(f(2),nx,nxv,iunit,nrec,lrec,name)
            endif
         else
            if (inorder==LINEAR) then
               call FCWRITE1(f(1),nx,nxv,iunit,nrec,lrec,noname)
            else if (inorder==CUBIC) then
               call FCWRITE1(f(3),nx,nxv,iunit,nrec,lrec,noname)
            else
               call FCWRITE1(f(2),nx,nxv,iunit,nrec,lrec,noname)
            endif
         endif
         end subroutine ifcwrite1
!
         subroutine ifcread1(f,nx,iunit,nrec,ierr,name,order)
! reads a subset of complex 1d scalar data from a direct access binary
! file
         implicit none
         integer :: nx, iunit, nrec, ierr
         complex, dimension(:), pointer :: f
         character(len=*), optional :: name
         integer, optional :: order
! local data
         integer :: lrec, nxv, inorder
         character(len=1) :: noname = ' '
         nxv = size(f,1)
         if (nrec <= 0) then
            inquire(iolength=lrec) f(1)
            lrec = nx*lrec
         endif
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (present(name)) then
            if (inorder==LINEAR) then
               call FCREAD1(f(1),nx,nxv,iunit,nrec,lrec,name,ierr)
            else if (inorder==CUBIC) then
               call FCREAD1(f(3),nx,nxv,iunit,nrec,lrec,name,ierr)
            else
               call FCREAD1(f(2),nx,nxv,iunit,nrec,lrec,name,ierr)
            endif
         else
            if (inorder==LINEAR) then
               call FCREAD1(f(1),nx,nxv,iunit,nrec,lrec,noname,ierr)
            else if (inorder==CUBIC) then
               call FCREAD1(f(3),nx,nxv,iunit,nrec,lrec,noname,ierr)
            else
               call FCREAD1(f(2),nx,nxv,iunit,nrec,lrec,noname,ierr)
            endif
         endif
         end subroutine ifcread1
!
         subroutine bfvcopen1(f,nx,iunit,nrec,fname)
! this subroutine opens direct access binary file
! for complex 1d vector data
! nrec = (0,-1) open (old, new file), reset to 1 if successful
         implicit none
         integer :: nx, iunit, nrec
         complex, dimension(:,:), pointer :: f
         character(len=*) :: fname
! local data
         integer :: lrec, nnx, ierr
         if (nrec > 0) return
         nnx = size(f,1)*nx
         inquire(iolength=lrec) f(1,1)
         lrec = nnx*lrec
         if (nrec==0) then
            open(unit=iunit,file=fname,form='unformatted',access='direct&
     &',recl=lrec,status='old',iostat=ierr)
         else
            open(unit=iunit,file=fname,form='unformatted',access='direct&
     &',recl=lrec,status='replace',iostat=ierr)
         endif
         if (ierr==0) nrec = 1
         end subroutine bfvcopen1
!
         subroutine ifvcwrite1(f,nx,iunit,nrec,name,order)
! writes a subset of complex 1d vector data to a direct access binary
! file
         implicit none
         integer :: nx, iunit, nrec
         complex, dimension(:,:), pointer :: f
         character(len=*), optional :: name
         integer, optional :: order
! local data
         integer :: nnx, nxv, lrec, inorder
         character(len=1) :: noname = ' '
         nnx = size(f,1)*nx
         nxv = size(f,1)*size(f,2)
         if (nrec <= 0) then
            inquire(iolength=lrec) f(1,1)
            lrec = nnx*lrec
            if (nrec < 0) iunit = get_funit(iunit)
         endif
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (present(name)) then
            if (inorder==LINEAR) then
               call FCWRITE1(f(1,1),nnx,nxv,iunit,nrec,lrec,name)
            else if (inorder==CUBIC) then
               call FCWRITE1(f(1,3),nnx,nxv,iunit,nrec,lrec,name)
            else
               call FCWRITE1(f(1,2),nnx,nxv,iunit,nrec,lrec,name)
            endif
         else
            if (inorder==LINEAR) then
               call FCWRITE1(f(1,1),nnx,nxv,iunit,nrec,lrec,noname)
            else if (inorder==CUBIC) then
               call FCWRITE1(f(1,3),nnx,nxv,iunit,nrec,lrec,noname)
            else
               call FCWRITE1(f(1,2),nnx,nxv,iunit,nrec,lrec,noname)
            endif
         endif
         end subroutine ifvcwrite1
!
         subroutine ifvcread1(f,nx,iunit,nrec,ierr,name,order)
! reads a subset of complex 1d vector data from a direct access binary
! file
         implicit none
         integer :: nx, iunit, nrec, ierr
         complex, dimension(:,:), pointer :: f
         character(len=*), optional :: name
         integer, optional :: order
! local data
         integer :: nnx, nxv, lrec, inorder
         character(len=1) :: noname = ' '
         nnx = size(f,1)*nx
         nxv = size(f,1)*size(f,2)
         if (nrec <= 0) then
            inquire(iolength=lrec) f(1,1)
            lrec = nnx*lrec
         endif
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (present(name)) then
            if (inorder==LINEAR) then
               call FCREAD1(f(1,1),nnx,nxv,iunit,nrec,lrec,name,ierr)
            else if (inorder==CUBIC) then
               call FCREAD1(f(1,3),nnx,nxv,iunit,nrec,lrec,name,ierr)
            else
               call FCREAD1(f(1,2),nnx,nxv,iunit,nrec,lrec,name,ierr)
            endif
         else
            if (inorder==LINEAR) then
               call FCREAD1(f(1,1),nnx,nxv,iunit,nrec,lrec,noname,ierr)
            else if (inorder==CUBIC) then
               call FCREAD1(f(1,3),nnx,nxv,iunit,nrec,lrec,noname,ierr)
            else
               call FCREAD1(f(1,2),nnx,nxv,iunit,nrec,lrec,noname,ierr)
            endif
         endif
         end subroutine ifvcread1
!
         subroutine ifwrite0(f,iunit,nrec,name)
         implicit none
         integer :: iunit, nrec
         real, dimension(:) :: f
         character(len=*) :: name
         integer :: nxp
         nxp = size(f)
         call FWRITE0(f,nxp,iunit,nrec,name)
         end subroutine ifwrite0
!
         subroutine ivdist1(part,fv,fvm,np,nmv)
! calculates 1d velocity distribution, velocity moments, and entropy
         integer :: np, nmv
         real, dimension(:,:), pointer :: part
         real, dimension(:,:), pointer :: fv, fvm
! local data
         integer :: idimp, nmvf, idimv
         idimp = size(part,1)
         nmvf = size(fv,1); idimv = size(fv,2)
         if (idimv==1) then
            call VDIST1(part,fv,fvm,idimp,np,nmv,nmvf)
         else if (idimv==3) then
            call VDIST13(part,fv,fvm,idimp,np,nmv,nmvf)
         endif
         end subroutine ivdist1
!
         subroutine ipsdist1(part,fps,fpsm,nx,np,nmv)
! calculates 1d phase space distribution, velocity moments, and entropy
         integer :: nx, np, nmv
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fps, fpsm
! local data
         double precision, dimension(2,size(fps,2),size(fps,3)) :: psm
         integer :: idimp, nmvf, idimv, nxb
         idimp = size(part,1)
         nmvf = size(fps,1); idimv = size(fps,2); nxb = size(fps,3)
         if (idimv==1) then
            call PSDIST1(part,fps,fpsm,psm,nx,nxb,idimp,np,nmv,nmvf)
         else if (idimv==3) then
            call PSDIST13(part,fps,fpsm,psm,nx,nxb,idimp,np,nmv,nmvf)
         endif
         end subroutine ipsdist1
!
!-----------------------------------------------------------------------
         subroutine settraj1(part,iprobt,nst,vtx,vtsx,dvtx,np,nprobt)
! this procedure sets test charge distribution
! input: all, output: iprobt, nprobt
! nst = type of test particle distribution, if ntt > 0
!   1 = uniformly distribution in real space
!   2 = uniform distribution in velocity space
!   3 = velocity slice at vtsx +- dvtx/2
! vtx = thermal velocity of particles in x direction, if nst = 2
! vtsx = center of velocity slice if nst = 3
! dvtx = width of velocity slice if nst = 3
! np = total number of particles in part
! nprobt = number of test charges whose trajectories will be stored.
         integer :: nst, np, nprobt
         real :: vtx, vtsx, dvtx
         real, dimension(:,:) :: part
         integer, dimension(:), pointer :: iprobt
! local data
         integer :: i, j, it, itt, idimp
         real :: st, at
         idimp = size(part,1)
         nprobt = max(0,min(nprobt,np))
         allocate(iprobt(nprobt))
! uniform distribution in real space
         if (nst==1) then
            it = np/nprobt
            do j = 1, nprobt
            iprobt(j) = 1 + it*(j - 1)
            enddo
! uniform distribution in velocity space in x direction
         else if (nst==2) then
            st = 4.0*vtx
            it = nprobt/2
            at = real(it)
            st = at/st
            at = at + 1.0 + 0.5*(nprobt - 2*it)
            do j = 1, nprobt
            iprobt(j) = 0
            enddo
            do j = 1, np
            it = part(2,j)*st + at
            if ((it > 0) .and. (it <= nprobt)) then
               if (iprobt(it)==0) iprobt(it) = j
            endif
            enddo
! remove any holes in probe array
            it = 0
            do j = 1, nprobt
            if (iprobt(j) > 0) then
               it = it + 1
               iprobt(it) = iprobt(j)
            endif
            enddo
            nprobt = it
! velocity slice in x direction
         else if (nst==3) then
            deallocate(iprobt)
            st = 1.0/dvtx
            itt = vtsx*st + 0.5
! determine how many particles have velocities near vtsx
            i = 0
            do j = 1, np
            it = part(2,j)*st + 0.5
            if (it==itt) i = i + 1
            enddo
            nprobt = i
            allocate(iprobt(nprobt))
! find particle locations
            i = 0
            do j = 1, np
            it = part(2,j)*st + 0.5
            if (it==itt) then
               i = i + 1
               iprobt(i) = j
            endif
            enddo
         endif
         end subroutine settraj1
!
         function open_graphs(nplot) result(irc)
! open graphics device
         integer :: nplot, irc
         call GROPEN
         call SETNPLT(nplot,irc)
         end function open_graphs
!
         subroutine close_graphs
! close graphics device
         call GRCLOSE
         end subroutine close_graphs
!
         subroutine reset_graphs
! reset graphics device
         call RSTSCRN
         npl = 1
         end subroutine reset_graphs
!
         subroutine idscaler1(f,label,itime,isc,ist,nx,irc,inorder)
! displays 1d scalar field in real space
! f = 1d scalar field in real space
! label = field label
! itime = current time step
! isc = power of 2 scale of range of values of f
! ist = flag for choosing positive and/or negative values
! nx = system length in x direction
! irc = return code (0 = normal return)
         implicit none
         integer :: itime, isc, ist, nx, irc
         integer, optional :: inorder
         character(len=*) label
         real, dimension(:), pointer :: f
! local data
         integer :: nxv, mx, order
         real :: xmin, xmax
         character(len=12) :: lbl
   91    format(' T = ',i7)
         if (npl==0) return
         nxv = size(f)
         xmin = 0.0
         write (lbl,91) itime
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            mx = min(nx+1,nxv); xmax = real(mx - 1)
            call DISPS(f(1),label,xmin,xmax,isc,ist,mx,lbl,irc)
         else if (inorder==CUBIC) then
            mx = min(nx+1,nxv-2); xmax = real(mx - 1)
            call DISPS(f(3),label,xmin,xmax,isc,ist,mx,lbl,irc)
         else
            mx = min(nx+1,nxv-1); xmax = real(mx - 1)
            call DISPS(f(2),label,xmin,xmax,isc,ist,mx,lbl,irc)
         endif
         if (irc > 127) then
            npl = irc - 128
            if (npl==0) call CLRSCRN
            irc = 0
         endif
         end subroutine idscaler1
!
         subroutine idvector1(f,label,itime,isc,ist,idm,nx,irc,inorder)
! displays 1d vector field in real space
! f = 1d vector field in real space
! label = field label
! itime = current time step
! isc = power of 2 scale of range of values of f
! ist = flag for choosing positive and/or negative values
! idm = (1,2,3) = display (components,amplitude,both)
! nx = system length in x direction
! irc = return code (0 = normal return)
         implicit none
         integer :: itime, isc, ist, idm, nx, irc
         integer, optional :: inorder
         character(len=*) label
         real, dimension(:,:), pointer :: f
! local data
         real, dimension(size(f,2)) :: g
         integer :: i, j, nxv, mx, order
         real :: xmin, xmax, sum1
         character(len=12) :: lbl
         character(len=2) :: c
   91    format(' T = ',i7)
         if (npl==0) return
         write (lbl,91) itime
         order = QUADRATIC
         if (present(inorder)) order = inorder
         nxv = size(f,2)
         xmin = 0.0
! display components
         if (idm /= 2) then
            do i = 1, size(f,1)
            g = f(i,:)
            if (i.eq.1) then
               c = ':Y'
            else if (i.eq.2) then
               c = ':Z'
            else
               write (c,'(":",i1)') i
            endif
            if (order==LINEAR) then
               mx = min(nx+1,nxv); xmax = real(mx - 1)
               call DISPS(g(1),label//c,xmin,xmax,isc,ist,mx,lbl,irc)
            else if (inorder==CUBIC) then
               mx = min(nx+1,nxv-2); xmax = real(mx - 1)
               call DISPS(g(3),label//c,xmin,xmax,isc,ist,mx,lbl,irc)
            else
               mx = min(nx+1,nxv-1); xmax = real(mx - 1)
               call DISPS(g(2),label//c,xmin,xmax,isc,ist,mx,lbl,irc)
            endif
            if (irc > 127) then
               npl = irc - 128
               if (npl==0) call CLRSCRN
               irc = 0
               return
            endif
            if (irc==1) return
            enddo
         endif
! display amplitude
         if (idm /= 1) then
            if (order==LINEAR) then
               mx = nx + 1
            else
               mx = nx + 3
            endif
            mx = min(mx,nxv)
            do j = 1, mx
            sum1 = 0.0
            do i = 1, size(f,1)
            sum1 = sum1 + f(i,j)**2
            enddo
            g(j) = sqrt(sum1)
            enddo
            if (order==LINEAR) then
               mx = min(nx+1,nxv); xmax = real(mx - 1)
               call DISPS(g(1),label,xmin,xmax,isc,ist,mx,lbl,irc)
            else if (inorder==CUBIC) then
               mx = min(nx+1,nxv-2); xmax = real(mx - 1)
               call DISPS(g(3),label,xmin,xmax,isc,ist,mx,lbl,irc)
            else
               mx = min(nx+1,nxv-1); xmax = real(mx - 1)
               call DISPS(g(2),label,xmin,xmax,isc,ist,mx,lbl,irc)
            endif
            if (irc > 127) then
               npl = irc - 128
               if (npl==0) call CLRSCRN
               irc = 0
            endif
         endif
         end subroutine idvector1
!
         subroutine idisplayfv1(fv,fvm,label,itime,nmv,idt,irc)
! displays velocity distribution functions
! fv = velocity distribution
! fvm = velocity moments
! label = long character string label for plot
! itime = current time step
! nmv = number of velocity intervals
! idt = (1,2,3) = display (individual,composite,both) functions
! irc = return code (0 = normal return)
         implicit none
         integer :: itime, nmv, idt, irc
         real, dimension(:,:), pointer :: fv, fvm
         character(len=*) :: label
! isc = 999 = use default display scale
! ist = 1 = display positive values only
! mks = 0 = cycle through line styles
         integer :: isc = 999, ist = 1, mks = 0
         integer :: i, nmvf, nmv2, idimv
         real :: vmax, vmin
         character(len=12) :: c
         character(len=2) :: cs
         character(len=54) :: lbl
         character(len=45) :: chr
         character(len=10), dimension(3) :: chrs
   91    format(', T =',i7)
   92    format(' VD =',f9.6,' VTH =',f9.5)
   93    format(' VTX =',f9.5)
   94    format(' VTX =',f9.5,' VTY =',f9.5,' VTZ =',f9.5)
! chrs = short array of characters to label individual line samples
         data chrs /'    VX    ','    VY    ','    VZ    '/
         if (npl==0) return
         idimv = size(fv,2)
         if ((idimv /= 1) .and. (idimv /= 3)) return
         nmvf = size(fv,1)
         nmv2 = 2*nmv + 1
         write (c,91) itime
! each velocity distributions on its own plot
         if (idt /= 2) then
            do i = 1, idimv
            cs = trim(adjustl(chrs(i)))
            lbl = trim(label)//' VELOCITY DISTR VS '//cs//c
            write (chr,92) fvm(1,i), fvm(2,i)
            vmax = fv(1,i)
            vmin = -vmax
            call DISPR(fv(2,i),lbl,vmin,vmax,isc,ist,mks,nmv2,nmvf,1,chr&
     &,chrs(i),irc)
            if (irc > 127) then
               npl = irc - 128
               if (npl==0) call CLRSCRN
               irc = 0
               return
            endif
            if (irc==1) return
            enddo
         endif
! all velocity distributions on common plot
         if (idt /= 1) then
            lbl = trim(label)//' VELOCITY DISTRS VS '//'V'//c
            if (idimv==1) then
               write (chr,93) fvm(2,1)
               vmax = fv(1,1)
               vmin = -vmax
            else if (idimv==3) then
               write (chr,94) fvm(2,1), fvm(2,2), fvm(2,3)
               vmax = max(fv(1,1),fv(1,2),fv(1,3))
               vmin = -vmax
            endif
            call DISPR(fv(2,1),lbl,vmin,vmax,isc,ist,mks,nmv2,nmvf,idimv&
     &,chr,chrs,irc)
            if (irc > 127) then
               npl = irc - 128
               if (npl==0) call CLRSCRN
               irc = 0
            endif
         endif
         end subroutine idisplayfv1
!
         subroutine idisplayt1(pt,t0,dtt,nt,it,isc,irc)
! displays time history of trajectories
! pt = time history array for trajectories
! t0 = initial trajectory time
! dtt = time between trajectories values
! nt = number of trajectory values to be displayed
! it = (1,2,3,4) plot (x,vx,vy,vz)
! isc = power of 2 scale of range of values of velocity
! irc = return code (0 = normal return)
         integer :: nt, it, isc, irc
         real :: t0, dtt
         real, dimension(:,:,:), pointer :: pt
! ist = 0 = display maximum range
! mks = 0 = cycle through line styles
         integer :: ist = 0, mks = 0
         integer :: i, nttd, ns
         real :: tmin, tmax
         character(len=36) :: lbl
         character(len=4), dimension(4) :: cs 
         character(len=10), dimension(size(pt,3)) :: chrs
! chrs = short array of characters to label individual line samples
         data cs /' X  ',' VX ',' VY ',' VZ '/
         if (npl==0) return
! quit if array is empty or incorrect
         if (nt <= 0) return
         nttd = size(pt,1)*size(pt,2)
         ns = size(pt,3)
! tmin/tmax = range of time values in plot
         tmin = t0
         tmax = t0 + dtt*(nt - 1)
! display individual trajectories
!        lbl = cs(it)//' VERSUS TIME'
         do i = 1, ns
         write (chrs(i),'(i10)') i
         chrs(i) = '    '//adjustl(chrs(i))
!        call DISPR(pt(1,it,1),lbl,tmin,tmax,isc,ist,mks,nt,nttd,1,' ', &
!    &chrs(i),irc)
!        if (irc==1) return
         enddo
! all trajectories on common plot
         lbl = cs(it)//' VERSUS TIME'
         call DISPR(pt(1,it,1),lbl,tmin,tmax,isc,ist,mks,nt,nttd,ns,    &
     &' TEST PARTICLES',chrs,irc)
         end subroutine idisplayt1
!
         subroutine idisplayw1(wt,t0,dtw,nt,irc)
! displays time history of electric field, kinetic, and
! total energies
! wt = time history array for energies
! t0 = initial energy time
! dtw = time between energy values
! nt = number of energy values to be displayed
! irc = return code (0 = normal return)
         integer :: nt, irc
         real :: t0, dtw
         real, dimension(:,:), pointer :: wt
! isc = 999 = use default display scale
! ist = 2 = display minimum range
! mks = 0 = cycle through line styles
         integer :: isc = 999, ist = 2, mks = 0
         integer :: i, ntwd, ns
         real :: tmin, tmax
         character(len=36) :: lbl
         character(len=20), dimension(7) :: cs 
         character(len=10), dimension(7) :: chrs
! chrs = short array of characters to label individual line samples
         data cs /' TOTAL FIELD        ',' ELECTRON KINETIC   ',' ION KI&
     &NETIC     ',' TOTAL              ',' ES FIELD           ',' ET FIE&
     &LD        ',' MAGNETIC FIELD     '/
         data chrs /'TOT FIELD ','ELECT ENRG',' ION ENRG ','TOTAL ENRG',&
     &' EL FIELD ',' ET FIELD ',' B FIELD  '/
         if (npl==0) return
! quit if array is empty or incorrect
         if (nt <= 0) return
         ntwd = size(wt,1)
         ns = min(size(wt,2),7)
! tmin/tmax = range of time values in plot
         tmin = t0
         tmax = t0 + dtw*(nt - 1)
! display individual energies
         do i = 1, ns
         lbl = trim(cs(i))//' ENERGY VERSUS TIME'
         call DISPR(wt(1,i),lbl,tmin,tmax,isc,ist,mks,nt,ntwd,1,' ',chrs&
     &(i),irc)
         if (irc==1) return
         enddo
! all energies on common plot
         lbl = ' ENERGIES VERSUS TIME'
         call DISPR(wt(1,1),lbl,tmin,tmax,isc,ist,mks,nt,ntwd,ns,' ',chr&
     &s,irc)
         end subroutine idisplayw1
!
         subroutine igrasp13(part,np,label,itime,isc,nx,iyp,ixp,npx,irc)
! displays phase space
         implicit none
         integer :: np, itime, isc, nx, iyp, ixp, npx, irc
         character(len=*) :: label
         real, dimension(:,:), pointer :: part
! local data
         integer :: idimp
         idimp = size(part,1)
         call GRASP13(part,label,itime,isc,nx,iyp,ixp,idimp,npx,np,irc)
         end subroutine igrasp13
!
      end module diag1d
