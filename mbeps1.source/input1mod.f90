!-----------------------------------------------------------------------
!
      module in1
!
! input1mod.f defines namelists containing input and output variables:
! readnml1 read namelist from unit iuin
! writnml1 write final diagnostic metafile to unit iudm
! written by viktor k. decyk, ucla
! copyright 2011, regents of the university of california
! update: october 6, 2016
!
      implicit none
!
! Basic Input Namelist
      save
! Identification Parameters:
! idrun/idrun0 = run identifier for current/old run
! idcode = code identifier
      integer :: idrun = 1, idrun0 = 0, idcode = 0
!
! Global Simulation Parameters:
! indx = exponent which determines length in x direction, nx=2**indx
      integer :: indx =  11
! psolve = type of poisson solver = (1,2,3)
!     integer :: psolve = PERIODIC_2D
! relativity = (no,yes) = (0,1) = relativity is used
      integer :: relativity = 0
! ci = reciprical of velocity of light
      real :: ci = 0.1
! xtras = fraction of extra particles needed for particle management
      real :: xtras = 0.2
! ndim = number of velocity dimensions = 1 or 3
      integer :: ndim = 3
! nvdist = velocity distribution type
! nvdist = (1,2) = (maxwellian,waterbag) distribution
      integer :: nvdist = 1
! treverse = (0,1) = (no,yes) reverse simulation at end back to start
      integer :: treverse = 0
!
! Background Electron Parameters:
! npx = number of background electrons distributed in x direction
      integer :: npx = 409600
! qme = charge on electron, in units of e
! vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
      real :: qme = -1.0, vtx = 1.0, vty = 1.0, vtz = 1.0
! vx0/vy0/vz0 = drift velocity of electrons in x/y/z direction
      real :: vx0 = 0.0, vy0 = 0.0, vz0 = 0.0
!
! Beam Electron Parameters:
! npxb = number of beam electrons in x direction
      integer :: npxb = 40960
! vtdx/vtdy/vtdz = thermal velocity of beam electrons in x/y/z direction
      real :: vtdx = 1.0, vtdy = 1.0, vtdz = 1.0
! vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
      real :: vdx = 0.0, vdy = 0.0, vdz = 0.0
!
! Time Parameters:
! tend = time at end of simulation, in units of plasma frequency
! dt = time interval between successive calculations
      real :: tend = 45.000, dt = 0.1
!
! Numerical Parameters:
! inorder = interpolation order
! popt = particle optimization scheme
! dopt = charge deposit optimization scheme
! djopt = current deposit optimization scheme
!     integer :: inorder = LINEAR, popt = STANDARD, dopt = LOOKAHEAD
!     integer :: djopt = STANDARD
! ax = half-width of particle in x direction
!     real :: ax = .816497
!     real :: ax = .866025
      real :: ax = .912871
! mx = number of grids in x in sorting tiles
      integer :: mx = 32
! nextrand = (0,N) = generate (default,Nth block) of random numbers
      integer :: nextrand = 0
! ndc = number of corrections in darwin iteration
!     integer :: ndc = 2
!
! Initial Electron Density Parameters:
! ndprof = profile type (uniform=0,linear=1,sinusoidal=2,gaussian=3,
!                        hyperbolic secant squared=4,exponential=5)
      integer :: ndprof = 0
! ampdx = amplitude of density compared to uniform in x
! scaledx = scale length for spatial coordinate in x
! shiftdx = shift of spatial coordinate in x
      real :: ampdx = 0.0, scaledx = 0.0, shiftdx = 0.0
!
! Zero Force Parameter:
! mzf = (0,1) = (no,yes) set forces to zero
      integer :: mzf = 0
!
! External Traveling Wave Driver:
! Ext(x) = e1*sin(k0*x + freq*time) + e2*cos(k0*x - freq*time)
! amodex = wave number which determines k0 = 2*pi*amodex/NX
! freq = frequency of external wave driver
! trmp = ramp-up time for external wave driver
! toff = shut-off time for external wave driver
      real :: amodex = 0.0, freq = 0.0, trmp = 0.0, toff = 0.0
! el0/er0 = external pump amplitude for left-going/right-going wave
!     e1 = el0*(time/trmp), e2 = er0*(time/trmp), if time < trmp
!     e1 = el0,             e2 = er0,             if trmp < time < toff
!     e1 = 0,               e2 = 0,               if time > toff
      real :: el0 = 0.0, er0 = 0.0
!
! Restart Parameters:
! nustrt = (0,1,2) = this is an (old start,new start,restart)
! ntr = number of time steps between restart routine
!     integer :: nustrt = 1, ntr = 0
!
! Energy Diagnostic Parameters:
! ntw = number of time steps between energy diagnostic
      integer :: ntw = 1
!
! Electron Density Diagnostic Parameters:
! ntde = number of time steps between electron density diagnostic
! modesxde = number of modes in x to keep for electron density
!            diagnostic
      integer :: ntde = 0, modesxde = 41
!
! Potential Diagnostic Parameters:
! ntp = number of time steps between potential diagnostic
! ndp = (0,1,2,3) = display (nothing,potential,spectrum,both)
! modesxp = number of modes in x to keep for potential diagnostic
      integer :: ntp = 0, ndp = 1, modesxp = 41
!
! Longitudinal Efield Diagnostic Parameters:
! ntel = number of time steps between longitudinal efield diagnostic
! modesxel = number of modes in x to keep for longitudinal efield
!            diagnostic
      integer :: ntel = 0, modesxel = 41
!
! Power Spectrum Diagnostic Parameters:
! wmin/wmax = minimum/maximum frequency used in power spectrum
! dw = frequency increment used in power spectrum
      real :: wmin = 0.0, wmax = 2.0, dw = 0.01
!
! Velocity-Space Diagnostic Parameter:
! ntv = number of time steps between velocity-space diagnostic
! nmv = number of segments in v for velocity distribution
      integer :: ntv = 0, nmv = 40
!
! Phase-Space Diagnostic
! nts = number of time steps between phase space diagnostic
! nsv = velocity component(s) for phase-space display(s), if nts > 0
! 1 = vx, 2 = vy, 3 = vx and vy, 4 = vz, 5 = vx and vz, 6 = vy and vz,
! 7 = vx and vy and vz
! ntsc = (0,1) = (no,yes) color beam particles
      integer :: nts = 0, nsv = 1, ntsc = 0
!
! Trajectory Diagnostic Parameters:
! ntt = number of time steps between trajectory diagnostic.
! nst = type of test particle distribution, if ntt > 0
! 1 = uniformly distribution in real space
! 2 = uniform distribution in velocity space
! 3 = velocity slice at vtsx +- dvtx/2
! nprobt = number of test charges whose trajectories will be stored.
      integer :: ntt = 0, nst = 0, nprobt = 0
! vtsx = center of velocity slice if nst = 3
! dvtx = width of velocity slice if nst = 3
      real :: vtsx = 0.0, dvtx = 0.1
!
! ntm = number of time steps between momentum diagnostic
!     integer :: ntm = 0
!
! Ion Parameter:
! movion = (0,N) = (no,N) number of moving ions
      integer :: movion = 0
!
! Display Parameter:
! nplot = maximum number of plots per page
      integer :: nplot = 4
!
! Multi-tasking Parameter:
! nvp = number of shared memory nodes (0=default)
      integer :: nvp = 0
!
! define namelist
!     namelist /input1/ idrun, idrun0, idcode, indx, npx, npxb, qme, vtx&
!    &, vty, vtz, vx0, vy0, vz0, vdx, vdy, vdz, vtdx, vtdy, vtdz,       &
!    &relativity, ci, ndim, tend, nextrand, dt, psolve, inorder, popt,  &
!    &dopt, djopt, ax, sortime, ndc, ndprof, ampdx, scaledx, shiftdx,   &
!    &omx, omy, omz, mzf, amodex, freq, trmp, toff, el0, er0, nustrt,   &
!    &ntr, ntw, ntp, nta, ntv, nts, ntm, nte, ntt, nsv, nst, nprobt,    &
!    &vtsx, dvtx, modesxp, modesxa, modesxe, movion, nplot, sntasks
!
! define namelist
      namelist /input1/ idrun, idrun0, idcode, indx, mx, npx, npxb, qme,&
     &vtx, vty, vtz, vx0, vy0, vz0, vdx, vdy, vdz, vtdx, vtdy, vtdz,    &
     &relativity, ci, xtras, ndim, nvdist, tend, dt, ax, nextrand, mzf, &
     &amodex, freq, trmp, toff, el0, er0, ntw, ntde, modesxde, ntp, ndp,&
     &modesxp, ntel, modesxel, wmin, wmax, dw, ntv, nmv, nts, nsv, ntsc,&
     &ntt, nst, nprobt, vtsx, dvtx, movion, nplot, nvp, treverse
!
! Electromagnetic Namelist
! External Magnetic Field Parameters:
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
      real :: omx = 0.0, omy = 0.0, omz = 0.0
!
! Vector Potential Diagnostic Parameters:
! nta = number of time steps between vector potential diagnostic
! nda = (0,1,2,3) = display (nothing,vector potential,spectrum,both)
! modesxa = number of modes in x to keep for vector potential diagnostic
      integer :: nta = 0, nda = 0, modesxa = 41
!
! Radiative Vector Potential Diagnostic Parameters:
! ntrvp = number of time steps between radiative vector potential
!         diagnostic
! ndrvp = (0,1,2,3) = display (nothing,radiative vector potential,
!                              spectrum,both)
! modesxrvp = number of modes in x to keep for radiative vector
!             potential diagnostic
      integer :: ntrvp = 0, ndrvp = 0, modesxrvp = 41
!
! define namelist
      namelist /input1b/ omx, omy, omz, nta, nda, modesxa, ntrvp, ndrvp,&
     &modesxrvp
!
! Ion Namelist
! Background Ion Parameters:
! npxi = number of background ions distributed in x direction
      integer :: npxi =  384
! qmi = charge on ion, in units of e
! rmass = ion/electron mass ratio
      real :: qmi = 1.0, rmass = 100.0
! rtempxi/rtempyi/rtempzi = electron/ion temperature ratio of background
! ions in x/y/z direction
      real :: rtempxi = 1.0, rtempyi = 1.0, rtempzi = 1.0
! vxi0/vyi0/vzi0 = drift velocity of ions in x/y/z direction
      real :: vxi0 = 0.0, vyi0 = 0.0, vzi0 = 0.0
!
! Beam Ion Parameters:
! npxbi = number of beam ions in x direction
      integer :: npxbi =   0
! vdxi/vdyi/vdzi = drift velocity of beam ions in x/y/z direction
      real :: vdxi = 0.0, vdyi = 0.0, vdzi = 0.0
! rtempdxi/rtempdyi/rtempdzi = electron/ion temperature ratio of beam
! ions in x/y/z direction
      real :: rtempdxi = 1.0, rtempdyi = 1.0, rtempdzi = 1.0
!
! Initial Ion Density Parameters:
! ndprofi = ion profile (uniform=0,linear=1,sinusoidal=2,gaussian=3,
!                        hyperbolic secant squared=4,exponential=5)
      integer :: ndprofi = 0
! ampdxi = amplitude of ion density compared to uniform in x
! scaledxi = scale length for spatial ion coordinate in x
! shiftdxi = shift of spatial ion coordinate in x
      real :: ampdxi = 0.0, scaledxi = 0.0, shiftdxi = 0.0
!
! Ion Density Diagnostic Parameters:
! ntdi = number of time steps between ion density diagnostic
! nddi = (0,1,2,3) = display (nothing,ion density,spectrum,both)
! modesxdi = number of modes in x to keep for ion density diagnostic
      integer :: ntdi = 0, nddi = 1, modesxdi = 41
!
! Ion Current Diagnostic Parameters:
! ntji = number of time steps between ion current diagnostic
! ndji = (0,1,2,3) = display (nothing,ion current,spectrum,both)
! modesxji = number of modes in x to keep for ion current diagnostic
      integer :: ntji = 0, ndji = 1, modesxji = 41
!
! Ion Power Spectrum Diagnostic Parameters:
! wimin/wimax = minimum/maximum frequency used in ion power spectrum
! dwi = frequency increment used in ion power spectrum
      real :: wimin = 0.0, wimax = 0.1, dwi = 0.001
!
! define namelist
      namelist /ions1/ npxi, npxbi, qmi, rmass, rtempxi, rtempyi,       &
     &rtempzi, vxi0, vyi0, vzi0, vdxi, vdyi, vdzi, rtempdxi, rtempdyi,  &
     &rtempdzi, ndprofi, ampdxi, scaledxi, shiftdxi, ntdi, nddi,        &
     &modesxdi, ntji, ndji, modesxji, wimin, wimax, dwi
!
! t0 = initial time value
! ceng = energy normalization
      real :: t0 = 0.0, ceng = 1.0
!
! Namelist output for electron density diagnostic
! nderec = current record number for potential writes
      integer :: nderec = 0
! fdename = file name for electron density diagnostic
      character(len=32) :: fdename = 'denek1.0'
! define namelist
      namelist /dene1d/ idrun, indx, ntde, modesxde, nderec, t0, tend,  &
     &dt, ceng, fdename
!
! Namelist output for potential diagnostic
! nprec = current record number for potential writes
      integer :: nprec = 0
! fpname = file name for potential diagnostic
      character(len=32) :: fpname = 'potk1.0'
! define namelist
!     namelist /pot1d/ idrun, indx, ntp, modesxp, psolve, omx, omy, omz,&
!    & nprec, t0, tend, dt, ceng, fpname
      namelist /pot1d/ idrun, indx, ntp, modesxp, nprec, t0, tend, dt,  &
     &ceng, fpname
!
! Namelist output for longitudinal efield diagnostic
! nelrec = current record number for longitudinal efield writes
      integer :: nelrec = 0
! felname = file name for longitudinal efield diagnostic
      character(len=32) :: felname = 'elk1.0'
! define namelist
      namelist /el1d/ idrun, indx, ntel, modesxel, nelrec, t0, tend, dt,&
     &ceng, felname
!
! Namelist output for vector potential diagnostic
! narec = current record number for vector potential writes
!     integer :: narec = 0
! faname = file name for vector potential diagnostic
!     character(len=32) :: faname = 'vpotk1.0'
! define namelist
!     namelist /vpot1d/ idrun, indx, nta, modesxa, psolve, ndim, omx,   &
!    &omy, omz, ci, narec, t0, tend, dt, ceng, faname
!
! Namelist output for electromagnetic diagnostic
! nerec = current record number for electromagnetic writes
!     integer :: nerec = 0
! fename = file name for electromagnetic diagnostic
!     character(len=32) :: fename = 'vpotrk1.0'
! define namelist
!     namelist /em1d/ idrun, indx, nte, modesxe, psolve, ndim, omx, omy,&
!    &omz, ci, nerec, t0, tend, dt, ceng, fename
!
! Namelist output for ion density diagnostic
! ndrec = current record number for ion density writes
      integer :: ndirec = 0
! fdname = file name for ion density diagnostic
      character(len=32) :: fdiname = 'denik1.0'
! define namelist
!     namelist /den1d/ idrun, indx, ntd, modesxd, psolve, ndrec, t0,    &
!    &tend, dt, ceng, fdname
      namelist /deni1d/ idrun, indx, ntdi, modesxdi, ndirec, t0, tend,  &
     &dt, ceng, fdiname 
!
! Namelist output for ion current diagnostic
! njirec = current record number for ion current writes
      integer :: njirec = 0
! fjiname = file name for ion current diagnostic
      character(len=32) :: fjiname = 'vcurk1.0'
! define namelist
!     namelist /vcur1d/ idrun, indx, ntj, modesxj, psolve, ndim, omx,   &
!    &omy, omz, ci, njrec, t0, tend, dt, ceng, fjname
      namelist /vicur1d/ idrun, indx, ntji, modesxji, ndim, omx, omy,   &
     &omz, ci, njirec, t0, tend, dt, ceng, fjiname
!
      contains
!
      subroutine readnml1(iuin)
! read namelist from unit iuin
      implicit none
      integer, intent(in) :: iuin
      open(unit=iuin,file='input1',form='formatted',status='old')
      read (iuin,input1)
      if (movion==1) then
         rewind iuin
         read (iuin,ions1)
      endif
      end subroutine
!
      subroutine writnml1(iudm)
! write final diagnostic metafile to unit iudm
      implicit none
      integer, intent(in) :: iudm
! local data
      character(len=10) :: cdrun
      character(len=32) :: fname
! create string from idrun
      write (cdrun,'(i10)') idrun
      cdrun = adjustl(cdrun)
      fname = 'diag1.'//cdrun
      open(unit=iudm,file=trim(fname),form='formatted',status='replace')
      write (iudm,input1)
      if (ntp > 0) then
         write (iudm,pot1d)
      endif
      if (movion==1) then
         write (iudm,ions1)
         if (ntdi > 0) then
            ceng = 0.0
            write (iudm,deni1d)
         endif
      endif
      end subroutine

      subroutine closenml1(iuin)
! read namelist from unit iuin
      implicit none
      integer, intent(in) :: iuin
      close(iuin)
      end subroutine
!
      end module
