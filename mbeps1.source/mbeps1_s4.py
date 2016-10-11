#-----------------------------------------------------------------------
# 1D Electrostatic OpenMP PIC code
# written by Viktor K. Decyk and Joshua Kelly, UCLA
# copyright 2016, regents of the university of california
"""

"""

import math
import numpy
from types import *
from libmpush1 import *
from fomplib import *
from fgraf1 import *
from dtimer import *

"""
This imports the gui code
"""
import sys
sys.path.append('./gui')
from ProceduralInterface import *

int_type = numpy.int32
double_type = numpy.float64
#if (minit1.fprecision()==0):
float_type = numpy.float32
complex_type = numpy.complex64
#else:
#  float_type = numpy.float64
#  complex_type = numpy.complex128
#  print "using double precision"

#Some boilderplate
def changeVarsCallback(obj, to):
   for key in to.var:
      setattr(obj, key, rightType(to.var[key]) )

def resetCallback(obj, to):
   obj.closenml1(obj.iuin)
   init(obj)

def rightType(val):
   ints, reals, complexs = int_type, float_type, complex_type
   if type(val) is IntType:
      return np.array([val], ints)
   elif type(val) is FloatType:
      return np.array([val],reals)
   elif type(val) is ComplexType:
      return np.array([val], complexs)

#initialize data structures
def init(s1):
    s1.timedirection = 0 #Set time forward.  Not a fortan variable
    s1.phikw = []
    s1.phikw_time = []
    s1.phikw_dta = []
    s1.phikwi = []
    s1.phikwi_time = []
    s1.phikwi_dta = []

    # s1.idimp = number of particle coordinates = 2
    # s1.ipbc = particle boundary condition: 1 = periodic
    s1.idimp = 2; s1.ipbc = 1
    # s1.wke/s1.wki/s1.we = electron/ion kinetic energies and electric field energy
    s1.wke = numpy.zeros((1),float_type)
    s1.wki = numpy.zeros((1),float_type)
    s1.we = numpy.zeros((1),float_type)
    # s1.list = (true,false) = s1.list of particles leaving tiles found in push
    s1.list = True

    # declare scalars for standard code
    s1.npi = 0
    s1.ws = numpy.zeros((1),float_type)

    # declare scalars for OpenMP code
    s1.nppmx = numpy.empty((1),int_type)
    s1.irc = numpy.zeros((1),int_type)

    # declare scalars for diagnostics
    s1.iuin = 8; s1.iudm = 19
    s1.iude = 10; s1.iup = 11; s1.iuel = 12
    s1.iudi = 16

    # declare and initialize timing data
    s1.tinit = 0.0; s1.tloop = 0.0
    s1.itime = numpy.empty((4),numpy.int32)
    s1.ltime = numpy.empty((4),numpy.int32)
    s1.tdpost = numpy.zeros((1),float_type)
    s1.tguard = numpy.zeros((1),float_type)
    s1.tfft = numpy.zeros((1),float_type)
    s1.tfield = numpy.zeros((1),float_type)
    s1.tpush = numpy.zeros((1),float_type)
    s1.tsort = numpy.zeros((1),float_type)
    s1.tdiag = numpy.zeros((1),float_type)
    s1.dtime = numpy.empty((1),double_type)

    # start timing initialization
    dtimer(s1.dtime,s1.itime,-1)
    # read namelist
    s1.readnml1(s1.iuin)
    # override input data
    s1.idcode = 1
    s1.ndim = 1

    # create string from idrun
    s1.cdrun = str(s1.idrun)
    # text output file
    s1.fname = "output1." + s1.cdrun
    s1.iuot = open(s1.fname,"w")
    #s1.nvp = int(input("enter number of nodes: "))
    # initialize for shared memory parallel processing
    omplib.init_omp(s1.nvp)

    # open graphics device
    s1.irc[0] = graf1.open_graphs(s1.nplot)

    if s1.nts > 0:
      pc.addGraph("EDRAWPHASE", "Electron Phase Plot") #Enable ion velocities
      if s1.movion > 0:
        pc.addGraph("IDRAWPHASE", "Ion Phase Plot") #Enable ion velocities

    # initialize scalars for standard code
    # increase number of coordinates for particle tag
    if ((s1.ntt > 0) or ((s1.nts > 0) and (s1.ntsc > 0))):
       s1.idimp += 1
    # s1.np = total number of particles in simulation
    # s1.nx = number of grid points in x direction
    s1.np = s1.npx + s1.npxb;
    s1.nx = int(math.pow(2,s1.indx)); nxh = int(s1.nx/2)
    # s1.npi = total number of ions in simulation
    if (s1.movion > 0):
       s1.npi = s1.npxi + s1.npxbi
    s1.nxe = s1.nx + 2; nxeh = s1.nxe/2
    # s1.mx1 = number of tiles in x direction
    s1.mx1 = int((s1.nx - 1)/s1.mx + 1)
    # s1.nloop = number of time steps in simulation
    # s1.ntime = current time step
    s1.nloop = int(s1.tend/s1.dt + .0001); s1.ntime = 0
    s1.qbme = s1.qme
    s1.affp = float(s1.nx)/float(s1.np)
    if (s1.movion==1):
       s1.qbmi = s1.qmi/s1.rmass
       s1.vtxi = s1.vtx/numpy.sqrt(s1.rmass*s1.rtempxi)
       s1.vtdxi = s1.vtdx/numpy.sqrt(s1.rmass*s1.rtempdxi)

    # check for unimplemented features
    if (s1.list):
       if (s1.ipbc != 1):
          print "s1.ipbc != 1 and s1.list = True not yet supported"
          s1.list = False
          print "s1.list reset to False"

    # allocate data for standard code
    # s1.part = particle array
    s1.part = numpy.empty((s1.idimp,s1.np),float_type,'F')
    # s1.qe = electron charge density with guard cells
    s1.qe = numpy.empty((s1.nxe),float_type,'F')
    # s1.qi = ion charge density with guard cells
    s1.qi = numpy.empty((s1.nxe),float_type,'F')
    # s1.fxe = smoothed electric field with guard cells
    s1.fxe = numpy.empty((s1.nxe),float_type,'F')
    # s1.ffc = form factor array for poisson solver
    s1.ffc = numpy.empty((nxh),complex_type,'F')
    # s1.mixup = bit reverse table for FFT
    s1.mixup = numpy.empty((nxh),int_type,'F')
    # s1.sct = sine/cosine table for FFT
    s1.sct = numpy.empty((nxh),complex_type,'F')
    # s1.kpic = number of electrons in each tile
    s1.kpic = numpy.empty((s1.mx1),int_type,'F')

    # prepare fft tables
    mfft1.mfft1_init(s1.mixup,s1.sct,s1.indx)
    # calculate form factors
    mfield1.mpois1_init(s1.ffc,s1.ax,s1.affp,s1.nx)
    # initialize different ensemble of random numbers
    if (s1.nextrand > 0):
       minit1.mnextran1(s1.nextrand,ndim,s1.np+s1.npi)

    # initialize electons
    # background electrons
    if (s1.npx > 0):
    #  minit1.mudistr1(s1.part,1,s1.npx,s1.nx,s1.ipbc)
       minit1.mfdistr1(s1.part,s1.ampdx,s1.scaledx,s1.shiftdx,1,s1.npx,s1.nx,
                       s1.ipbc,s1.ndprof)
       minit1.wmvdistr1(s1.part,1,s1.vtx,s1.vx0,s1.npx,s1.nvdist)
    # beam electrons
    if (s1.npxb > 0):
       s1.it = s1.npx + 1
    #  minit1.mudistr1(s1.part,s1.it,s1.npxb,s1.nx,s1.ipbc)
       minit1.mfdistr1(s1.part,s1.ampdx,s1.scaledx,s1.shiftdx,s1.it,s1.npxb,
                       s1.nx,s1.ipbc,s1.ndprof)
       minit1.wmvdistr1(s1.part,s1.it,s1.vtdx,s1.vdx,s1.npxb,s1.nvdist)

    # marks electron beam particles
    if ((s1.nts > 0) and (s1.ntsc > 0)):
       mdiag1.setmbeam1(s1.part,s1.npx)

    # find number of electrons in each of mx, tiles: updates s1.kpic, s1.nppmx
    minit1.mdblkp2(s1.part,s1.kpic,s1.nppmx,s1.mx,s1.irc)

    # allocate vector electron data
    s1.nppmx0 = int((1.0 + s1.xtras)*s1.nppmx)
    s1.ntmax = int(s1.xtras*s1.nppmx)
    s1.npbmx = int(s1.xtras*s1.nppmx)
    # s1.ppart = tiled electron array
    s1.ppart = numpy.empty((s1.idimp,s1.nppmx0,s1.mx1),float_type,'F')
    # s1.ppbuff = buffer array for reordering tiled particle array
    s1.ppbuff = numpy.empty((s1.idimp,s1.npbmx,s1.mx1),float_type,'F')
    # s1.ncl = number of particles departing tile in each direction
    s1.ncl = numpy.empty((2,s1.mx1),int_type,'F')
    # s1.ihole = location/destination of each particle departing tile
    s1.ihole = numpy.empty((2,s1.ntmax+1,s1.mx1),int_type,'F')

    # copy ordered electron data for OpenMP: updates s1.ppart and s1.kpic
    mpush1.mpmovin1(s1.part,s1.ppart,s1.kpic,s1.mx,s1.irc)

    # sanity check for electrons
    mpush1.mcheck1(s1.ppart,s1.kpic,s1.nx,s1.mx,s1.irc)

    # initialize background charge density: updates s1.qi
    if (s1.movion==0):
       s1.qi.fill(0.0)
       s1.qmi = -s1.qme
       mpush1.mpost1(s1.ppart,s1.qi,s1.kpic,s1.qmi,s1.tdpost,s1.mx)
       mgard1.maguard1(s1.qi,s1.tguard,s1.nx)

    # initialize ions
    if (s1.movion==1):
       s1.part = numpy.empty((s1.idimp,s1.npi),float_type,'F')
    # s1.kipic = number of ions in each tile
       s1.kipic = numpy.empty((s1.mx1),int_type,'F')
       s1.it = s1.npxi + 1
    # background ions
       if (s1.npxi > 0):
    #     minit1.mudistr1(s1.part,1,s1.npxi,s1.nx,s1.ipbc)
          minit1.mfdistr1(s1.part,s1.ampdxi,s1.scaledxi,s1.shiftdxi,1,
                          s1.npxi,s1.nx,s1.ipbc,s1.ndprofi)
          minit1.wmvdistr1(s1.part,1,s1.vtxi,s1.vxi0,s1.npxi,s1.nvdist)
    # beam ions
       if (s1.npxbi > 0):
    #     minit1.mudistr1(s1.part,s1.it,s1.npxbi,s1.nx,s1.ipbc)
          minit1.mfdistr1(s1.part,s1.ampdxi,s1.scaledxi,s1.shiftdxi,s1.it,
                          s1.npxbi,s1.nx,s1.ipbc,s1.ndprofi)
          minit1.wmvdistr1(s1.part,s1.it,s1.vtdxi,s1.vdxi,s1.npxbi,s1.nvdist)

    # find number of ions in each of mx, tiles: updates s1.kipic, s1.nppmx
       minit1.mdblkp2(s1.part,s1.kipic,s1.nppmx,s1.mx,s1.irc)

    # allocate vector ion data
       s1.nppmx1 = int((1.0 + s1.xtras)*s1.nppmx)
       s1.pparti = numpy.empty((s1.idimp,s1.nppmx1,s1.mx1),float_type,'F')

    # copy ordered ion data for OpenMP: updates s1.pparti and s1.kipic
       mpush1.mpmovin1(s1.part,s1.pparti,s1.kipic,s1.mx,s1.irc)

    # sanity check for ions
       mpush1.mcheck1(s1.pparti,s1.kipic,s1.nx,s1.mx,s1.irc)

    # allocate diagnostic arrays
    # reverse simulation at end back to start
    if (s1.treverse==1):
       s1.nloop = 2*s1.nloop

    # energy time history
    if (s1.ntw > 0):
       s1.mtw = int((s1.nloop - 1)/s1.ntw + 1); s1.itw = 0
    # s1.wt = energy time history array
       s1.wt = numpy.zeros((s1.mtw,4),float_type,'F')
       s1.s = numpy.zeros((4),double_type,'F')
       pc.addGraph("ENERGY", "Energy") #Enable electron velocity

    # allocate scratch arrays for scalar fields
    if ((s1.ntde > 0) or (s1.ntp > 0) or (s1.ntel > 0) or (s1.ntdi > 0)):
       s1.sfieldc = numpy.empty((nxh),complex_type,'F')
       s1.sfield = numpy.empty((s1.nxe),float_type,'F')
    # allocate and initialize scratch array for spectral analysis
       if (s1.ndp > 1):
          s1.iw = int((s1.wmax - s1.wmin)/s1.dw + 1.5)
          s1.wm = numpy.empty((s1.iw),float_type,'F')
          s1.wm[:] = s1.wmin + s1.dw*numpy.linspace(0,s1.iw,s1.iw)
          pc.addGraph("DRAWPHI","Phi(k,w)")
    # s1.cwk = labels for power spectrum display
          s1.cwk = numpy.array([" W > 0"," W < 0"],'S6')
    # allocate and initialize scratch array for ion spectral analysis
       if (s1.movion==1):
          if (s1.nddi > 1):
             s1.iwi = int((s1.wimax - s1.wimin)/s1.dwi + 1.5)
             s1.wmi = numpy.empty((s1.iwi),float_type,'F')
             s1.wmi[:] = s1.wimin + s1.dwi*numpy.linspace(0,s1.iwi,s1.iwi)

    # initialize electron density diagnostic
    if (s1.ntde > 0):
       s1.fdename = "denek1." + s1.cdrun
       s1.modesxde = int(min(s1.modesxde,nxh+1))
    # s1.denet/s1.denit = store selected fourier modes for electron density
       s1.denet = numpy.empty((s1.modesxde),complex_type,'F')
    # open file: updates nderec and possibly s1.iude
       if (s1.nderec==0):
          mdiag1.dafopennc1(s1.denet,s1.iude,s1.nderec,s1.fdename)

    # initialize ion density diagnostic
    if (s1.ntdi > 0):
       pc.addGraph("IONSPECT","Ion Density Spectral Analysis(k,w)")
       s1.fdiname[:] = "denik1." + s1.cdrun
       s1.modesxdi = int(min(s1.modesxdi,nxh+1))
    # s1.denit = store selected fourier modes for ion density
       s1.denit = numpy.empty((s1.modesxdi),complex_type,'F')
    # open file: updates ndirec and possibly s1.iudi
       if (s1.ndirec==0):
          mdiag1.dafopennc1(s1.denit,s1.iudi,s1.ndirec,s1.fdiname)
    # ion spectral analysis
       if ((s1.nddi==2) or (s1.nddi==3)):
          s1.mtdi = int((s1.nloop - 1)/s1.ntdi) + 1; s1.itdi = 0
    # s1.pkwdi = power spectrum for potential
          s1.pkwdi = numpy.empty((s1.modesxdi,s1.iwi,2),float_type,'F')
    # s1.pksdi = accumulated complex spectrum for potential
          s1.pksdi = numpy.zeros((4,s1.modesxdi,s1.iwi),double_type,'F')
    # s1.wkdi = maximum frequency as a function of k for potential
          s1.wkdi = numpy.empty((s1.modesxdi,2),float_type,'F')

    # initialize potential diagnostic
    if (s1.ntp > 0):
       pc.addGraph("DRAWPOT","Potential")
       s1.fpname[:] = "potk1." + s1.cdrun
       s1.modesxp = int(min(s1.modesxp,nxh+1))
    # s1.pott = store selected fourier modes for potential
       s1.pott = numpy.empty((s1.modesxp),complex_type,'F')
    # open file: updates nprec and possibly s1.iup
       if (s1.nprec==0):
          mdiag1.dafopennc1(s1.pott,s1.iup,s1.nprec,s1.fpname)
    # spectral analysis
       if ((s1.ndp==2) or (s1.ndp==3)):
          s1.mtp = int((s1.nloop - 1)/s1.ntp) + 1; s1.itp = 0
    # s1.pkw = power spectrum for potential
          s1.pkw = numpy.empty((s1.modesxp,s1.iw,2),float_type,'F')
    # s1.pks = accumulated complex spectrum for potential
          s1.pks = numpy.zeros((4,s1.modesxp,s1.iw),double_type,'F')
    # s1.wk = maximum frequency as a function of k for potential
          s1.wk = numpy.empty((s1.modesxp,2),float_type,'F')

    # initialize longitudinal efield diagnostic
    if (s1.ntel > 0):
       s1.felname = "elk1." + s1.cdrun
       s1.modesxel = int(min(s1.modesxel,nxh+1))
    # s1.elt = store selected fourier modes for longitudinal efield
       s1.elt = numpy.empty((s1.modesxel),complex_type,'F')
    # open file: updates nelrec and possibly s1.iuel
       if (s1.nelrec==0):
          mdiag1.dafopennc1(s1.elt,s1.iuel,s1.nelrec,s1.felname)

    # initialize velocity diagnostic
    if (s1.ntv > 0):
       pc.addGraph("EVELOCITY", "Electron Velocity") #Enable electron velocity
    # s1.sfv = electron velocity distribution functions in tile
       s1.sfv = numpy.empty((2*s1.nmv+2,1,s1.mx1+1),float_type,'F')
    # s1.fvm = electron vdrift, vth, entropy for global distribution
       s1.fvm = numpy.empty((3,1),float_type,'F')
       s1.mtv = int((s1.nloop - 1)/s1.ntv) + 1; s1.itv = 0
    # s1.fvtm = time history of electron vdrift, vth, and entropy
       s1.fvtm = numpy.zeros((s1.mtv,3,1),float_type,'F')
       s1.sfv[0,:,:] = 2.0*max(4.0*s1.vtx+abs(s1.vx0),
                            4.0*s1.vtdx+abs(s1.vdx))
    # ions
       if (s1.movion==1):
          pc.addGraph("IVELOCITY", "Ion Velocity") #Enable ion velocities
    # s1.sfvi = ion velocity distribution functions in tile
          s1.sfvi = numpy.empty((2*s1.nmv+2,1,s1.mx1+1),float_type,'F')
    # s1.fvmi = ion vdrift, vth, entropy for global distribution
          s1.fvmi = numpy.empty((3,1),float_type,'F')
    # s1.fvtmi = time history of ion vdrift, vth, and entropy
          s1.fvtmi = numpy.zeros((s1.mtv,3,1),float_type,'F')
          s1.sfvi[0,:,:] = 2.0*max(4.0*s1.vtxi+abs(s1.vxi0),
                                4.0*s1.vtdxi+abs(s1.vdxi))

    # initialize trajectory diagnostic
    if (s1.ntt > 0):
    # s1.iprobt = scratch array 
       s1.iprobt = numpy.empty((s1.nprobt),numpy.int32)
       mdiag1.setptraj1(s1.ppart,s1.kpic,s1.iprobt,s1.nst,s1.vtx,s1.vtsx,s1.dvtx,
                        s1.np,s1.nprobt)
       if (s1.nprobt > 16777215):
          print "nprobt overflow = ", s1.nprobt
          exit(1)
    # s1.partt = particle trajectories tracked
       s1.partt = numpy.empty((s1.idimp,s1.nprobt),float_type,'F')
       if ((s1.nst==1) or (s1.nst==2)):
          s1.it = int((s1.nloop - 1)/s1.ntt + 1); itt = 0
    # s1.partd = trajectory time history array
          s1.partd = numpy.empty((s1.it,s1.idimp,s1.nprobt),float_type,'F')
       elif (s1.nst==3):
    # s1.fvtp = velocity distribution function for test particles
          s1.fvtp = numpy.empty((2*s1.nmv+2,1),float_type,'F')
    # s1.fvmtp = vdrift, vth, and entropy for test particles
          s1.fvmtp = numpy.empty((3,1),float_type,'F')
          s1.fvtp[0,:] = 2.0*int(max(4.0*s1.vtx+abs(s1.vx0),
                                  4.0*s1.vtdx+abs(s1.vdx)))

    # initialization time
    dtimer(s1.dtime,s1.itime,1)
    s1.tinit = s1.tinit + float(s1.dtime)
    # start timing loop
    dtimer(s1.dtime,s1.ltime,-1)

    print >> s1.iuot, "program mbeps1\n"

#One time step
def step(s1):
   s1.curtime = s1.ntime*s1.dt
   pc.getEvents(s1)
   pc.setTime(s1.curtime)
   pc.fastForward(s1.curtime, s1)
   # deposit charge with OpenMP: updates s1.qe
   dtimer(s1.dtime,s1.itime,-1)
   s1.qe.fill(0.0)
   dtimer(s1.dtime,s1.itime,1)
   s1.tdpost[0] = s1.tdpost[0] + float(s1.dtime)
   mpush1.mpost1(s1.ppart,s1.qe,s1.kpic,s1.qme,s1.tdpost,s1.mx)
    # add guard cells: updates s1.qe
   mgard1.maguard1(s1.qe,s1.tguard,s1.nx)

    # electron density diagnostic
   if (s1.ntde > 0):
      s1.it = int(s1.ntime/s1.ntde)
      if (s1.ntime==s1.ntde*s1.it):
         s1.sfield[:] = -numpy.copy(s1.qe)
         # transform electron density to fourier space: updates s1.sfield
         isign = -1
         mfft1.mfft1r(s1.sfield,isign,s1.mixup,s1.sct,s1.tfft,s1.indx)
         # calculate smoothed density in fourier space: updates s1.sfieldc
         mfield1.msmooth1(s1.sfield,s1.sfieldc,s1.ffc,s1.tfield,s1.nx)
         # store selected fourier modes: updates s1.denet
         mfield1.mrdmodes1(s1.sfieldc,s1.denet,s1.tfield,s1.nx,s1.modesxde)
         # write diagnostic output: updates nderec
         mdiag1.dafwritec1(s1.denet,s1.tdiag,s1.iude,s1.nderec,s1.modesxde)
         # transform smoothed electron density to real space: updates s1.sfield
         mfft1.mfft1cr(s1.sfieldc,s1.sfield,s1.mixup,s1.sct,s1.tfft,s1.indx)
         mgard1.mdguard1(s1.sfield,s1.tguard,s1.nx)
         # display smoothed electron density
         graf1.dscaler1(s1.sfield,' EDENSITY',s1.ntime,999,0,s1.nx,s1.irc)
         if (s1.irc[0]==1):
            return
         s1.irc[0] = 0

    # deposit ion charge with OpenMP: updates s1.qi
   if (s1.movion==1):
      dtimer(s1.dtime,s1.itime,-1)
      s1.qi.fill(0.0)
      dtimer(s1.dtime,s1.itime,1)
      s1.tdpost[0] = s1.tdpost[0] + float(s1.dtime)
      mpush1.mpost1(s1.pparti,s1.qi,s1.kipic,s1.qmi,s1.tdpost,s1.mx)
      # add guard cells: updates s1.qi
      mgard1.maguard1(s1.qi,s1.tguard,s1.nx)

   # ion density diagnostic
   if (s1.movion==1):
      if (s1.ntdi > 0):
         s1.it = int(s1.ntime/s1.ntdi)
         if (s1.ntime==s1.ntdi*s1.it):
            s1.sfield[:] = numpy.copy(s1.qi)
            # transform ion density to fourier space: updates s1.sfield
            isign = -1
            mfft1.mfft1r(s1.sfield,isign,s1.mixup,s1.sct,s1.tfft,s1.indx)
            # calculate smoothed density in fourier space: updates s1.sfieldc
            mfield1.msmooth1(s1.sfield,s1.sfieldc,s1.ffc,s1.tfield,s1.nx)
            # store selected fourier modes: updates s1.denit
            mfield1.mrdmodes1(s1.sfieldc,s1.denit,s1.tfield,s1.nx,s1.modesxdi)
            # write diagnostic output: updates ndirec
            mdiag1.dafwritec1(s1.denit,s1.tdiag,s1.iudi,s1.ndirec,s1.modesxdi)
            # transform smoothed ion density to real space: updates s1.sfield
            if ((s1.nddi==1) or (s1.nddi==3)):
               mfft1.mfft1cr(s1.sfieldc,s1.sfield,s1.mixup,s1.sct,s1.tfft,s1.indx)
               mgard1.mdguard1(s1.sfield,s1.tguard,s1.nx)
               # display smoothed ion density
               graf1.dscaler1(s1.sfield,' IDENSITY',s1.ntime,999,1,s1.nx,s1.irc)
               if (s1.irc[0]==1):
                  return
               s1.irc[0] = 0
               # ion spectral analysis
               if ((s1.nddi==2) or (s1.nddi==3)):
                  s1.itdi += 1
                  ts = s1.dt*float(s1.ntime)
                  mdiag1.micspect1(s1.denit,s1.wmi,s1.pkwdi,s1.pksdi,ts,s1.t0,
                                   s1.tdiag,s1.mtdi,s1.iwi,s1.modesxdi,s1.nx,-1)
                  # performs frequency analysis of accumulated complex time series
                  s1.wkdi[:,0] = s1.wm[numpy.argmax(s1.pkwdi[:,:,0],axis=1)]
                  s1.wkdi[:,1] = s1.wm[numpy.argmax(s1.pkwdi[:,:,1],axis=1)]
                  # display frequency spectrum
                  s1.phikwi.append( numpy.array(s1.denit, copy=True) )  #add data to phikw
                  s1.phikwi_time.append( s1.ntime*s1.dt )
                  s1.phikwi_dta.append( s1.dt )
                  pc.showPhi(s1.phikwi_time, s1.phikwi, s1.phikwi_dta, s1.dt, plottype="IONSPECT")
                  graf1.dmscaler1(s1.wkdi,'IDENSITY OMEGA VS MODE',s1.ntime,
                                  999,1,s1.modesxdi,s1.cwk,s1.irc)
                  if (s1.irc[0]==1):
                     return
                  s1.irc[0] = 0

    # add electron and ion densities: updates s1.qe
   mfield1.maddqei1(s1.qe,s1.qi,s1.tfield,s1.nx)

    # transform charge to fourier space: updates s1.qe
   isign = -1
   mfft1.mfft1r(s1.qe,isign,s1.mixup,s1.sct,s1.tfft,s1.indx)

    # calculate force/charge in fourier space: updates s1.fxe, s1.we
   mfield1.mpois1(s1.qe,s1.fxe,s1.ffc,s1.we,s1.tfield,s1.nx)

    # transform force to real space: updates s1.fxe
   isign = 1
   mfft1.mfft1r(s1.fxe,isign,s1.mixup,s1.sct,s1.tfft,s1.indx)

    # add external traveling wave field
   ts = s1.dt*float(s1.ntime)
   mfield1.meaddext1(s1.fxe,s1.tfield,s1.amodex,s1.freq,ts,s1.trmp,
                     s1.toff,s1.el0,s1.er0,s1.nx)

    # copy guard cells: updates s1.fxe
   mgard1.mdguard1(s1.fxe,s1.tguard,s1.nx)

    # potential diagnostic
   if (s1.ntp > 0):
      s1.it = int(s1.ntime/s1.ntp)
      if (s1.ntime==s1.ntp*s1.it):
         # calculate potential in fourier space: updates s1.sfieldc
         mfield1.mpot1(s1.qe,s1.sfieldc,s1.ffc,s1.ws,s1.tfield,s1.nx)
         # store selected fourier modes: updates s1.pott
         mfield1.mrdmodes1(s1.sfieldc,s1.pott,s1.tfield,s1.nx,s1.modesxp)
         # write diagnostic output: updates nprec
         mdiag1.dafwritec1(s1.pott,s1.tdiag,s1.iup,s1.nprec,s1.modesxp)
         # transform potential to real space: updates s1.sfield
         if ((s1.ndp==1) or (s1.ndp==3)):
            mfft1.mfft1cr(s1.sfieldc,s1.sfield,s1.mixup,s1.sct,s1.tfft,s1.indx)
            mgard1.mdguard1(s1.sfield,s1.tguard,s1.nx)
            # display potential
            pc.showPotential(s1.sfield)
            graf1.dscaler1(s1.sfield,' POTENTIAL',s1.ntime,999,0,s1.nx,s1.irc)
            if (s1.irc[0]==1):
               return
            s1.irc[0] = 0
         # spectral analysis
         if ((s1.ndp==2) or (s1.ndp==3)):
            s1.itp += 1
            ts = s1.dt*float(s1.ntime)
            mdiag1.micspect1(s1.pott,s1.wm,s1.pkw,s1.pks,ts,s1.t0,s1.tdiag,s1.mtp,s1.iw,
                             s1.modesxp,s1.nx,1)
            # performs frequency analysis of accumulated complex time series
            s1.wk[:,0] = s1.wm[numpy.argmax(s1.pkw[:,:,0],axis=1)]
            s1.wk[:,1] = s1.wm[numpy.argmax(s1.pkw[:,:,1],axis=1)]
            # display frequency spectrum
            s1.phikw.append( numpy.array(s1.pott, copy=True) )  #add data to phikw
            s1.phikw_time.append( s1.ntime*s1.dt )
            s1.phikw_dta.append( s1.dt )
            pc.showPhi(s1.phikw_time, s1.phikw, s1.phikw_dta, s1.dt)
            graf1.dmscaler1(s1.wk,'POTENTIAL OMEGA VS MODE',s1.ntime,999,1,
                            s1.modesxp,s1.cwk,s1.irc)
            if (s1.irc[0]==1):
               return
            s1.irc[0] = 0

    # longitudinal efield diagnostic
   if (s1.ntel > 0):
      s1.it = int(s1.ntime/s1.ntel)
      if (s1.ntime==s1.ntel*s1.it):
         # calculate longitudinal efield in fourier space: updates s1.sfieldc
         mfield1.melfield1(s1.qe,s1.sfieldc,s1.ffc,s1.ws,s1.tfield,s1.nx)
         # store selected fourier modes: updates s1.elt
         mfield1.mrdmodes1(s1.sfieldc,s1.elt,s1.tfield,s1.nx,s1.modesxel)
         # write diagnostic output: updates nelrec
         mdiag1.dafwritec1(s1.elt,s1.tdiag,s1.iuel,s1.nelrec,s1.modesxel)
         # transform longitudinal efield to real space: updates s1.sfield
         if ((s1.ndel==1) or (s1.ndel==3)):
            mfft1.mfft1cr(s1.sfieldc,s1.sfield,s1.mixup,s1.sct,s1.tfft,s1.indx)
            mgard1.mdguard1(s1.sfield,s1.tguard,s1.nx)
            # display longitudinal efield
            graf1.dscaler1(s1.sfield,' ELFIELD',s1.ntime,999,0,s1.nx,s1.irc)
            if (s1.irc[0]==1):
               return
            s1.irc[0] = 0

    # velocity diagnostic
   if (s1.ntv > 0):
      s1.it = int(s1.ntime/s1.ntv)
      if (s1.ntime==s1.ntv*s1.it):
         # calculate electron distribution function and moments
         mdiag1.mvpdist1(s1.ppart,s1.kpic,s1.sfv,s1.fvm,s1.tdiag,s1.np,s1.nmv)
         # store time history electron vdrift, vth, and entropy
         s1.fvtm[s1.itv,:,:] = s1.fvm
         # display electron velocity distributions
         pc.showVelocity(s1.sfv[:,:,s1.mx1], fvm=s1.fvm, plottype="EVELOCITY")
         graf1.displayfv1(s1.sfv[:,:,s1.mx1],s1.fvm,' ELECTRON',s1.ntime,s1.nmv,1,
                          s1.irc)
         if (s1.irc[0]==1):
            return
         s1.irc[0] = 0
         # ion distribution function
         if (s1.movion==1):
            mdiag1.mvpdist1(s1.pparti,s1.kipic,s1.sfvi,s1.fvmi,s1.tdiag,s1.npi,s1.nmv)
            # store time history ion vdrift, vth, and entropy
            s1.fvtmi[s1.itv,:,:] = s1.fvmi
            # display ion velocity distributions
            pc.showVelocity(s1.sfvi[:,:,s1.mx1], fvm=s1.fvmi, plottype="IVELOCITY")
            graf1.displayfv1(s1.sfvi[:,:,s1.mx1],s1.fvmi,' ION',s1.ntime,s1.nmv,1,
                             s1.irc)
            if (s1.irc[0]==1):
               return
            s1.irc[0] = 0
         s1.itv += 1

    # trajectory diagnostic
   if (s1.ntt > 0):
      s1.it = int(s1.ntime/s1.ntt)
      if (s1.ntime==s1.ntt*s1.it):
        # copies trajectories to array s1.partt
         mdiag1.mptraj1(s1.ppart,s1.kpic,s1.partt,s1.tdiag)
         if ((s1.nst==1) or (s1.nst==2)):
            s1.partd[itt,:,:] = s1.partt
            itt += 1
         elif (s1.nst==3):
            # calculate particle distribution function and moments
            mdiag1.mvdist1(s1.partt,s1.fvtp,s1.fvmtp,s1.tdiag,s1.nprobt,s1.nmv)
            # display velocity distributions
            graf1.displayfv1(s1.fvtp,s1.fvmtp,' ELECTRON',s1.ntime,s1.nmv,1,s1.irc)
            if (s1.irc[0]==1):
               return
            s1.irc[0] = 0

    # phase space diagnostic
   if (s1.nts > 0):
      s1.it = int(s1.ntime/s1.nts)
      if (s1.ntime==s1.nts*s1.it):
         # plot electrons vx versus x
         pc.showPhase(s1.ppart, s1.kpic, plottype = "EDRAWPHASE")
         graf1.dpmgrasp1(s1.ppart,s1.kpic,' ELECTRON',s1.ntime,999,s1.nx,2,1,
                         s1.ntsc,s1.irc)
         if (s1.irc[0]==1):
            return
         s1.irc[0] = 0
         # ion phase space
         if (s1.movion==1):
            # plot electrons vx versus x
            pc.showPhase(s1.pparti, s1.kipic, plottype = "IDRAWPHASE")
            graf1.dpmgrasp1(s1.pparti,s1.kipic,' ION',s1.ntime,999,s1.nx,2,1,
                            s1.ntsc,s1.irc)
            if (s1.irc[0]==1):
               return
            s1.irc[0] = 0
         
    # push electrons with OpenMP:
   s1.wke[0] = 0.0
    # updates s1.part, s1.wke and possibly s1.ncl, s1.ihole, and s1.irc
   if (s1.mzf==0):
      mpush1.wmpush1(s1.ppart,s1.fxe,s1.kpic,s1.ncl,s1.ihole,s1.qbme,s1.dt,s1.ci,s1.wke,
                     s1.tpush,s1.nx,s1.mx,s1.ipbc,s1.relativity,s1.list,s1.irc)
      # zero force: updates s1.part, s1.wke and possibly s1.ncl, s1.ihole, and s1.irc
   else:
      mpush1.wmpush1zf(s1.ppart,s1.kpic,s1.ncl,s1.ihole,s1.dt,s1.ci,s1.wke,s1.tpush,s1.nx,
                       s1.mx,s1.ipbc,s1.relativity,s1.list,s1.irc)

    # reorder electrons by tile with OpenMP:
    # updates s1.ppart, s1.ppbuff, s1.kpic, s1.ncl, s1.irc, and possibly s1.ihole
   msort1.wmporder1(s1.ppart,s1.ppbuff,s1.kpic,s1.ncl,s1.ihole,s1.tsort,s1.nx,s1.mx,s1.list,
                    s1.irc)

    # sanity check for electrons
   mpush1.mcheck1(s1.ppart,s1.kpic,s1.nx,s1.mx,s1.irc)

    # push ions with OpenMP:
   if (s1.movion==1):
      s1.wki[0] = 0.0
      if (s1.mzf==0):
         # updates s1.pparti, s1.wki and possibly s1.ncl, s1.ihole, and s1.irc
         mpush1.wmpush1(s1.pparti,s1.fxe,s1.kipic,s1.ncl,s1.ihole,s1.qbmi,s1.dt,s1.ci,
                        s1.wki,s1.tpush,s1.nx,s1.mx,s1.ipbc,s1.relativity,s1.list,s1.irc)
        # zero force: updates s1.pparti, s1.wki and possibly s1.ncl, s1.ihole, and s1.irc
      else:
         mpush1.wmpush1zf(s1.pparti,s1.kipic,s1.ncl,s1.ihole,s1.dt,s1.ci,s1.wki,
                          s1.tpush,s1.nx,s1.mx,s1.ipbc,s1.relativity,s1.list,s1.irc)
      s1.wki[0] = s1.wki[0]*s1.rmass

    # reorder ions by tile with OpenMP:
    # updates s1.pparti, s1.ppbuff, s1.kipic, s1.ncl, s1.irc, and possibly s1.ihole
      msort1.wmporder1(s1.pparti,s1.ppbuff,s1.kipic,s1.ncl,s1.ihole,s1.tsort,s1.nx,s1.mx,
                       s1.list,s1.irc)

    # sanity check for ions
      mpush1.mcheck1(s1.pparti,s1.kipic,s1.nx,s1.mx,s1.irc)

    # start running simulation backwards:
    # need to reverse time lag in leap-frog integeration scheme
   if (s1.treverse==1):
      if (((s1.ntime+1)==(s1.nloop/2)) or ((s1.ntime+1)==s1.nloop)):
         dt = -dt
         s1.ws[0] = 0.0
         mpush1.wmpush1zf(s1.ppart,s1.kpic,s1.ncl,s1.ihole,s1.dt,s1.ci,s1.ws,s1.tpush,
                          s1.nx,s1.mx,s1.ipbc,s1.relativity,s1.list,s1.irc)
         msort1.wmporder1(s1.ppart,s1.ppbuff,s1.kpic,s1.ncl,s1.ihole,s1.tsort,s1.nx,s1.mx,
                          s1.list,s1.irc)
         if (s1.movion==1):
            mpush1.wmpush1zf(s1.pparti,s1.kipic,s1.ncl,s1.ihole,s1.dt,s1.ci,s1.ws,
                             s1.tpush,s1.nx,s1.mx,s1.ipbc,s1.relativity,s1.list,
                             s1.irc)
            msort1.wmporder1(s1.pparti,s1.ppbuff,s1.kipic,s1.ncl,s1.ihole,s1.tsort,s1.nx,
                             s1.mx,s1.list,s1.irc)

    # energy diagnostic
   if (s1.ntw > 0):
      s1.it = int(s1.ntime/s1.ntw)
      if (s1.ntime==s1.ntw*s1.it):
         s1.ws[0] = s1.we[0] + s1.wke[0] + s1.wki[0]
         print >> s1.iuot, "Field, Kinetic and Total Energies:"
         if (s1.movion==0):
            s1.iuot.write("%14.7e %14.7e %14.7e\n" % (s1.we[0],s1.wke[0],s1.ws[0])) 
         else: 
            s1.iuot.write("%14.7e %14.7e %14.7e %14.7e\n" % (s1.we[0],s1.wke[0],
                       s1.wki[0],s1.ws[0])) 
         s1.wt[s1.itw,:] = [s1.we[0],s1.wke[0],s1.wki[0],s1.ws[0]]
         pc.showEnergy(numpy.array(range(s1.itw))*in1.dt, s1.wt, s1.itw )
         s1.itw += 1
         s1.s[0] += s1.we[0]
         s1.s[1] += s1.wke[0]

#write results to file
def write(s1):
    print >> s1.iuot
    print >> s1.iuot, "s1.ntime, relativity = ", s1.ntime, ",", s1.relativity
    if (s1.treverse==1):
       print >> s1.iuot, "treverse = ", s1.treverse

    if ((s1.ntw > 0) or (s1.ntt > 0)):
       graf1.reset_graphs()

    # velocity diagnostic
    if (s1.ntv > 0):
       ts = s1.t0 + s1.dt*float(s1.ntv)
       graf1.displayfvt1(s1.fvtm,' ELECTRON',ts,s1.dt*float(s1.ntv),s1.itv,s1.irc)
       if (s1.irc[0]==1):
          exit(1)
    # ions
       if (s1.movion==1):
          graf1.displayfvt1(s1.fvtmi,' ION',ts,s1.dt*float(s1.ntv),s1.itv,s1.irc)
          if (s1.irc[0]==1):
             exit(1)

    # trajectory diagnostic
    if (s1.ntt > 0):
       if ((s1.nst==1) or (s1.nst==2)):
          if (s1.nplot > 0):
             s1.irc[0] = graf1.open_graphs(1)
          ts = s1.t0 + s1.dt*float(s1.ntt)
          graf1.displayt1(s1.partd,ts,s1.dt*float(s1.ntt),itt,2,3,s1.irc)
          if (s1.irc[0]==1):
             exit(1)

    # energy diagnostic
    if (s1.ntw > 0):
       s1.irc[0] = graf1.open_graphs(s1.nplot)
       ts = s1.t0 + s1.dt*float(s1.ntw)
       graf1.displayw1(s1.wt,ts,s1.dt*float(s1.ntw),s1.itw,s1.irc)
       if (s1.irc[0]==1):
          exit(1)
       swe = s1.s[0]; swke = s1.s[1]
       swe = swe/float(s1.itw)
       print >> s1.iuot, "Average Field Energy <WE> = ", float(swe)
       swke = swke/float(s1.itw)
       print >> s1.iuot, "Average Electron Kinetic Energy <WKE> = ",float(swke)
       print >> s1.iuot, "Ratio <WE>/<WKE>= ", float(swe/swke)

    # display final spectral analysis for ion density
    if (s1.ntdi > 0):
       if ((s1.nddi==2) or (s1.nddi==3)):
    # performs frequency analysis of accumulated complex time series
          s1.wkdi[:,0] = s1.wm[numpy.argmax(s1.pkwdi[:,:,0],axis=1)]
          s1.wkdi[:,1] = s1.wm[numpy.argmax(s1.pkwdi[:,:,1],axis=1)]
    # display frequency spectrum
          graf1.dmscaler1(s1.wkdi,'IDENSITY OMEGA VS MODE',s1.ntime,999,1,
                          s1.modesxdi,s1.cwk,s1.irc)

    # display final spectral analysis for potential
    if (s1.ntp > 0):
       if ((s1.ndp==2) or (s1.ndp==3)):
    # performs frequency analysis of accumulated complex time series
          s1.wk[:,0] = s1.wm[numpy.argmax(s1.pkw[:,:,0],axis=1)]
          s1.wk[:,1] = s1.wm[numpy.argmax(s1.pkw[:,:,1],axis=1)]
    # display frequency spectrum
          graf1.dmscaler1(s1.wk,'POTENTIAL OMEGA VS MODE',s1.ntime,999,1,
                          s1.modesxp,s1.cwk,s1.irc)

    print >> s1.iuot
    print >> s1.iuot, "initialization time = ", s1.tinit
    print >> s1.iuot, "deposit time = ", s1.tdpost[0]
    print >> s1.iuot, "guard time = ", s1.tguard[0]
    print >> s1.iuot, "solver time = ", s1.tfield[0]
    print >> s1.iuot, "fft time = ", s1.tfft[0]
    print >> s1.iuot, "push time = ", s1.tpush[0]
    print >> s1.iuot, "sort time = ", s1.tsort[0]
    s1.tfield[0] = s1.tfield[0] + s1.tguard[0] + s1.tfft[0]
    print >> s1.iuot, "total solver time = ", s1.tfield[0]
    time = s1.tdpost[0] + s1.tpush[0] + s1.tsort[0]
    print >> s1.iuot, "total particle time = ", time
    s1.ws[0] = time + s1.tfield[0] + s1.tdiag[0]
    s1.tloop = s1.tloop - s1.ws[0]
    print >> s1.iuot, "total and additional time = ", s1.ws[0], ",", s1.tloop
    print >> s1.iuot

    s1.ws[0] = 1.0e+09/(float(s1.nloop)*float(s1.np))
    print >> s1.iuot, "Push Time (nsec) = ", s1.tpush[0]*s1.ws[0]
    print >> s1.iuot, "Deposit Time (nsec) = ", s1.tdpost[0]*s1.ws[0]
    print >> s1.iuot, "Sort Time (nsec) = ", s1.tsort[0]*s1.ws[0]
    print >> s1.iuot, "Total Particle Time (nsec) = ", time*s1.ws[0]
    print >> s1.iuot

    s1.writnml1(s1.iudm)
    print >> s1.iuot, " * * * q.e.d. * * *"
    s1.iuot.close()
    # close graphics device
    graf1.close_graphs()

#init GUI
pc = PlasmaContext()  #Create GUI
pc.showGraphs(True)   #enable graphics.  Setting to false will disable graphics
pc.asyncMode(0)     #Run in synchornos mode.
pc.callbacks["VARCHANGE"] = changeVarsCallback  #Set a callback
pc.callbacks["RESET"] = resetCallback  #define the callback
pc.clearGraphList()  #remove all default graph options

#Now init Fortran simulation
init(in1)

#Begin main loop
for ntime in xrange(0,in1.nloop):
    in1.ntime = ntime
    print >> in1.iuot, "ntime = ", in1.ntime
    step(in1)

in1.ntime = in1.ntime + 1
# loop time
dtimer(in1.dtime,in1.ltime,1)
in1.tloop = in1.tloop + float(in1.dtime)

#write to file
write(in1)

#Wait in a loop to maintain interactivity
pc.wait()
