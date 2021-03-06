!-----------------------------------------------------------------------
! Fortran Library for diagnostics
! 1D OpenMP PIC Codes:
! CSPECT1 performs frequency analysis of complex time series
! ICSPECT1 performs incremental frequency analysis of complex time
!          series for one time step
! VPDIST1 calculates 1 component velocity distribution, velocity moments,
!         and entropy, for segmented particle array
! VDIST1 calculates 1 component velocity distribution, velocity moments,
!        and entropy for standard particle array
! STPTRAJ1 sets test charge distribution by setting a particle id
!          in particle location 3 for 1d code
! STPTRAJ13 sets test charge distribution by setting a particle id
!           in particle location 5 for 1-2/2d code
! PTRAJ1 copies tagged particles in ppart to array partt for 1d code
! STPBEAM1 for 1d code, marks beam particles by setting a particle id in
!          particle location 3
! STPBEAM13 for 1-2/2d code, marks beam particles by setting a particle
!           id in particle location 5
! written by viktor k. decyk, ucla
! copyright 2016, regents of the university of california
! update: october 6, 2016
!-----------------------------------------------------------------------
      subroutine CSPECT1(fc,wm,pkw,t0,dt,nt,iw,modesx,ntd,iwd,modesxd)
! this subroutine performs frequency analysis of complex time series,
! pkw(w,k,1:2) = |(1/nt)*sum {fc(t,k)*exp(sqrt(-1)*w*(t-t0))}|**2
! where wm(j) stores the frequency values w
! it is an sft (slow fourier transform), but you can pick your frequency
! on input, fc contains the data to be analyzed, real and imaginary
! parts stored adjacent, and wm(w) contains the (positive) frequencies.
! on output, pkw(:,:,1) contains result for positive frequencies,
! pkw(:,:,2) the negative frequencies.
! t0 = starting time value
! dt = time step
! nt = number of input data points, 
! iw = number of (positive) frequencies
! modesx = number of modes in x direction
! ntd = dimension of input array, ntd >= nt
! iwd = dimension of frequency array iw >= iw
! modesxd = dimension of output array pott, modesxd >= modesx
      implicit none
      integer nt, iw, modesx, ntd, iwd, modesxd
      real t0, dt
      real wm, pkw
      dimension wm(iwd), pkw(modesxd,iwd,2)
      complex fc
      dimension fc(ntd,modesxd)
! local data
      integer i, j, k
      real anl, fr, fi
      double precision at1, at2, at3, sum1, sum2, sum3, sum4, cwdt, swdt
      anl = 1.0/real(nt)
! loop over frequencies
      do 30 j = 1, iw
      at3 = wm(j)*dt
      cwdt = dcos(at3)
      swdt = dsin(at3)
! loop over modes
      do 20 k = 1, modesx
      at3 = wm(j)*t0
      at1 = dcos(at3)
      at2 = -dsin(at3)
      sum1 = 0.0d0
      sum2 = 0.0d0
      sum3 = 0.0d0
      sum4 = 0.0d0
! loop over time
      do 10 i = 1, nt
      fr = real(fc(i,k))
      fi = aimag(fc(i,k))
      sum1 = sum1 + fr*at1
      sum2 = sum2 + fi*at1
      sum3 = sum3 + fi*at2
      sum4 = sum4 + fr*at2
      at3 = at1*cwdt - at2*swdt
      at2 = at2*cwdt + at1*swdt
      at1 = at3
   10 continue
      at1 = anl*(sum1 - sum3)
      at2 = anl*(sum2 + sum4)
      pkw(k,j,1) = at1*at1 + at2*at2
      at1 = anl*(sum1 + sum3)
      at2 = anl*(sum2 - sum4)
      pkw(k,j,2) = at1*at1 + at2*at2
   20 continue
   30 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine ICSPECT1(fc,wm,pkw,pks,time,t0,nt,iw,modesx,nx,norm,iwd&
     &,modesxd)
! this subroutine performs incremental frequency analysis of complex
! time series for one time step
! pkw(w,k,1:2) = |(anorm/nt)*sum {fc(k)*exp(sqrt(-1)*w*(time-t0)}|**2
! where wm(j) stores the frequency values w
! it is an sft (slow fourier transform), but you can pick your frequency
! on input, fc contains the data to be analyzed for one time step,
! real and imaginary parts stored adjacent, and
! wm(w) contains the (positive) frequencies.
! on output, pkw(:,:,1) contains result for positive frequencies,
! pkw(:,:,2) the negative frequencies.
! pks = accumulated complex spectrum up to current time,
! should be initialized to zero
! time = current time value
! t0 = starting time value
! itw current time index
! nt = number of input data points (for normalization)
! iw = number of (positive) frequencies
! modesx = number of modes in x direction
! nx = system length in x direction
! norm = (-1,0,1) = normalize with (inverse gradient,null,gradient) op
! norm = 1 for potential or norm = -1 for density gives spectrum as
!        electric field energy
! ntd = dimension of input array, ntd >= nt
! iwd = dimension of frequency array iw >= iw
! modesxd = dimension of output array pott, modesxd >= modesx
      implicit none
      integer nt, iw, modesx, nx, norm, iwd, modesxd
      real time, t0
      real wm, pkw
      dimension wm(iwd), pkw(modesxd,iwd,2)
      complex fc
      dimension fc(modesxd)
      double precision pks
      dimension pks(4,modesxd,iwd)
! local data
      integer j, k
      real anl, dnx, anorm, fr, fi
      double precision at1, at2, at3, sum1, sum2, sum3, sum4
      anl = 1.0/real(nt)
      dnx = 6.28318530717959/real(nx)
      anorm = anl
! loop over frequencies
      do 20 j = 1, iw
! loop over modes
      do 10 k = 1, modesx
      at3 = wm(j)*(time - t0)
      at1 = dcos(at3)
      at2 = dsin(at3)
! add contribution for current time
      fr = real(fc(k))
      fi = aimag(fc(k))
      sum1 = pks(1,k,j) + fr*at1
      sum2 = pks(2,k,j) + fi*at1
      sum3 = pks(3,k,j) + fi*at2
      sum4 = pks(4,k,j) + fr*at2
! save accumulation for next time
      pks(1,k,j) = sum1
      pks(2,k,j) = sum2
      pks(3,k,j) = sum3
      pks(4,k,j) = sum4
! normalize
      if (k.gt.1) then
         if (norm==1) then
            anorm = anl*(dnx*real(k - 1))
         else if (norm==(-1)) then
            anorm = anl/(dnx*real(k - 1))
         endif
      endif
! calculate spectrum for accumulated data
      at1 = anorm*(sum1 - sum3)
      at2 = anorm*(sum2 + sum4)
      pkw(k,j,1) = at1*at1 + at2*at2
      at1 = anorm*(sum1 + sum3)
      at2 = anorm*(sum2 - sum4)
      pkw(k,j,2) = at1*at1 + at2*at2
   10 continue
   20 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine VPDIST1(ppart,kpic,sfv,fvm,idimp,nppmx,mx1,np,nmv,nmvf)
! for 1d code, this subroutine calculates 1d velocity distribution,
! velocity moments, and entropy
! particles stored segmented array
! input: all except fvm, output: fv, fvm
! ppart(2,n,m) = velocity vx of particle n in tile m
! kpic = number of particles per tile
! sfv = distribution function particles in each velocity range in tile
! maximum velocity (used for scaling) is contained in first element sfv.
! vdrift is contained in fvm(1)
! vth is contained in fvm(2)
! entropy is contained in fvm(3), defined to be:
! s/k = -sum(f(v)/np)*log(f(v)/(np*delta_v*delta_x)).
! idimp = size of phase space = 2
! nppmx = maximum number of particles in tile
! mx1 = (system length in x direction - 1)/mx + 1
! np = number of particles
! nmvf = dimension of fv
! the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer idimp, nppmx, mx1, np, nmv, nmvf
      real ppart, sfv, fvm
      dimension ppart(idimp,nppmx,mx1)
      dimension sfv(nmvf,mx1+1), fvm(3)
      integer kpic
      dimension kpic(mx1)
! local data
      integer j, k, npp, nvx
      real anmv, svx, svxx, vx
      double precision sumvx, sumvx2, ssumvx, ssumvx2, anp, sum1
! velocity scaling, same scaling used for all tiles
      anmv = real(nmv)
      svx = anmv/sfv(1,mx1+1)
! normalization constant for entropy
      svxx = svx*real(mx1)
! zero out distribution
      do 20 k = 1, mx1+1
      do 10 j = 2, nmvf
      sfv(j,k) = 0.0
   10 continue
   20 continue
! count particles in each velocity region
      anmv = anmv + 2.5
      sumvx = 0.0d0
      sumvx2 = 0.0d0
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(j,k,npp,nvx,vx,ssumvx,ssumvx2)
!$OMP& REDUCTION(+:sumvx) REDUCTION(+:sumvx2)
      do 40 k = 1, mx1
      npp = kpic(k)
      ssumvx = 0.0d0
      ssumvx2 = 0.0d0
! loop over particles in tile
      do 30 j = 1, npp
      vx = ppart(2,j,k)
      nvx = vx*svx + anmv
      ssumvx = ssumvx + vx
      ssumvx2 = ssumvx2 + vx*vx
      if ((nvx.ge.2).and.(nvx.le.nmvf)) sfv(nvx,k) = sfv(nvx,k) + 1.0
   30 continue
! calculate global sums
      sumvx = sumvx + ssumvx
      sumvx2 = sumvx2 + ssumvx2
   40 continue
!$OMP END PARALLEL DO
! calculate global distribution
      do 60 j = 2, nmvf
      sum1 = 0.0d0
      do 50 k = 1, mx1
      sum1 = sum1 + sfv(j,k)
   50 continue
      sfv(j,mx1+1) = sum1
   60 continue
! calculate global velocity moments
      anp = 0.0d0
      if (np.gt.0) anp = 1.0d0/dble(np)
      sumvx = sumvx*anp
      fvm(1) = sumvx
      fvm(2) = dsqrt(sumvx2*anp - sumvx**2)
! count number of particles in global distribution
      sumvx = 0.0d0
      do 70 j = 2, nmvf
      sumvx = sumvx + sfv(j,mx1+1)
   70 continue
! calculate entropy
      sumvx2 = 0.0d0
      do 90 k = 1, mx1
      do 80 j = 2, nmvf
      if (sfv(j,k).gt.0.0) then
         sumvx2 = sumvx2 + sfv(j,k)*dlog(dble(sfv(j,k)*svxx))
      endif
   80 continue
   90 continue
      if (sumvx.gt.0.0d0) sumvx = -sumvx2/sumvx + dlog(sumvx)
      fvm(3) = sumvx
      return
      end
!-----------------------------------------------------------------------
      subroutine VDIST1(part,fv,fvm,idimp,np,nmv,nmvf)
! for 1d code, this subroutine calculates 1d velocity distribution,
! velocity moments, and entropy
! input: all except fvm, output: fv, fvm
! part(2,n) = velocity vx of particle n
! fv = distribution function, number of particles in each velocity range
! maximum velocity (used for scaling) is contained in first element fv.
! vdrift is contained in fvm(1)
! vth is contained in fvm(2)
! entropy is contained in fvm(3), defined to be:
! s/k = -sum(f(v)/np)*log(f(v)/(np*delta_v)).  Assumes that distribution
! is uniform in space
! idimp = size of phase space = 2
! np = number of particles
! the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer idimp, np, nmv, nmvf
      real part, fv, fvm
      dimension part(idimp,np), fv(nmvf), fvm(3)
! local data
      double precision sumvx, sumvx2, anp
      real anmv, svx
      integer j, nvx
      anmv = real(nmv)
      svx = anmv/fv(1)
! zero out distribution
      do 10 j = 2, nmvf
      fv(j) = 0.0
   10 continue
! count particles in each velocity region
      anmv = anmv + 2.5
      sumvx = 0.0d0
      sumvx2 = 0.0d0
      do 20 j = 1, np
      nvx = part(2,j)*svx + anmv
      sumvx = sumvx + part(2,j)
      sumvx2 = sumvx2 + part(2,j)**2
      if ((nvx.ge.2).and.(nvx.le.nmvf)) fv(nvx) = fv(nvx) + 1.0
   20 continue
! calculate velocity moments
      anp = 1.0d0/dble(np)
      sumvx = sumvx*anp
      fvm(1) = sumvx
      fvm(2) = dsqrt(sumvx2*anp - sumvx**2)
! calculate entropy
      sumvx = 0.0d0
      sumvx2 = 0.0d0
      do 30 j = 2, nmvf
      if (fv(j).gt.0.0) then
         sumvx = sumvx + fv(j)
         sumvx2 = sumvx2 + fv(j)*dlog(dble(fv(j)*svx))
      endif
   30 continue
      if (sumvx.gt.0.0d0) sumvx = -sumvx2/sumvx + dlog(sumvx)
      fvm(3) = sumvx
      return
      end
!-----------------------------------------------------------------------
      subroutine STPTRAJ1(ppart,kpic,iprobt,nst,vtx,vtsx,dvtx,idimp,    &
     &nppmx,mx1,np,nprobt)
! for 1d code, this procedure sets test charge distribution by setting
! a particle id in particle location 3, whose values are between 1 and
! nprobt
! particles stored segmented array
! input: all, output: iprobt, nprobt
! ppart(2,n,m) = velocity vx of particle n in tile m
! ppart(3,n,m) = particle id of tagged particle n in tile m
! kpic = number of particles per tile
! iprobt = scratch array of size nprobt, used for nst = 2
! nst = type of test particle distribution
!   1 = uniformly distribution in real space
!   2 = uniform distribution in velocity space
!   3 = velocity slice at vtsx +- dvtx/2
! vtx = thermal velocity of particles in x direction, if nst = 2
! vtsx = center of velocity slice if nst = 3
! dvtx = width of velocity slice if nst = 3
! idimp = size of phase space = 3
! nppmx = maximum number of particles in tile
! mx1 = (system length in x direction - 1)/mx + 1
! np = total number of particles in part
! nprobt = number of test charges whose trajectories will be stored.
! particle id should be <= 16777215
      implicit none
      integer nst, idimp, nppmx, mx1, np, nprobt
      real vtx, vtsx, dvtx
      real ppart
      dimension ppart(idimp,nppmx,mx1)
      integer kpic, iprobt
      dimension kpic(mx1)
      dimension iprobt(nprobt)
! local data
      integer j, k, npp, joff, it, nt, itt
      real st, at
      if (idimp < 3) return
! set up constants
      itt = 0; at = 0.0; st = 0.0
! uniform distribution in number space
      if (nst.eq.1) then
         it = np/nprobt
         itt = 1
! uniform distribution in velocity space in x direction
      else if (nst.eq.2) then
         st = 4.0*vtx
         it = nprobt/2
         at = real(it)
         st = at/st
         at = at + 1.0 + 0.5*(nprobt - 2*it)
         do 10 j = 1, nprobt
         iprobt(j) = 0
   10    continue
! velocity slice in x direction
      else if (nst.eq.3) then
         st = 1.0/dvtx
         itt = vtsx*st + 0.5
      endif
      joff = 0
      nt = 0
! loop over tiles
      do 30 k = 1, mx1
      npp = kpic(k)
! loop over particles in tile
      do 20 j = 1, npp
! clear tag
      ppart(3,j,k) = 0.0
! uniform distribution in number space
      if (nst.eq.1) then
         if ((j+joff).eq.itt) then
            nt = nt + 1
            ppart(3,j,k) = real(nt)
            itt = itt + it
         endif
! uniform distribution in velocity space in x direction
      else if (nst.eq.2) then
         it = ppart(2,j,k)*st + at
         if ((it.gt.0).and.(it.le.nprobt)) then
            if (iprobt(it).eq.0) then
               nt = nt + 1
               iprobt(it) = j + joff
               ppart(3,j,k) = real(nt)
            endif
         endif
! velocity slice in x direction
      else if (nst.eq.3) then
         it = ppart(2,j,k)*st + 0.5
         if (it.eq.itt) then
            nt = nt + 1
            ppart(3,j,k) = real(nt)
         endif
      endif
   20 continue
      joff = joff + npp
   30 continue
      nprobt = nt
      return
      end
!-----------------------------------------------------------------------
      subroutine STPTRAJ13(ppart,kpic,iprobt,nst,vtx,vtsx,dvtx,idimp,   &
     &nppmx,mx1,np,nprobt)
! for 1-2/2d code, this procedure sets test charge distribution by
! setting a particle id in particle location 5, whose values are between
! 1 and nprobt
! particles stored segmented array
! input: all, output: iprobt, nprobt
! ppart(2,n,m) = velocity vx of particle n in tile m
! ppart(5,n,m) = particle id of tagged particle n in tile m
! kpic = number of particles per tile
! iprobt = scratch array of size nprobt, used for nst = 2
! nst = type of test particle distribution
!   1 = uniformly distribution in real space
!   2 = uniform distribution in velocity space
!   3 = velocity slice at vtsx +- dvtx/2
! vtx = thermal velocity of particles in x direction, if nst = 2
! vtsx = center of velocity slice if nst = 3
! dvtx = width of velocity slice if nst = 3
! idimp = size of phase space = 5
! nppmx = maximum number of particles in tile
! mx1 = (system length in x direction - 1)/mx + 1
! np = total number of particles in part
! nprobt = number of test charges whose trajectories will be stored.
! particle id should be <= 16777215
      implicit none
      integer nst, idimp, nppmx, mx1, np, nprobt
      real vtx, vtsx, dvtx
      real ppart
      dimension ppart(idimp,nppmx,mx1)
      integer kpic, iprobt
      dimension kpic(mx1)
      dimension iprobt(nprobt)
! local data
      integer j, k, npp, joff, it, nt, itt
      real st, at
      if (idimp < 5) return
! set up constants
      itt = 0; at = 0.0; st = 0.0
! uniform distribution in number space
      if (nst.eq.1) then
         it = np/nprobt
         itt = 1
! uniform distribution in velocity space in x direction
      else if (nst.eq.2) then
         st = 4.0*vtx
         it = nprobt/2
         at = real(it)
         st = at/st
         at = at + 1.0 + 0.5*(nprobt - 2*it)
         do 10 j = 1, nprobt
         iprobt(j) = 0
   10    continue
! velocity slice in x direction
      else if (nst.eq.3) then
         st = 1.0/dvtx
         itt = vtsx*st + 0.5
      endif
      joff = 0
      nt = 0
! loop over tiles
      do 30 k = 1, mx1
      npp = kpic(k)
! loop over particles in tile
      do 20 j = 1, npp
! clear tag
      ppart(5,j,k) = 0.0
! uniform distribution in number space
      if (nst.eq.1) then
         if ((j+joff).eq.itt) then
            nt = nt + 1
            ppart(5,j,k) = real(nt)
            itt = itt + it
         endif
! uniform distribution in velocity space in x direction
      else if (nst.eq.2) then
         it = ppart(2,j,k)*st + at
         if ((it.gt.0).and.(it.le.nprobt)) then
            if (iprobt(it).eq.0) then
               nt = nt + 1
               iprobt(it) = j + joff
               ppart(5,j,k) = real(nt)
            endif
         endif
! velocity slice in x direction
      else if (nst.eq.3) then
         it = ppart(2,j,k)*st + 0.5
         if (it.eq.itt) then
            nt = nt + 1
            ppart(5,j,k) = real(nt)
         endif
      endif
   20 continue
      joff = joff + npp
   30 continue
      nprobt = nt
      return
      end
!-----------------------------------------------------------------------
      subroutine PTRAJ1(ppart,kpic,partt,idimp,nppmx,mx1,nprobt)
! this procedure copies tagged particles in ppart to array partt
! input: all except partt, output: partt
! ppart(2,n,m) = velocity vx of particle n in tile m
! ppart(3,n,m) = particle id of tagged particle n in tile m
! kpic = number of particles per tile
! idimp = size of phase space = 3
! nppmx = maximum number of particles in tile
! mx1 = (system length in x direction - 1)/mx + 1
! nprobt = number of test charges whose trajectories will be stored.
! particle id should be <= 16777215
      implicit none
      integer idimp, nppmx, mx1, nprobt
      real ppart, partt
      dimension ppart(idimp,nppmx,mx1), partt(idimp,nprobt)
      integer kpic
      dimension kpic(mx1)
! local data
      integer i, j, k, npp, nt
      real tn
! loop over tiles
!$OMP PARALLEL DO
!$OMP& PRIVATE(i,j,k,npp,tn)
      do 30 k = 1, mx1
      npp = kpic(k)
! loop over particles in tile
      do 20 j = 1, npp
      tn = ppart(3,j,k)
      if (tn.gt.0.0) then
         nt = tn
         do 10 i = 1, idimp
         partt(i,nt) = ppart(i,j,k)
   10    continue
      endif
   20 continue
   30 continue
!$OMP END PARALLEL DO
      return
      end
!-----------------------------------------------------------------------
      subroutine STPBEAM1(part,npx,idimp,nop)
! for 1d code, this procedure marks initial beam distribution by setting
! a particle id in particle location 3, whose values are negative for
! particle n if n > npx, zero otherwise
! part(3,n) = particle id of tagged particle n
! npx = number of background particles distributed in x direction
! idimp = size of phase space = 3
! nop = number of particles
      implicit none
      integer npx, idimp, nop
      real part
      dimension part(idimp,nop)
! local data
      integer j, js
      if (idimp < 3) return
! zero tag for initial background particles
      do 10 j = 1, npx
      part(3,j) = 0.0
   10 continue
! negative tag for initial beam particles
      js = npx + 1
      do 20 j = js, nop
      part(3,j) = -1.0
   20 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine STPBEAM13(part,npx,idimp,nop)
! for 1-2/2d code, this procedure marks initial beam distribution by
! setting a particle id in particle location 3, whose values are
! negative for particle n if n > npx, zero otherwise
! part(3,n) = particle id of tagged particle n
! npx = number of background particles distributed in x direction
! idimp = size of phase space = 5
! nop = number of particles
      implicit none
      integer npx, idimp, nop
      real part
      dimension part(idimp,nop)
! local data
      integer j, js
      if (idimp < 5) return
! zero tag for initial background particles
      do 10 j = 1, npx
      part(5,j) = 0.0
   10 continue
! negative tag for initial beam particles
      js = npx + 1
      do 20 j = js, nop
      part(5,j) = -1.0
   20 continue
      return
      end
