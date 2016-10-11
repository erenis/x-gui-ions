c-----------------------------------------------------------------------
c 1d PIC library for diagnostics
c diag1lib.f contains diagnostic procedures:
c VDIST1 calculates 1 component velocity distribution, velocity moments,
c        and entropy.
c VDIST13 calculates 3 component velocity distribution, velocity
c         moments, and entropy.
c PSDIST1 calculates 1d phase space distribution, velocity moments, and
c         entropy.
c PSDIST13 calculates 3d phase space distribution, velocity moments, and
c          entropy.
c FWRITE1 writes real binary data to file
c FCWRITE1 writes complex binary data to file.
c FREAD1 writes real binary data to file.
c FCREAD1 writes complex binary data to file.
c written by viktor k. decyk, ucla
c copyright 1994, regents of the university of california
c update: july 1, 2013
c-----------------------------------------------------------------------
      subroutine VDIST1(part,fv,fvm,idimp,np,nmv,nmvf)
c for 1d code, this subroutine calculates 1d velocity distribution,
c velocity moments, and entropy
c input: all except fvm, output: fv, fvm
c part(2,n) = velocity vx of particle n
c fv = distribution function, number of particles in each velocity range
c maximum velocity (used for scaling) is contained in first element fv.
c vdrift is contained in fvm(1)
c vth is contained in fvm(2)
c entropy is contained in fvm(3), defined to be:
c s/k = -sum(f(v)/np)*log(f(v)/(np*delta_v)).  Assumes that distribution
c is uniform in space
c idimp = size of phase space = 2
c np = number of particles
c the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer idimp, np, nmv, nmvf
      real part, fv, fvm
      dimension part(idimp,np), fv(nmvf), fvm(3)
c local data
      double precision sumvx, sumvx2, anp
      real anmv, svx
      integer j, nvx
      anmv = real(nmv)
      svx = anmv/fv(1)
c zero out distribution
      do 10 j = 2, nmvf
      fv(j) = 0.0
   10 continue
c count particles in each velocity region
      anmv = anmv + 2.5
      sumvx = 0.0d0
      sumvx2 = 0.0d0
      do 20 j = 1, np
      nvx = part(2,j)*svx + anmv
      sumvx = sumvx + part(2,j)
      sumvx2 = sumvx2 + part(2,j)**2
      if ((nvx.ge.2).and.(nvx.le.nmvf)) fv(nvx) = fv(nvx) + 1.0
   20 continue
c calculate velocity moments
      anp = 1.0d0/dble(np)
      sumvx = sumvx*anp
      fvm(1) = sumvx
      fvm(2) = dsqrt(sumvx2*anp - sumvx**2)
c calculate entropy
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
      subroutine VDIST13(part,fv,fvm,idimp,np,nmv,nmvf)
c for 1-2/2d code, this subroutine calculates 3d velocity distribution,
c velocity moments, and entropy
c input: all except fvm, output: fv, fvm
c part(2,n) = velocity vx of particle n
c part(3,n) = velocity vy of particle n
c part(4,n) = velocity vz of particle n
c fv = distribution function, number of particles in each velocity range
c maximum velocity (used for scaling) is contained in first element fv.
c vdrift for i-th dimension is contained in fvm(1,i)
c vth for i-th dimension is contained in fvm(2,i)
c entropy for i-th dimension is contained in fvm(3,i), defined to be:
c s/k = -sum(f(v)/np)*log(f(v)/(np*delta_v)).  Assumes that distribution
c is uniform in space and distributions in each dimension are
c independent.
c idimp = size of phase space = 4
c np = number of particles
c the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer idimp, np, nmv, nmvf
      real part, fv, fvm
      dimension part(idimp,np), fv(nmvf,3), fvm(3,3)
c local data
      double precision sumvx, sumvy, sumvz, sumvx2, sumvy2, sumvz2, anp
      real anmv, svx, svy, svz
      integer j, nvx, nvy, nvz
      anmv = real(nmv)
      svx = anmv/fv(1,1)
      svy = anmv/fv(1,2)
      svz = anmv/fv(1,3)
c zero out distribution
      do 10 j = 2, nmvf
      fv(j,1) = 0.0
      fv(j,2) = 0.0
      fv(j,3) = 0.0
   10 continue
c count particles in each velocity region
      anmv = anmv + 2.5
      sumvx = 0.0d0
      sumvy = 0.0d0
      sumvz = 0.0d0
      sumvx2 = 0.0d0
      sumvy2 = 0.0d0
      sumvz2 = 0.0d0
      do 20 j = 1, np
      nvx = part(2,j)*svx + anmv
      sumvx = sumvx + part(2,j)
      sumvx2 = sumvx2 + part(2,j)**2
      nvy = part(3,j)*svy + anmv
      sumvy = sumvy + part(3,j)
      sumvy2 = sumvy2 + part(3,j)**2
      nvz = part(4,j)*svz + anmv
      sumvz = sumvz + part(4,j)
      sumvz2 = sumvz2 + part(4,j)**2
      if ((nvx.ge.2).and.(nvx.le.nmvf)) fv(nvx,1) = fv(nvx,1) + 1.0
      if ((nvy.ge.2).and.(nvy.le.nmvf)) fv(nvy,2) = fv(nvy,2) + 1.0
      if ((nvz.ge.2).and.(nvz.le.nmvf)) fv(nvz,3) = fv(nvz,3) + 1.0
   20 continue
c calculate velocity moments
      anp = 1.0d0/dble(np)
      sumvx = sumvx*anp
      fvm(1,1) = sumvx
      fvm(2,1) = dsqrt(sumvx2*anp - sumvx**2)
      sumvy = sumvy*anp
      fvm(1,2) = sumvy
      fvm(2,2) = dsqrt(sumvy2*anp - sumvy**2)
      sumvz = sumvz*anp
      fvm(1,3) = sumvz
      fvm(2,3) = dsqrt(sumvz2*anp - sumvz**2)
c calculate entropy
      sumvx = 0.0d0
      sumvy = 0.0d0
      sumvz = 0.0d0
      sumvx2 = 0.0d0
      sumvy2 = 0.0d0
      sumvz2 = 0.0d0
      do 30 j = 2, nmvf
      if (fv(j,1).gt.0.0) then
         sumvx = sumvx + fv(j,1)
         sumvx2 = sumvx2 + fv(j,1)*dlog(dble(fv(j,1)*svx))
      endif
      if (fv(j,2).gt.0.0) then
         sumvy = sumvy + fv(j,2)
         sumvy2 = sumvy2 + fv(j,2)*dlog(dble(fv(j,2)*svy))
      endif
      if (fv(j,3).gt.0.0) then
         sumvz = sumvz + fv(j,3)
         sumvz2 = sumvz2 + fv(j,3)*dlog(dble(fv(j,3)*svz))
      endif
   30 continue
      if (sumvx.gt.0.0d0) sumvx = -sumvx2/sumvx + dlog(sumvx)
      if (sumvy.gt.0.0d0) sumvy = -sumvy2/sumvy + dlog(sumvy)
      if (sumvz.gt.0.0d0) sumvz = -sumvz2/sumvz + dlog(sumvz)
      fvm(3,1) = sumvx
      fvm(3,2) = sumvy
      fvm(3,3) = sumvz
      return
      end
      subroutine PSDIST1(part,fps,fpsm,psm,nx,nxb,idimp,np,nmv,nmvf)
c for 1d code, this subroutine calculates 1d phase space distribution,
c velocity moments, and entropy
c input: all except fpsm and psm, output: fps, fpsm, psm
c part(1,n) = position x of particle n
c part(2,n) = velocity vx of particle n
c fps = distribution function, number of particles in each velocity
c and spatial range
c maximum velocity (used for scaling) is contained in first element fps.
c vdrift is contained in fpsm(1)
c vth is contained in fpsm(2)
c entropy is contained in fpsm(3), defined to be:
c s/k = -sum(f(v)/np)*log(f(v)/(np*delta_v)).
c psm = double precision scratch array
c nx = system length in x direction
c nxb = number of spatial regions for distribution
c idimp = size of phase space = 2
c np = number of particles
c the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer nx, nxb, idimp, np, nmv, nmvf
      real part, fps, fpsm
      dimension part(idimp,np), fps(nmvf,nxb), fpsm(3,nxb)
      double precision psm
      dimension psm(2,nxb)
c local data
      double precision sumvx, sumvx2, anp
      real anmv, svx, at1
      integer j, k, n, nvx
      anmv = real(nmv)
c zero out distribution
      do 20 k = 1, nxb
      do 10 j = 2, nmvf
      fps(j,k) = 0.0
   10 continue
      fpsm(1,k) = anmv/fps(1,k)
      psm(1,k) = 0.0d0
      psm(2,k) = 0.0d0
   20 continue
c count particles in each velocity and spatial region
      at1 = real(nxb)/real(nx)
      anmv = anmv + 2.5
      do 30 j = 1, np
      n = part(1,j)*at1
      n = n + 1
      if ((n.ge.1).and.(n.le.nxb)) then
         nvx = part(2,j)*fpsm(1,n) + anmv
         psm(1,n) = psm(1,n) + part(2,j)
         psm(2,n) = psm(2,n) + part(2,j)**2
         if ((nvx.ge.2).and.(nvx.le.nmvf)) fps(nvx,n) = fps(nvx,n) + 1.0
      endif
   30 continue
c calculate velocity moments
      anp = 1.0d0/dble(np)
      do 50 k = 1, nxb
      svx = fpsm(1,k)
      sumvx = psm(1,k)*anp
      fpsm(1,k) = sumvx
      fpsm(2,k) = dsqrt(psm(2,k)*anp - sumvx**2)
c calculate entropy
      sumvx = 0.0d0
      sumvx2 = 0.0d0
      do 40 j = 2, nmvf
      if (fps(j,k).gt.0.0) then
         sumvx = sumvx + fps(j,k)
         sumvx2 = sumvx2 + fps(j,k)*dlog(dble(fps(j,k)*svx))
      endif
   40 continue
      if (sumvx.gt.0.0d0) sumvx = -sumvx2/sumvx + dlog(sumvx)
      fpsm(3,k) = sumvx
   50 continue
      return
      end
      subroutine PSDIST13(part,fps,fpsm,psm,nx,nxb,idimp,np,nmv,nmvf)
c for 1-2/2d code, this subroutine calculates 3d phase space
c distribution, velocity moments, and entropy
c input: all except fpsm and psm, output: fps, fpsm, psm
c part(1,n) = position x of particle n
c part(2,n) = velocity vx of particle n
c part(3,n) = velocity vy of particle n
c part(4,n) = velocity vz of particle n
c fps = distribution function, number of particles in each velocity
c and spatial range
c maximum velocity (used for scaling) is contained in first element fps.
c vdrift for i-th dimension is contained in fpsm(1,i)
c vth for i-th dimension is contained in fpsm(2,i)
c entropy for i-th dimension is contained in fpsm(3,i), defined to be:
c s/k = -sum(f(v)/np)*log(f(v)/(np*delta_v)).  Assumes that
c distribution distributions in each dimension are independent.
c psm = double precision scratch array
c nx = system length in x direction
c nxb = number of spatial regions for distribution
c idimp = size of phase space = 4
c np = number of particles
c the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer nx, nxb, idimp, np, nmv, nmvf
      real part, fps, fpsm
      dimension part(idimp,np), fps(nmvf,3,nxb), fpsm(3,3,nxb)
      double precision psm
      dimension psm(2,3,nxb)
c local data
      double precision sumvx, sumvy, sumvz, sumvx2, sumvy2, sumvz2, anp
      real anmv, svx, svy, svz, at1
      integer j, k, n, nvx, nvy, nvz
      anmv = real(nmv)
c zero out distribution
      do 20 k = 1, nxb
      do 10 j = 2, nmvf
      fps(j,1,k) = 0.0
      fps(j,2,k) = 0.0
      fps(j,3,k) = 0.0
   10 continue
      fpsm(1,1,k) = anmv/fps(1,1,k)
      fpsm(1,2,k) = anmv/fps(1,2,k)
      fpsm(1,3,k) = anmv/fps(1,3,k)
      psm(1,1,k) = 0.0d0
      psm(2,1,k) = 0.0d0
      psm(1,2,k) = 0.0d0
      psm(2,2,k) = 0.0d0
      psm(1,3,k) = 0.0d0
      psm(2,3,k) = 0.0d0
   20 continue
c count particles in each velocity and spatial region
      at1 = real(nxb)/real(nx)
      anmv = anmv + 2.5
      do 30 j = 1, np
      n = part(1,j)*at1
      n = n + 1
      if ((n.ge.1).and.(n.le.nxb)) then
         nvx = part(2,j)*fpsm(1,1,n) + anmv
         psm(1,1,n) = psm(1,1,n) + part(2,j)
         psm(2,1,n) = psm(2,1,n) + part(2,j)**2
         nvy = part(3,j)*fpsm(1,2,n) + anmv
         psm(1,2,n) = psm(1,2,n) + part(3,j)
         psm(2,2,n) = psm(2,2,n) + part(3,j)**2
         nvz = part(4,j)*fpsm(1,3,n) + anmv
         psm(1,3,n) = psm(1,3,n) + part(4,j)
         psm(2,3,n) = psm(2,3,n) + part(4,j)**2
         if ((nvx.ge.2).and.(nvx.le.nmvf)) fps(nvx,1,n) = fps(nvx,1,n) +
     11.0
         if ((nvy.ge.2).and.(nvy.le.nmvf)) fps(nvy,2,n) = fps(nvy,2,n) +
     11.0
         if ((nvz.ge.2).and.(nvz.le.nmvf)) fps(nvz,3,n) = fps(nvz,3,n) +
     11.0
      endif
   30 continue
c calculate velocity moments
      anp = 1.0d0/dble(np)
      do 50 k = 1, nxb
      svx = fpsm(1,1,k)
      svy = fpsm(1,2,k)
      svz = fpsm(1,3,k)
      sumvx = psm(1,1,k)*anp
      fpsm(1,1,k) = sumvx
      fpsm(2,1,k) = dsqrt(psm(2,1,k)*anp - sumvx**2)
      sumvy = psm(1,2,k)*anp
      fpsm(1,2,k) = sumvy
      fpsm(2,2,k) = dsqrt(psm(2,2,k)*anp - sumvy**2)
      sumvz = psm(1,3,k)*anp
      fpsm(1,3,k) = sumvz
      fpsm(2,3,k) = dsqrt(psm(2,3,k)*anp - sumvz**2)
c calculate entropy
      sumvx = 0.0d0
      sumvy = 0.0d0
      sumvz = 0.0d0
      sumvx2 = 0.0d0
      sumvy2 = 0.0d0
      sumvz2 = 0.0d0
      do 40 j = 2, nmvf
      if (fps(j,1,k).gt.0.0) then
         sumvx = sumvx + fps(j,1,k)
         sumvx2 = sumvx2 + fps(j,1,k)*dlog(dble(fps(j,1,k)*svx))
      endif
      if (fps(j,2,k).gt.0.0) then
         sumvy = sumvy + fps(j,2,k)
         sumvy2 = sumvy2 + fps(j,2,k)*dlog(dble(fps(j,2,k)*svy))
      endif
      if (fps(j,3,k).gt.0.0) then
         sumvz = sumvz + fps(j,3,k)
         sumvz2 = sumvz2 + fps(j,3,k)*dlog(dble(fps(j,3,k)*svz))
      endif
   40 continue
      if (sumvx.gt.0.0d0) sumvx = -sumvx2/sumvx + dlog(sumvx)
      if (sumvy.gt.0.0d0) sumvy = -sumvy2/sumvy + dlog(sumvy)
      if (sumvz.gt.0.0d0) sumvz = -sumvz2/sumvz + dlog(sumvz)
      fpsm(3,1,k) = sumvx
      fpsm(3,2,k) = sumvy
      fpsm(3,3,k) = sumvz
   50 continue
      return
      end
      subroutine FWRITE1(f,nx,nxv,iunit,nrec,lrec,name)
c this subroutine writes real 1d data f to a direct access file
c f = input data to be written
c nx = length of data f in x to write
c nxv = first dimension of data array f, must be >= nx
c iunit = fortran unit number
c nrec = current record number for write, if nrec > 0
c if nrec < 0, open new file and write first record
c if nrec = 0, open old file, do not write
c lrec = record length (used only if nrec <= 0)
c name = file name (used only if nrec <= 0)
c input: f, nx, nxv, iunit, nrec, lrec, fname
c output: nrec
      implicit none
      integer nx, nxv, iunit, nrec, lrec
      real f
      character*(*) name
      dimension f(nxv)
c local data
      integer j
c open new file
      if (nrec.lt.0) then
         open(unit=iunit,file=name,form='unformatted',access='direct',re
     1cl=lrec,status='replace')
         nrec = 1
c open old file
      else if (nrec.eq.0) then
         open(unit=iunit,file=name,form='unformatted',access='direct',re
     1cl=lrec,status='old')
         return
      endif
      write (unit=iunit,rec=nrec) (f(j),j=1,nx)
      nrec = nrec + 1
      return
      end
      subroutine FREAD1(f,nx,nxv,iunit,nrec,lrec,name,ierr)
c this subroutine reads real 1d data f from a file
c f = input data to be read
c nx = length of data f in x to read
c nxv = first dimension of data array f, must be >= nx
c iunit = fortran unit number
c nrec = current record number for read, if nrec >  0
c if nrec < 0, open old file and read first record
c if nrec = 0, open old file, do not read
c lrec = record length (used only if nrec <= 0)
c name = file name (used only if nren <= 0)
c input: nx, nxv, iunit, nrec, lrec, fname
c output: f, nrec, ierr
      implicit none
      integer nx, nxv, iunit, nrec, lrec, ierr
      real f
      character*(*) name
      dimension f(nxv)
c local data
      integer j
      ierr = 0
      if (nrec.le.0) then
         open(unit=iunit,file=name,form='unformatted',access='direct',re
     1cl=lrec,status='old')
         if (nrec.eq.0) return
         nrec = 1
      endif
      read (unit=iunit,rec=nrec,err=10) (f(j),j=1,nx)
      nrec = nrec + 1
      return
   10 ierr = 1
      return
      end
      subroutine FCWRITE1(f,nx,nxv,iunit,nrec,lrec,name)
c this subroutine writes complex 1d data f to a direct access file
c f = input data to be written
c nx = length of data f in x to write
c nxv = first dimension of data array f, must be >= nx
c iunit = fortran unit number
c nrec = current record number for write, if nrec > 0
c if nrec < 0, open new file and write first record
c if nrec = 0, open old file, do not write
c lrec = record length (used only if nrec <= 0)
c name = file name (used only if nrec <= 0)
c input: f, nx, nxv, iunit, nrec, lrec, fname
c output: nrec
      implicit none
      integer nx, nxv, iunit, nrec, lrec
      complex f
      character*(*) name
      dimension f(nxv)
c local data
      integer j
c open new file
      if (nrec.lt.0) then
         open(unit=iunit,file=name,form='unformatted',access='direct',re
     1cl=lrec,status='replace')
         nrec = 1
c open old file
      else if (nrec.eq.0) then
         open(unit=iunit,file=name,form='unformatted',access='direct',re
     1cl=lrec,status='old')
         return
      endif
      write (unit=iunit,rec=nrec) (f(j),j=1,nx)
      nrec = nrec + 1
      return
      end
      subroutine FCREAD1(f,nx,nxv,iunit,nrec,lrec,name,ierr)
c this subroutine reads complex 1d data f from a file
c f = input data to be read
c nx = length of data f in x to read
c nxv = first dimension of data array f, must be >= nx
c iunit = fortran unit number
c nrec = current record number for read, if nrec >  0
c if nrec < 0, open old file and read first record
c if nrec = 0, open old file, do not read
c lrec = record length (used only if nrec <= 0)
c name = file name (used only if nren <= 0)
c input: nx, nxv, iunit, nrec, lrec, fname
c output: f, nrec, ierr
      implicit none
      integer nx, nxv, iunit, nrec, lrec, ierr
      complex f
      character*(*) name
      dimension f(nxv)
c local data
      integer j
      ierr = 0
      if (nrec.le.0) then
         open(unit=iunit,file=name,form='unformatted',access='direct',re
     1cl=lrec,status='old')
         if (nrec.eq.0) return
         nrec = 1
      endif
      read (unit=iunit,rec=nrec,err=10) (f(j),j=1,nx)
      nrec = nrec + 1
      return
   10 ierr = 1
      return
      end
      subroutine FWRITE0(f,nxp,iunit,nrec,name)
c this subroutine write real data f to a file
c f = input data to be written
c nxp = size of file f
c iunit = fortran unit number
c nrec = current record number for write (if negative, open file with
c recl=-nren)
c name = file name (used only if nren < 0)
c input: f, nxp, iunit, nrec, fname
c output: nrec
      implicit none
      integer nxp, iunit, nrec
      real f
      character*(*) name
      dimension f(nxp)
c local data
      integer lrec
      if (nrec.lt.0) then
         lrec = -nrec
         open(unit=iunit,file=name,form='unformatted',access='direct',re
     1cl=lrec,status='unknown')
         nrec = 1
      endif
      write (unit=iunit,rec=nrec) f
      nrec = nrec + 1
      return
      end
