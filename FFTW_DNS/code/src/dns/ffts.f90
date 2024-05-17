!
!=======================================================================
!
      subroutine saxpy(n,sa,sx,incx,sy,incy)

!
!     constant times a vector plus a vector.
!     uses unrolled loop for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!
      real sx(1),sy(1),sa
      integer i,incx,incy,ix,iy,m,mp1,n
!
      if(n.le.0)return
      if (sa .eq. 0.0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sy(iy) + sa*sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sy(i) + sa*sx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        sy(i) = sy(i) + sa*sx(i)
        sy(i + 1) = sy(i + 1) + sa*sx(i + 1)
        sy(i + 2) = sy(i + 2) + sa*sx(i + 2)
        sy(i + 3) = sy(i + 3) + sa*sx(i + 3)
   50 continue
      return
      end subroutine
!
!=======================================================================
!
      subroutine cftfax(n,ifax,trigs)

      dimension ifax(13),trigs(1)
!
!     this routine was modified from temperton's original
!     by dave fulker.  it no longer produces factors in ascending
!     order, and there are none of the original 'mode' options.
!
! on input     n
!               the length of each complex transform to be performed
!
!               n must be greater than 1 and contain no prime
!               factors greater than 5.
!
! on output    ifax
!               ifax(1)
!                 the number of factors chosen or -99 in case of error
!               ifax(2) thru ifax( ifax(1)+1 )
!                 the factors of n in the followin order:  appearing
!                 first are as many factors of 4 as can be obtained.
!                 subsequent factors are primes, and appear in
!                 ascending order, except for multiple factors.
!
!              trigs
!               2n sin and cos values for use by the transform routine
!
      call fact(n,ifax)
      k = ifax(1)
      if (k .lt. 1 .or. ifax(k+1) .gt. 5) ifax(1) = -99
      if (ifax(1) .le. 0 )then
        write(*,1900)n
1900  format(' fftfax - invalid n=',i20)
        return
        endif
      call cftrig (n, trigs)
      return
      end
!
!=======================================================================
!
      subroutine cftrig(n,trigs)

      dimension trigs(1)
      pi=2.0*asin(1.0)
      del=(pi+pi)/float(n)
      l=n+n
      do 10 i=1,l,2
      angle=0.5*float(i-1)*del
      trigs(i)=cos(angle)
      trigs(i+1)=sin(angle)
   10 continue
      return
      end
!
!=======================================================================
!
      subroutine rfftmlt(a,work,trigs,ifax,inc,jump,n,lot,isgn)

!
!     "rfftmlt" - multiple real/half-complex periodic
!     fast fourier transform
!
!     procedure used to convert to half-length complex transform
!     is given by cooley, lewis and welch (j. sound vib., vol. 12
!     (1970), 315-337)
!
!     a is the array containing input and output data
!     work is an area of size (n+1)*lot
!     trigs is a previously prepared list of trig function values
!     ifax is a previously prepared list of factors of n/2
!     inc is the increment within each data 'vector'
!         (e.g. inc=1 for consecutively stored data)
!     jump is the increment between the start of each data vector
!     n is the length of the data vectors
!     lot is the number of data vectors
!     isgn = +1 for transform from spectral to gridpoint
!           = -1 for transform from gridpoint to spectral
!
!     ordering of coefficients:
!         a(0),b(0),a(1),b(1),a(2),b(2),...,a(n/2),b(n/2)
!         where b(0)=b(n/2)=0; (n+2) locations required
!
!     ordering of data:
!         x(0),x(1),x(2),...,x(n-1)
!
!     vectorization is achieved on cray by doing the transforms in
!     parallel
!
!     *** n.b. n is assumed to be an even number
!
!     definition of transforms:
!     -------------------------
!
!     isgn=+1: x(j)=sum(k=0,...,n-1)(c(k)*exp(2*i*j*k*pi/n))
!         where c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)
!
!     isgn=-1: a(k)=(1/n)*sum(j=0,...,n-1)(x(j)*cos(2*j*k*pi/n))
!               b(k)=-(1/n)*sum(j=0,...,n-1)(x(j)*sin(2*j*k*pi/n))
!
      dimension a(n),work(n),trigs(n),ifax(1)
      nfax=ifax(1)
      nx=n+1
      nh=n/2
      ink=inc+inc
      if (isgn.eq.+1) go to 30
!
!     if necessary, transfer data to work area
      igo=50
      if (mod(nfax,2).eq.1) goto 40
      ibase=1
      jbase=1
      do 20 l=1,lot
      i=ibase
      j=jbase
! 
      do 10 m=1,n
      work(j)=a(i)
      i=i+inc
      j=j+1
   10 continue
      ibase=ibase+jump
      jbase=jbase+nx
   20 continue
!
      igo=60
      go to 40
!
!     preprocessing (isgn=+1)
!     ------------------------
!
   30 continue
      call fft99a(a,work,trigs,inc,jump,n,lot)
      igo=60
!
!     complex transform
!     -----------------
!
   40 continue
      ia=1
      la=1
      do 80 k=1,nfax
      if (igo.eq.60) go to 60
   50 continue
      call vpassm(a(ia),a(ia+inc),work(1),work(2),trigs,ink,2,jump,nx,lot,nh,ifax(k+1),la)
      igo=60
      go to 70
   60 continue
      call vpassm(work(1),work(2),a(ia),a(ia+inc),trigs,2,ink,nx,jump,lot,nh,ifax(k+1),la)
      igo=50
   70 continue
      la=la*ifax(k+1)
   80 continue
!
      if (isgn.eq.-1) go to 130
!
!     if necessary, transfer data from work area
      if (mod(nfax,2).eq.1) go to 110
      ibase=1
      jbase=1
      do 100 l=1,lot
      i=ibase
      j=jbase
! 
      do 90 m=1,n
      a(j)=work(i)
      i=i+1
      j=j+inc
   90 continue
      ibase=ibase+nx
      jbase=jbase+jump
  100 continue
!
!     fill in zeros at end
  110 continue
      ib=n*inc+1
! 
      do 120 l=1,lot
      a(ib)=0.0
      a(ib+inc)=0.0
      ib=ib+jump
  120 continue
      go to 140
!
!     postprocessing (isgn=-1):
!     --------------------------
!
  130 continue
      call fft99b(work,a,trigs,inc,jump,n,lot)
!
  140 continue
      return
      end
!
!=======================================================================
!
      subroutine fftfax(n,ifax,trigs)

      dimension ifax(13),trigs(1)
!
! mode 3 is used for real/half-complex transforms.  it is possible
! to do complex/complex transforms with other values of mode, but
! documentation of the details were not available when this routine
! was written.
!
      data mode /3/

      call fax (ifax, n, mode)
      i = ifax(1)
      if (ifax(i+1) .gt. 5 .or. n .le. 4) ifax(1) = -99
      if (ifax(1) .le. 0 )then
        write(*,1900)n
1900    format(' fftfax - invalid n=',i20)
        return
        endif
      call fftrig (trigs, n, mode)
      return
      end
!
!=======================================================================
!
      subroutine fftrig(trigs,n,mode)

      dimension trigs(1)
      pi=2.0*asin(1.0)
      imode=iabs(mode)
      nn=n
      if (imode.gt.1.and.imode.lt.6) nn=n/2
      del=(pi+pi)/float(nn)
      l=nn+nn
      do 10 i=1,l,2
      angle=0.5*float(i-1)*del
      trigs(i)=cos(angle)
      trigs(i+1)=sin(angle)
   10 continue
      if (imode.eq.1) return
      if (imode.eq.8) return
      del=0.5*del
      nh=(nn+1)/2
      l=nh+nh
      la=nn+nn
      do 20 i=1,l,2
      angle=0.5*float(i-1)*del
      trigs(la+i)=cos(angle)
      trigs(la+i+1)=sin(angle)
   20 continue
      if (imode.le.3) return
      del=0.5*del
      la=la+nn
      if (mode.eq.5) go to 40
      do 30 i=2,nn
      angle=float(i-1)*del
      trigs(la+i)=2.0*sin(angle)
   30 continue
      return
   40 continue
      del=0.5*del
      do 50 i=2,n
      angle=float(i-1)*del
      trigs(la+i)=sin(angle)
   50 continue
      return
      end
!
!=======================================================================
!
      subroutine cfftmlt(ar,ai,work,trigs,ifax,inc,jump,n,lot,isgn)

!
! purpose      performs multiple fast fourier transforms.  this package
!              will perform a number of simultaneous complex periodic
!              fourier transforms or corresponding inverse transforms.
!              that is, given a set of complex gridpoint vectors, the
!              package returns a set of complex fourier
!              coefficient vectors, or vice versa.  the length of the
!              transforms must be a number greater than 1 that has
!              no prime factors other than 2, 3, and 5.
!
!              the package cfft99 contains several user-level routines:
!
!            subroutine cftfax
!                an initialization routine that must be called once
!                before a sequence of calls to cfft99
!                (provided that n is not changed).
!
!            subroutine cfft99
!                the actual transform routine routine, cabable of
!                performing both the transform and its inverse.
!                however, as the transforms are not normalized,
!                the application of a transform followed by its
!                inverse will yield the original values multiplied
!                by n.
!
!
! access       *fortran,p=xlib,sn=cfft99
!
!
! usage        let n be of the form 2**p * 3**q * 5**r, where p .ge. 0,
!              q .ge. 0, and r .ge. 0.  then a typical sequence of
!              calls to transform a given set of complex vectors of
!              length n to a set of (unscaled) complex fourier
!              coefficient vectors of length n is
!
!                   dimension ifax(13),trigs(2*n)
!                   complex a(...), work(...)
!
!                   call cftfax (n, ifax, trigs)
!                   call cfft99 (a,work,trigs,ifax,inc,jump,n,lot,isgn)
!
!              the output vectors overwrite the input vectors, and
!              these are stored in a.  with appropriate choices for
!              the other arguments, these vectors may be considered
!              either the rows or the columns of the array a.
!              see the individual write-ups for cftfax and
!              cfft99 below, for a detailed description of the
!              arguments.
!
! history      the package was written by clive temperton at ecmwf in
!              november, 1978.  it was modified, documented, and tested
!              for ncar by russ rew in september, 1980.  it was
!              further modified for the fully complex case by dave
!              fulker in november, 1980.
!
!-----------------------------------------------------------------------
!
! subroutine cftfax (n,ifax,trigs)
!
! purpose      a set-up routine for cfft99.  it need only be
!              called once before a sequence of calls to cfft99,
!              provided that n is not changed.
!
! argument     ifax(13),trigs(2*n)
! dimensions
!
! arguments
!
! on input     n
!               an even number greater than 1 that has no prime factor
!               greater than 5.  n is the length of the transforms (see
!               the documentation for cfft99 for the definition of
!               the transforms).
!
!              ifax
!               an integer array.  the number of elements actually used
!               will depend on the factorization of n.  dimensioning
!               ifax for 13 suffices for all n less than 1 million.
!
!              trigs
!               a real array of dimension 2*n
!
! on output    ifax
!               contains the factorization of n.  ifax(1) is the
!               number of factors, and the factors themselves are stored
!               in ifax(2),ifax(3),...  if n has any prime factors
!               greater than 5, ifax(1) is set to -99.
!
!              trigs
!               an array of trigonometri! function values subsequently
!               used by the cft routines.
!
!-----------------------------------------------------------------------
!
! subroutine cfft99 (a,work,trigs,ifax,inc,jump,n,lot,isgn)
!
! purpose      perform a number of simultaneous (unnormalized) complex
!              periodic fourier transforms or corresponding inverse
!              transforms.  given a set of complex gridpoint
!              vectors, the package returns a set of
!              complex fourier coefficient vectors, or vice
!              versa.  the length of the transforms must be a
!              number having no prime factors other than
!              2, 3, and 5.  this routine is
!              optimized for use on the cray-1.
!
! argument     complex a(n*inc+(lot-1)*jump), work(n*lot)
! dimensions   real trigs(2*n), integer ifax(13)
!
! arguments
!
! on input     a
!               a complex array of length n*inc+(lot-1)*jump containing
!               the input gridpoint or coefficient vectors.  this array
!               overwritten by the results.
!
!               n.b. although the array a is usually considered to be of
!               type complex in the calling program, it is treated as
!               real within the transform package.  this requires that
!               such type conflicts are permitted in the user"s
!               environment, and that the storage of complex numbers
!               matches the assumptions of this routine.  this routine
!               assumes that the real and imaginary portions of a
!               complex number occupy adjacent elements of memory.  if
!               these conditions are not met, the user must treat the
!               array a as real (and of twice the above length), and
!               write the calling program to treat the real and
!               imaginary portions explicitly.
!
!              work
!               a complex work array of length n*lot or a real array
!               of length 2*n*lot.  see n.b. above.
!
!              trigs
!               an array set up by cftfax, which must be called first.
!
!              ifax
!               an array set up by cftfax, which must be called first.
!
!
!               n.b. in the following arguments, increments are measured
!               in word pairs, because each complex element is assumed
!               to occupy an adjacent pair of words in memory.
!
!              inc
!               the increment (in word pairs) between successive element
!               of each (complex) gridpoint or coefficient vector
!               (e.g.  inc=1 for consecutively stored data).
!
!              jump
!               the increment (in word pairs) between the first elements
!               of successive data or coefficient vectors.  on the cray-
!               try to arrange data so that jump is not a multiple of 8
!               (to avoid memory bank conflicts).  for clarification of
!               inc and jump, see the examples below.
!
!              n
!               the length of each transform (see definition of
!               transforms, below).
!
!              lot
!               the number of transforms to be done simultaneously.
!
!              isgn
!               = -1 for a transform from gridpoint values to fourier
!                    coefficients.
!               = +1 for a transform from fourier coefficients to
!                    gridpoint values.
!
! on output    a
!               if isgn = -1, and lot gridpoint vectors are supplied,
!               each containing the complex sequence:
!
!               g(0),g(1), ... ,g(n-1)  (n complex values)
!
!               then the result consists of lot complex vectors each
!               containing the corresponding n coefficient values:
!
!               c(0),c(1), ... ,c(n-1)  (n complex values)
!
!               defined by:
!                 c(k) = sum(j=0,...,n-1)( g(j)*exp(-2*i*j*k*pi/n) )
!                 where i = sqrt(-1)
!
!
!               if isgn = +1, and lot coefficient vectors are supplied,
!               each containing the complex sequence:
!
!               c(0),c(1), ... ,c(n-1)  (n complex values)
!
!               then the result consists of lot complex vectors each
!               containing the corresponding n gridpoint values:
!
!               g(0),g(1), ... ,g(n-1)  (n complex values)
!
!               defined by:
!                 g(j) = sum(k=0,...,n-1)( g(k)*exp(+2*i*j*k*pi/n) )
!                 where i = sqrt(-1)
!
!
!               a call with isgn=-1 followed by a call with isgn=+1
!               (or vice versa) returns the original data, multiplied
!               by the factor n.
!
!
! example       given a 64 by 9 grid of complex values, stored in
!               a 66 by 9 complex array, a, compute the two dimensional
!               fourier transform of the grid.  from transform theory,
!               it is known that a two dimensional transform can be
!               obtained by first transforming the grid along one
!               direction, then transforming these results along the
!               orthogonal direction.
!
!               complex a(66,9), work(64,9)
!               real trigs1(128), trigs2(18)
!               integer ifax1(13), ifax2(13)
!
!               set up the ifax and trigs arrays for each direction:
!
!               call cftfax(64, ifax1, trigs1)
!               call cftfax( 9, ifax2, trigs2)
!
!               in this case, the complex values of the grid are
!               stored in memory as follows (using u and v to
!               denote the real and imaginary components, and
!               assuming conventional fortran storage):
!
!   u(1,1), v(1,1), u(2,1), v(2,1),  ...  u(64,1), v(64,1), 4 nulls,
!
!   u(1,2), v(1,2), u(2,2), v(2,2),  ...  u(64,2), v(64,2), 4 nulls,
!
!   .       .       .       .         .   .        .        .
!   .       .       .       .         .   .        .        .
!   .       .       .       .         .   .        .        .
!
!   u(1,9), v(1,9), u(2,9), v(2,9),  ...  u(64,9), v(64,9), 4 nulls.
!
!               we choose (arbitrarily) to transorm first along the
!               direction of the first subscript.  thus we define
!               the length of the transforms, n, to be 64, the
!               number of transforms, lot, to be 9, the increment
!               between elements of each transform, inc, to be 1,
!               and the increment between the starting points
!               for each transform, jump, to be 66 (the first
!               dimension of a).
!
!               call cfft99( a, work, trigs1, ifax1, 1, 66, 64, 9, -1)
!
!               to transform along the direction of the second subscript
!               the roles of the increments are reversed.  thus we defin
!               the length of the transforms, n, to be 9, the
!               number of transforms, lot, to be 64, the increment
!               between elements of each transform, inc, to be 66,
!               and the increment between the starting points
!               for each transform, jump, to be 1
!
!               call cfft99( a, work, trigs2, ifax2, 66, 1, 9, 64, -1)
!
!               these two sequential steps results in the two-dimensiona
!               fourier coefficient array overwriting the input
!               gridpoint array, a.  the same two steps applied again
!               with isgn = +1 would result in the reconstruction of
!               the gridpoint array (multiplied by a factor of 64*9).
!
!
!-----------------------------------------------------------------------
      dimension ar(n),ai(n),work(n),trigs(n),ifax(n)
!
!     subroutine "cfft99" - multiple fast complex fourier transform
!
!     a is the array containing input and output data
!     work is an area of size n*lot
!     trigs is a previously prepared list of trig function values
!     ifax is a previously prepared list of factors of n
!     inc is the increment within each data 'vector'
!         (e.g. inc=1 for consecutively stored data)
!     jump is the increment between the start of each data vector
!     n is the length of the data vectors
!     lot is the number of data vectors
!     isgn = +1 for transform from spectral to gridpoint
!           = -1 for transform from gridpoint to spectral
!
!
!     vectorization is achieved on cray by doing the transforms in
!     parallel.
!
      nn = n+n
      ink=inc+inc
      jum = jump+jump
      nfax=ifax(1)
      jnk = 2
      jst = 2
      if (isgn.ge.0) go to 30
!
!     the innermost temperton routines have no facility for the
!     forward (isgn = -1) transform.  therefore, the input must be
!     rearranged as follows:
!
!     the order of each input vector,
!
!     g(0), g(1), g(2), ... , g(n-2), g(n-1)
!
!     is reversed (excluding g(0)) to yield
!
!     g(0), g(n-1), g(n-2), ... , g(2), g(1).
!
!     within the transform, the corresponding exponential multiplier
!     is then precisely the conjugate of that for the normal
!     ordering.  thus the forward (isgn = -1) transform is
!     accomplished
!
!     for nfax odd, the input must be transferred to the work array,
!     and the rearrangement can be done during the move.
!
      jnk = -2
      jst = nn-2
      if (mod(nfax,2).eq.1) goto 40
!
!     for nfax even, the rearrangement must be applied directly to
!     the input array.  this can be done by swapping elements.
!
      ibase = 1
      ilast = (n-1)*inc
      nh = n/2
      do 20 l=1,lot
      i1 = ibase+inc
      i2 = ibase+ilast
! 
      do 10 m=1,nh
!     swap real and imaginary portions
      hreal = ar(i1)
      himag = ai(i1)
      ar(i1) = ar(i2)
      ai(i1) = ai(i2)
      ar(i2) = hreal
      ai(i2) = himag
      i1 = i1+inc
      i2 = i2-inc
   10 continue
      ibase = ibase+jump
   20 continue
      goto 100
!
   30 continue
      if (mod(nfax,2).eq.0) goto 100
!
   40 continue
!
!     during the transform process, nfax steps are taken, and the
!     results are stored alternately in work and in a.  if nfax is
!     odd, the input data are first moved to work so that the final
!     result (after nfax steps) is stored in array a.
!
      ibase=1
      jbase=1
      do 60 l=1,lot
!     move real and imaginary portions of element zero
      work(jbase) = ar(ibase)
      work(jbase+1) = ai(ibase)
      i=ibase+inc
      j=jbase+jst
! 
      do 50 m=2,n
!     move real and imaginary portions of other elements (possibly in
!     reverse order, depending on jst and jnk)
      work(j) = ar(i)
      work(j+1) = ai(i)
      i=i+inc
      j=j+jnk
   50 continue
      ibase=ibase+jump
      jbase=jbase+nn
   60 continue
!
  100 continue
!
!     perform the transform passes, one pass for each factor.  during
!     each pass the data are moved from a to work or from work to a.
!
!     for nfax even, the first pass moves from a to work
      igo = 110
!     for nfax odd, the first pass moves from work to a
      if (mod(nfax,2).eq.1) igo = 120
      la=1
      do 140 k=1,nfax
      if (igo.eq.120) go to 120
  110 continue
      call vpassm(ar,ai,work(1),work(2),trigs,inc,2,jump,nn,lot,n,ifax(k+1),la)
      igo=120
      go to 130
  120 continue
      call vpassm(work(1),work(2),ar,ai,trigs,2,inc,nn,jump,lot,n,ifax(k+1),la)
      igo=110
  130 continue
      la=la*ifax(k+1)
  140 continue
!
!     at this point the final transform result is stored in a.
!
      return
      end
!
!=======================================================================
!
      subroutine  scopy(n,sx,incx,sy,incy)

!
!     copies a vector, x, to a vector, y.
!     uses unrolled loops for increments equal to 1.
!     jack dongarra, linpack, 3/11/78.
!
      real sx(1),sy(1)
      integer i,incx,incy,ix,iy,m,mp1,n
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        sy(i) = sx(i)
        sy(i + 1) = sx(i + 1)
        sy(i + 2) = sx(i + 2)
        sy(i + 3) = sx(i + 3)
        sy(i + 4) = sx(i + 4)
        sy(i + 5) = sx(i + 5)
        sy(i + 6) = sx(i + 6)
   50 continue
      return
      end
!
!=======================================================================
!
      real function ssum(n,sx,incx)
!
!        takes the sum of the values of a vector.
!        uses unrolled loops for increment equal to one.
!        jack dongarra, linpack, 3/11/78.
!
      real sx(1),stemp
      integer i,incx,m,mp1,n,ix
!
      ssum = 0.0e0
      stemp = 0.0e0
      if (n .lt. 0) stop
      if (n .eq. 0) return
      if (incx .eq. 1) go to 20
!
!        code for increment not equal to 1
!
      ix = 1
      if (incx .le. 0) ix = (-n+1)*incx + 1
      do 10 i = 1,n
        stemp = stemp + sx(ix)
        ix = ix + incx
   10 continue
      ssum = stemp
      return
!
!        code for increment equal to 1
!
!        clean-up loop
!
   20 m = mod(n,6)
      if (m .eq. 0) go to 40
      do 30 i = 1,m
        stemp = stemp + sx(i)
   30 continue
      if (n .lt. 6) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        stemp = stemp + sx(i) + sx(i+1) + sx(i+2) + sx(i+3) + sx(i+4) + sx(i+5)
   50 continue
   60 ssum = stemp
      return
      end
!
!=======================================================================
!
      subroutine ccopy(n,cx,incx,cy,incy)

!
!     copies a vector, x, to a vector, y.
!     jack dongarra, linpack, 3/11/78.
!
      complex cx(1),cy(1)
      integer i,incx,incy,ix,iy,n
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        cy(iy) = cx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        code for both increments equal to 1
!
   20 do 30 i = 1,n
        cy(i) = cx(i)
   30 continue
      return
      end
!
!=======================================================================
!
      real function snrm2 ( n, sx, incx)
      integer          next
      real   sx(1),  cutlo, cuthi, hitest, sum, xmax, zero, one
      data   zero, one /0.0e0, 1.0e0/
!
!     euclidean norm of the n-vector stored in sx() with storage
!     increment incx .
!     if    n .le. 0 return with result = 0.
!     if n .ge. 1 then incx must be .ge. 1
!
!           c.l.lawson, 1978 jan 08
!
!     four phase method     using two built-in constants that are
!     hopefully applicable to all machines.
!         cutlo = maximum of  sqrt(u/eps)  over all known machines.
!         cuthi = minimum of  sqrt(v)      over all known machines.
!     where
!         eps = smallest no. such that eps + 1. .gt. 1.
!         u   = smallest positive no.   (underflow limit)
!         v   = largest  no.            (overflow  limit)
!
!     brief outline of algorithm..
!
!     phase 1    scans zero components.
!     move to phase 2 when a component is nonzero and .le. cutlo
!     move to phase 3 when a component is .gt. cutlo
!     move to phase 4 when a component is .ge. cuthi/m
!     where m = n for x() real and m = 2*n for complex.
!
!     values for cutlo and cuthi..
!     from the environmental parameters listed in the imsl converter
!     document the limiting values are as follows..
!     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
!                   univac and dec at 2**(-103)
!                   thus cutlo = 2**(-51) = 4.44089e-16
!     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
!                   thus cuthi = 2**(63.5) = 1.30438e19
!     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
!                   thus cutlo = 2**(-33.5) = 8.23181d-11
!     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
!     data cutlo, cuthi / 8.232d-11,  1.304d19 /
!     data cutlo, cuthi / 4.441e-16,  1.304e19 /
      data cutlo, cuthi / 4.441e-16,  1.304e19 /
!
      if(n .gt. 0) go to 10
         snrm2  = zero
         go to 300
!
   10 assign 30 to next
      sum = zero
      nn = n * incx
!                                                 begin main loop
      i = 1
   20    go to next,(30, 50, 70, 110)
   30 if( abs(sx(i)) .gt. cutlo) go to 85
      assign 50 to next
      xmax = zero
!
!                        phase 1.  sum is zero
!
   50 if( sx(i) .eq. zero) go to 200
      if( abs(sx(i)) .gt. cutlo) go to 85
!
!                                prepare for phase 2.
      assign 70 to next
      go to 105
!
!                                prepare for phase 4.
!
  100 i = j
      assign 110 to next
      sum = (sum / sx(i)) / sx(i)
  105 xmax = abs(sx(i))
      go to 115
!
!                   phase 2.  sum is small.
!                             scale to avoid destructive underflow.
!
   70 if( abs(sx(i)) .gt. cutlo ) go to 75
!
!                     common code for phases 2 and 4.
!                     in phase 4 sum is large.  scale to avoid overflow.
!
  110 if( abs(sx(i)) .le. xmax ) go to 115
         sum = one + sum * (xmax / sx(i))**2
         xmax = abs(sx(i))
         go to 200
!
  115 sum = sum + (sx(i)/xmax)**2
      go to 200
!
!
!                  prepare for phase 3.
!
   75 sum = (sum * xmax) * xmax
!
!
!     for real or d.p. set hitest = cuthi/n
!     for complex      set hitest = cuthi/(2*n)
!
   85 hitest = cuthi/float( n )
!
!                   phase 3.  sum is mid-range.  no scaling.
!
      do 95 j =i,nn,incx
      if(abs(sx(j)) .ge. hitest) go to 100
   95    sum = sum + sx(j)**2
      snrm2 = sqrt( sum )
      go to 300
!
  200 continue
      i = i + incx
      if ( i .le. nn ) go to 20
!
!              end of main loop.
!
!              compute square root and adjust for scaling.
!
      snrm2 = xmax * sqrt(sum)
  300 continue
      return
      end
!
!=======================================================================
!
      subroutine fact(n,ifax)

!     factorization routine that first extracts all factors of 4
      dimension ifax(13)
      if (n.gt.1) go to 10
      ifax(1) = 0
      if (n.lt.1) ifax(1) = -99
      return
   10 nn=n
      k=1
!     test for factors of 4
   20 if (mod(nn,4).ne.0) go to 30
      k=k+1
      ifax(k)=4
      nn=nn/4
      if (nn.eq.1) go to 80
      go to 20
!     test for extra factor of 2
   30 if (mod(nn,2).ne.0) go to 40
      k=k+1
      ifax(k)=2
      nn=nn/2
      if (nn.eq.1) go to 80
!     test for factors of 3
   40 if (mod(nn,3).ne.0) go to 50
      k=k+1
      ifax(k)=3
      nn=nn/3
      if (nn.eq.1) go to 80
      go to 40
!     now find remaining factors
   50 l=5
      max = sqrt(float(nn))
      inc=2
!     inc alternately takes on values 2 and 4
   60 if (mod(nn,l).ne.0) go to 70
      k=k+1
      ifax(k)=l
      nn=nn/l
      if (nn.eq.1) go to 80
      go to 60
   70 if (l.gt.max) go to 75
      l=l+inc
      inc=6-inc
      go to 60
   75 k = k+1
      ifax(k) = nn
   80 ifax(1)=k-1
!     ifax(1) now contains number of factors
      return
      end
!
!=======================================================================
!
      subroutine fax(ifax,n,mode)

      dimension ifax(10)
      nn=n
      if (iabs(mode).eq.1) go to 10
      if (iabs(mode).eq.8) go to 10
      nn=n/2
      if ((nn+nn).eq.n) go to 10
      ifax(1)=-99
      return
   10 k=1
!     test for factors of 4
   20 if (mod(nn,4).ne.0) go to 30
      k=k+1
      ifax(k)=4
      nn=nn/4
      if (nn.eq.1) go to 80
      go to 20
!     test for extra factor of 2
   30 if (mod(nn,2).ne.0) go to 40
      k=k+1
      ifax(k)=2
      nn=nn/2
      if (nn.eq.1) go to 80
!     test for factors of 3
   40 if (mod(nn,3).ne.0) go to 50
      k=k+1
      ifax(k)=3
      nn=nn/3
      if (nn.eq.1) go to 80
      go to 40
!     now find remaining factors
   50 l=5
      inc=2
!     inc alternately takes on values 2 and 4
   60 if (mod(nn,l).ne.0) go to 70
      k=k+1
      ifax(k)=l
      nn=nn/l
      if (nn.eq.1) go to 80
      go to 60
   70 l=l+inc
      inc=6-inc
      go to 60
   80 ifax(1)=k-1
!     ifax(1) contains number of factors
      nfax=ifax(1)
!     sort factors into ascending order
      if (nfax.eq.1) go to 110
      do 100 ii=2,nfax
      istop=nfax+2-ii
      do 90 i=2,istop
      if (ifax(i+1).ge.ifax(i)) go to 90
      item=ifax(i)
      ifax(i)=ifax(i+1)
      ifax(i+1)=item
   90 continue
  100 continue
  110 continue
      return
      end
!
!=======================================================================
!
      subroutine fft99a(a,work,trigs,inc,jump,n,lot)

      dimension a(n),work(n),trigs(n)
!
!     fft99a - preprocessing step, isgn=+1
!     (spectral to gridpoint transform)
!
      nh=n/2
      nx=n+1
      ink=inc+inc
!
!     a(0) and a(n/2)
      ia=1
      ib=n*inc+1
      ja=1
      jb=2
! 
      do 10 l=1,lot
      work(ja)=a(ia)+a(ib)
      work(jb)=a(ia)-a(ib)
      ia=ia+jump
      ib=ib+jump
      ja=ja+nx
      jb=jb+nx
   10 continue
!
!     remaining wavenumbers
      iabase=2*inc+1
      ibbase=(n-2)*inc+1
      jabase=3
      jbbase=n-1
!
      do 30 k=3,nh,2
      ia=iabase
      ib=ibbase
      ja=jabase
      jb=jbbase
      c=trigs(n+k)
      s=trigs(n+k+1)
! 
      do 20 l=1,lot
      work(ja)=(a(ia)+a(ib))-(s*(a(ia)-a(ib))+c*(a(ia+inc)+a(ib+inc)))
      work(jb)=(a(ia)+a(ib))+(s*(a(ia)-a(ib))+c*(a(ia+inc)+a(ib+inc)))
      work(ja+1)=(c*(a(ia)-a(ib))-s*(a(ia+inc)+a(ib+inc)))+(a(ia+inc)-a(ib+inc))
      work(jb+1)=(c*(a(ia)-a(ib))-s*(a(ia+inc)+a(ib+inc)))-(a(ia+inc)-a(ib+inc))
      ia=ia+jump
      ib=ib+jump
      ja=ja+nx
      jb=jb+nx
   20 continue
      iabase=iabase+ink
      ibbase=ibbase-ink
      jabase=jabase+2
      jbbase=jbbase-2
   30 continue
!
      if (iabase.ne.ibbase) go to 50
!     wavenumber n/4 (if it exists)
      ia=iabase
      ja=jabase
! 
      do 40 l=1,lot
      work(ja)=2.0*a(ia)
      work(ja+1)=-2.0*a(ia+inc)
      ia=ia+jump
      ja=ja+nx
   40 continue
!
   50 continue
      return
      end
!
!=======================================================================
!
      subroutine fft99b(work,a,trigs,inc,jump,n,lot)

      dimension work(n),a(n),trigs(n)
!
!     fft99b - postprocessing step, isgn=-1
!     (gridpoint to spectral transform)
!
      nh=n/2
      nx=n+1
      ink=inc+inc
!
!     a(0) and a(n/2)
      scale=1.0/float(n)
      ia=1
      ib=2
      ja=1
      jb=n*inc+1
! 
      do 10 l=1,lot
      a(ja)=scale*(work(ia)+work(ib))
      a(jb)=scale*(work(ia)-work(ib))
      a(ja+inc)=0.0
      a(jb+inc)=0.0
      ia=ia+nx
      ib=ib+nx
      ja=ja+jump
      jb=jb+jump
   10 continue
!
!     remaining wavenumbers
      scale=0.5*scale
      iabase=3
      ibbase=n-1
      jabase=2*inc+1
      jbbase=(n-2)*inc+1
!
      do 30 k=3,nh,2
      ia=iabase
      ib=ibbase
      ja=jabase
      jb=jbbase
      c=trigs(n+k)
      s=trigs(n+k+1)
! 
      do 20 l=1,lot
      a(ja)=scale*((work(ia)+work(ib))+(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
      a(jb)=scale*((work(ia)+work(ib))-(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
      a(ja+inc)=scale*((c*(work(ia)-work(ib))-s*(work(ia+1)+work(ib+1)))+(work(ib+1)-work(ia+1)))
      a(jb+inc)=scale*((c*(work(ia)-work(ib))-s*(work(ia+1)+work(ib+1)))-(work(ib+1)-work(ia+1)))
      ia=ia+nx
      ib=ib+nx
      ja=ja+jump
      jb=jb+jump
   20 continue
      iabase=iabase+2
      ibbase=ibbase-2
      jabase=jabase+ink
      jbbase=jbbase-ink
   30 continue
!
      if (iabase.ne.ibbase) go to 50
!     wavenumber n/4 (if it exists)
      ia=iabase
      ja=jabase
      scale=2.0*scale
! 
      do 40 l=1,lot
      a(ja)=scale*work(ia)
      a(ja+inc)=-scale*work(ia+1)
      ia=ia+nx
      ja=ja+jump
   40 continue
!
   50 continue
      return
      end
!
!=======================================================================
!
      subroutine vpassm(a,b,c,d,trigs,inc1,inc2,inc3,inc4,lot,n,ifac,la)

      dimension a(n),b(n),c(n),d(n),trigs(n)
!
!     "vpassm" - multiple version of "vpassa"
!     performs one pass through data
!     as part of multiple complex fft routine
!     a is first real input vector
!     b is first imaginary input vector
!     c is first real output vector
!     d is first imaginary output vector
!     trigs is precalculated table of sines " cosines
!     inc1 is addressing increment for a and b
!     inc2 is addressing increment for c and d
!     inc3 is addressing increment between a"s & b"s
!     inc4 is addressing increment between c"s & d"s
!     lot is the number of vectors
!     n is length of vectors
!     ifac is current factor of n
!     la is product of previous factors
!
      data sin36/0.587785252292473/,cos36/0.809016994374947/,sin72/0.951056516295154/,cos72/0.309016994374947/,   &
           sin60/0.866025403784437/
!
      m=n/ifac
      iink=m*inc1
      jink=la*inc2
      jump=(ifac-1)*jink
      ibase=0
      jbase=0
      igo=ifac-1
      if (igo.gt.4) return
      go to (10,50,90,130),igo
!
!     coding for factor 2
!
   10 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      do 20 l=1,la
      i=ibase
      j=jbase
! 
      do 15 ijk=1,lot
      c(ja+j)=a(ia+i)+a(ib+i)
      d(ja+j)=b(ia+i)+b(ib+i)
      c(jb+j)=a(ia+i)-a(ib+i)
      d(jb+j)=b(ia+i)-b(ib+i)
      i=i+inc3
      j=j+inc4
   15 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   20 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 40 k=la1,m,la
      kb=k+k-2
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      do 30 l=1,la
      i=ibase
      j=jbase
! 
      do 25 ijk=1,lot
      c(ja+j)=a(ia+i)+a(ib+i)
      d(ja+j)=b(ia+i)+b(ib+i)
      c(jb+j)=c1*(a(ia+i)-a(ib+i))-s1*(b(ia+i)-b(ib+i))
      d(jb+j)=s1*(a(ia+i)-a(ib+i))+c1*(b(ia+i)-b(ib+i))
      i=i+inc3
      j=j+inc4
   25 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   30 continue
      jbase=jbase+jump
   40 continue
      return
!
!     coding for factor 3
!
   50 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      do 60 l=1,la
      i=ibase
      j=jbase
! 
      do 55 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
      c(jb+j)=(a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))
      c(jc+j)=(a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))
      d(jb+j)=(b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i)))
      d(jc+j)=(b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i)))
      i=i+inc3
      j=j+inc4
   55 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   60 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 80 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      do 70 l=1,la
      i=ibase
      j=jbase
! 
      do 65 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
      c(jb+j)=c1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i))))   &
         -s1*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
      d(jb+j)=s1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i))))   &
         +c1*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
      c(jc+j)=c2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i))))   &
         -s2*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i))))
      d(jc+j)=s2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i))))   &
         +c2*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i))))
      i=i+inc3
      j=j+inc4
   65 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   70 continue
      jbase=jbase+jump
   80 continue
      return
!
!     coding for factor 4
!
   90 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      id=ic+iink
      jd=jc+jink
      do 100 l=1,la
      i=ibase
      j=jbase
! 
      do 95 ijk=1,lot
      c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
      c(jc+j)=(a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))
      d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
      d(jc+j)=(b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i))
      c(jb+j)=(a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))
      c(jd+j)=(a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))
      d(jb+j)=(b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i))
      d(jd+j)=(b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i))
      i=i+inc3
      j=j+inc4
   95 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  100 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 120 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      kd=kc+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      c3=trigs(kd+1)
      s3=trigs(kd+2)
      do 110 l=1,la
      i=ibase
      j=jbase
! 
      do 105 ijk=1,lot
      c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
      d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
      c(jc+j)=c2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i)))   &
         -s2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
      d(jc+j)=s2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i)))   &
         +c2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
      c(jb+j)=c1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i)))   &
         -s1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
      d(jb+j)=s1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i)))   &
         +c1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
      c(jd+j)=c3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i)))   &
         -s3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
      d(jd+j)=s3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i)))   &
         +c3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
      i=i+inc3
      j=j+inc4
  105 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  110 continue
      jbase=jbase+jump
  120 continue
      return
!
!     coding for factor 5
!
  130 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      id=ic+iink
      jd=jc+jink
      ie=id+iink
      je=jd+jink
      do 140 l=1,la
      i=ibase
      j=jbase
! 
      do 135 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
      c(jb+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))   &
        -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
      c(je+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))   &
        +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
      d(jb+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))   &
        +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
      d(je+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))   &
        -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
      c(jc+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))   &
        -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
      c(jd+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))   &
        +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
      d(jc+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))   &
        +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
      d(jd+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))   &
        -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
      i=i+inc3
      j=j+inc4
  135 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  140 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 160 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      kd=kc+kb
      ke=kd+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      c3=trigs(kd+1)
      s3=trigs(kd+2)
      c4=trigs(ke+1)
      s4=trigs(ke+2)
      do 150 l=1,la
      i=ibase
      j=jbase
! 
      do 145 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
      c(jb+j)=c1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))   &
            -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))               &
         -s1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))       &
            +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      d(jb+j)=s1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))   &
            -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))               & 
         +c1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))       & 
            +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      c(je+j)=c4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))   &
            +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))               &
         -s4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))       &
            -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      d(je+j)=s4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))   &
            +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))               &
         +c4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))       &
            -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      c(jc+j)=c2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))   &
            -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))               &
         -s2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))       &
            +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      d(jc+j)=s2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))   &
            -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))               &
         +c2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))       &
            +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      c(jd+j)=c3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))   &
            +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))               &  
         -s3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))       &
            -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      d(jd+j)=s3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))   &
            +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))               &
         +c3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))       &
            -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      i=i+inc3
      j=j+inc4
  145 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  150 continue
      jbase=jbase+jump
  160 continue
      return
      end

      subroutine rcsexp(lfft,ifax,trig)
      dimension trig(1)
      integer ifax(1)
      call fftfax(lfft,ifax,trig)
      return
      end subroutine

      subroutine ccexp(lfft,ifax,trig)
      dimension trig(1)
      integer ifax(1)
      call cftfax(lfft,ifax,trig)
      return
      end subroutine
      subroutine rcosexp(lfft,sine,cosine,ifax,trig)
      dimension sine(1),cosine(1)
      dimension trig(1)
      integer ifax(1)
      
      tpi=8.*atan(1.)
      xnum = tpi/float(2*lfft)

      call fftfax(lfft,ifax,trig)
      do i=1,lfft
        sine(i)   = sin(float(i-1)*xnum)
        cosine(i) = cos(float(i-1)*xnum)
      end do 
      return
      end subroutine
      subroutine ccosexp(lfft,sine,cosine,ifax,trig)
      dimension sine(1),cosine(1)
      dimension trig(1)
      integer ifax(1) 

      tpi=8.*atan(1.)
      xnum = tpi/float(2*lfft)

      call cftfax(lfft,ifax,trig)
      do i=1,lfft
        sine(i)   = sin(float(i-1)*xnum)
        cosine(i) = cos(float(i-1)*xnum)
      end do 
      return
      end subroutine
      subroutine rsinexp(lfft,sine,cosine,ifax,trig)
      dimension sine(1),cosine(1)
      dimension trig(1)
      integer ifax(1) 

      tpi=8.*atan(1.)
      xnum = tpi/float(2*lfft)

      call fftfax(lfft,ifax,trig)
      do i=1,lfft
        sine(i)   = sin(float(i-1)*xnum)
        cosine(i) = cos(float(i-1)*xnum)
      end do 
       
      end subroutine
      subroutine csinexp(lfft,sine,cosine,ifax,trig)
      dimension sine(1),cosine(1)
      dimension trig(1)
      integer ifax(1)

      tpi=8.*atan(1.)
      xnum = tpi/float(2*lfft)

      call cftfax(lfft,ifax,trig)
      do i=1,lfft
        sine(i)   = sin(float(i-1)*xnum)
        cosine(i) = cos(float(i-1)*xnum)
      end do 
       
      end subroutine


      subroutine csr(a,b,work,lfft,n,ifax,trig)
      dimension a(1),b(1),work(1)
      integer ifax(1)
      dimension trig(1)
!cccccccccc   symmetric complex to real  fft  
!     note: lfft must be even. n is the number of ffts to be 
!     performed. the input , a , is a packed array and the 
!     output is returned in packed form to a . 
!       the inverse of csr is rcs. is = +1 always for csr. 
!ccccccc    first , unpack the data  cccccccccccccccccccccccccccc

      do k=1,n 

         loca = (k-1)*lfft       ! starting location or a
         locb = (k-1)*(lfft+2)   ! starting location for b
         loc1 = locb  + 2        ! imag part of first mode 
         loc2 = locb  + (lfft+1) ! real part of last mode 

         do i=1,lfft

            b(i+locb) = a(i+loca)

         end do


         b(loc2) = b(loc1)  !real part of last mode equals
         !imaginary part of first mode

         b(loc1) = 0.0   !imag part of first mode is zero
         b(loc2+1) = 0.0 !imag part of last mode is zero

      end do
!ccccccc    now perform the complex to real transform on b  cccccc
            is =   +1
            inc =   1
            jump =   lfft +2
            lot =   n 

         call rfftmlt(b,work,trig,ifax,inc,jump,lfft,lot,is)

!ccccccc   now pack the result  ccccccccccccccccccccccccccccccccccc
         do k =1,n 

            loca = (k-1)*lfft       ! starting location for a  
            locb = (k-1)*(lfft+2)   ! starting location for b 

            do i=1,lfft

               a(i+loca)  =  b(i+locb) 

            end do
         end do
           
      end subroutine 
      subroutine ccfft(a,br,bi,work,lfft,n,is,ifax,trig)
      dimension a(1),br(1),bi(1),work(1)
      integer ifax(1)
      dimension trig(1)
!cccccccccc  complex to complex fft  
!    note: lfft must be even and designates the fft length in (complex)
!    word pairs. n is the number of ffts to be performed.
!    performed. 
!    is = -1 or +1 for forward or inverse . 
!ccccccc    must avoid jump being a multiple of 2 ccccccccccccccccccc
              lfft2 = 2*lfft
              
              
              if (is.eq.-1) then 
                 fac = 1.0
              else
                 fac = (1./float(lfft))
              endif


              do k=1,n 
                 loca = (k-1)*lfft2     ! starting location or a
                 locb = (k-1)*(lfft+1)  ! starting location for br and bi

                 do i=1,lfft
                    locare  = loca + (2*i-1)
                    locaim   = locare+1

                    br(locb +i) = a(locare)
                    bi(locb +i) = a(locaim)
                 end do
              end do
!ccccccc    now perform complex to complex transform on b  cccccc
            inc =   1
            jump =   lfft +1       
            lot = n

         call cfftmlt(br,bi,work,trig,ifax,inc,jump,lfft,lot,is)
                           
!ccccccc   now pack the result  ccccccccccccccccccccccccccccccccccc
         do k=1,n 
            loca = (k-1)*lfft2     ! starting location or a
            locb = (k-1)*(lfft+1)  ! starting location for br and bi 

            do i=1,lfft
               locare  = loca + (2*i-1)
               locaim   = locare+1

               a(locare)=br(locb +i)*fac
               a(locaim)=bi(locb +i)*fac

            end do
         end do
        
      end subroutine 

      subroutine ccheb(a,b,br,bi,work,sumre,sumim,lfft,n,is,ifax,trig,sine,cosine)
        dimension a(1),b(1),br(1),bi(1),work(1)
        dimension trig(1),sine(1),cosine(1)
        dimension sumre(1),sumim(1)
        integer ifax(1) 

!       is = -1 real to chebyshev
!       is = +1 chebyshev to real                    
        lfftd2 = lfft/2
        lfft2  = 2*lfft
        fac1 = 4.0
        fac2 = 1./(8.*(float(lfft)))
!ccccccccc     preprocess the complex chebyshev coefficients  cc
        if (is.eq.+1) then 
           
           do k=1,n

              loca = (k-1)*(lfft2+2)

              a(loca+1) = 2.0*a(loca+1)
              a(loca+2) = 2.0*a(loca+2)
              a(loca+lfft2+1) = 2.0*a(loca+lfft2+1)
              a(loca+lfft2+2) = 2.0*a(loca+lfft2+2)

              do  i=1, lfft+1

                 locr =  loca + 2*i-1 
                 locim = locr+1

                 a(locr)= a(locr)*fac2
                 a(locim) = a(locim)*fac2

              end do
           end do
        endif

!ccccccccc compute complex cosine transform  cccccccccccccccccccccccc

        ! initialize sums 
        do k=1,n 
           sumr = 0.0 
           sumi = 0.0
           loca = (k-1)*(lfft2+2)  ! locations for a 
           loc1r = loca +1          ! first real loc for a 
           loc2r = loca +(lfft2+1) ! last real location for a   
           loc1i = loca +2          ! first imag. loc for a 
           loc2i = loca +(lfft2+2) ! last imag loc for a

           do i=2,lfft

              locar = loca + (2*i-1)  ! real locs for a 
              locai = locar +1        ! imag locs for a 

              sumr = sumr + a(locar)*cosine(i)
              sumi = sumi + a(locai)*cosine(i)
           end do

           sumre(k) = sumr + a(loc1r)/2. - a(loc2r)/2. 
           sumim(k) = sumi + a(loc1i)/2. - a(loc2i)/2. 
        end do

        !  pre-processing 

        do k=1,n 

           loca = (k-1)*(lfft2+2)    ! locs for a 
           locb = (k-1)*(lfft+1)    ! locs for br and bi 

           do i=1,lfft
              ir = 2*i-1
              im = ir+1

              ilowr  =  loca +ir                  ! 1:2l-1
              iupr   =  loca +(lfft2+2) - ir     ! 2l+1:3
              ilowim =  loca + im                 ! 2:2l
              iupim  =  loca +(lfft2+2) - (im-2) ! 2l+2:4


              br(locb +i) = (a(ilowr)+a(iupr))-2.0*sine(i)*(a(ilowr)-a(iupr))

              bi(locb +i) = (a(ilowim)+a(iupim))-2.0*sine(i)*(a(ilowim)-a(iupim))

           end do
        end do

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        ! complex to complex transform 
        isign = -1
        inc = 1 
        jump = lfft +1 
        lot = n

        call cfftmlt(br,bi,work,trig,ifax,inc,jump,lfft,lot,isign)    

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       
       do k=1,n 
          loca =   (k-1)*(lfft2+2)     ! starting location or a
          locb = (k-1)*(lfft+1)  ! starting location for real part of b

          do i=1,lfft
             locare  = loca + (2*i-1)
             locaim   = locare+1

             b(locare)=br(locb +i)
             b(locaim)=bi(locb +i)
          end do
       end do

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       ! post-processing 
       do k=1,n 

          loca = (k-1)*(lfft2+2)
          locb = (k-1)*(lfft2+2)

          ! fill real part of a 

          a(loca+1) = 2.* b(locb+1)
          a(loca+3) = fac1*sumre(k)
          a(loca+(lfft2+1)) = 2.*b(locb+lfft+1)

          ! fill imag part of a 

          a(loca+2) = 2.* b(locb+2)
          a(loca+4) = fac1*sumim(k)
          a(loca+(lfft2+2)) = 2.*b(locb+lfft+2)
       end do

       ! fill rest of a 

       do i=2,lfftd2
          ir = 2*i-1
          im = ir+1

          ir1 = 2*ir-1    ! 5:2l-3
          im1 = ir1+2     ! 7:2l-1
          ir2 = 2*ir      ! 6:2l-2
          im2 = ir2+2     ! 8:2l

          do k=1,n 

             loca = (k-1)*(lfft2+2)
             locb = (k-1)*(lfft2+2)

             ilowr  = locb +ir                  ! 3:l-1
             iupr   = locb +lfft2 -(ir-2) ! 2l-1:l-1
             ilowim = locb +im                  ! 4:l
             iupim  = locb +lfft2 - (ir-3)  ! 2l:l

             a(loca+ir1) = b(ilowr)+b(iupr)
             a(loca+im1) = a(loca+im1-4) - (b(ilowim)-b(iupim))
             a(loca+ir2) = b(ilowim)+b(iupim)
             a(loca+im2) = a(loca+im2-4) + (b(ilowr)-b(iupr))

          end do
       end do

!cccccc postprocessing on the cosine coefficients to get chebyshev coefficients 

       if (is.eq.-1) then 


          do k=1,n
             loca = (k-1)*(lfft2+2)

             a(loca+1)=  a(loca+1)/2.0
             a(loca+2) = a(loca+2)/2.0
             a(loca+lfft2+1) = a(loca+lfft2+1)/2.0
             a(loca+lfft2+2) = a(loca+lfft2+2)/2.0

          end do

       endif
         
      end subroutine 


      subroutine rcs(a,b,work,lfft,n,ifax,trig)
      dimension a(1),b(1),work(1)
      integer ifax(1)
      dimension trig(1) 
!ccccccc      real to symmetric complex fft  
!     note: lfft must be even. n is the number of ffts to be 
!     performed. the input , a , is a packed array and the 
!     output is returned in packed form to a . 
!       the inverse of rcs is csr. is = -1 always for rcs.  
!ccccccc    first , unpack the data  cccccccccccccccccccccccccccc

      do k=1,n 
         loca = (k-1)*lfft     ! starting location or a
         locb = (k-1)*(lfft+2) ! starting location for b
         do i=1,lfft
            b(i+locb) = a(i+loca)
         end do
      end do
 
!ccccccc    now perform the real to complex transform on b  cccccc
            is =   -1
            inc =   1
            jump =   lfft +2
            lot =   n 


          call rfftmlt(b,work,trig,ifax,inc,jump,lfft,lot,is)

!ccccccc   now pack the result  ccccccccccccccccccccccccccccccccccc
          
      
          do k =1,n

             loca = (k-1)*lfft       !starting location for a
             locb = (k-1)*(lfft+2)   !starting location for b

             loc1 = locb + (lfft+1)   !location for real part of last mode
             loc2 = locb + 2          !location for imag part of first mode        

             b(loc2) = b(loc1)


             do i=1,lfft

                a(i+loca)  =  b(i+locb) 

             end do
          end do
        
      end subroutine 
