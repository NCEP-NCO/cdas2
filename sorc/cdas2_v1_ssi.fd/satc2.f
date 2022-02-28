       subroutine satc(t,nsig,jcap,nlon,nlath,
     *            pln,trigs,ifax,cshat,isatv)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    satc      sat. error correlation operator
c   prgmmr: parrish          org: w/nmc22    date: 91-08-15
c
c abstract: operates satellite error correlation on fields
c
c program history log:
c   91-08-15  parrish
c
c   input argument list:
c     t        - gridded input temperature field
c     nsig     - number of sigma levels
c     jcap     - triangular truncation
c     nlon     - number of longitudes
c     nlath    - number of gaussian lats in one hemisphere
c     pln      - spherical harmonics
c     trigs,ifax - used by fft
c     cshat    - spectral weights
c     isatv    - flags for indicating whether or not to apply 
c                correlation, >0 yes, <0 no
c
c   output argument list:
c     t        - gridded output temperature field
c
c attributes:
c   language: cft77
c   machine:  cray ymp
c
c$$$
c
C-CRA          dimension cshat((jcap+1)*(jcap+2))
C-CRA          dimension pln((jcap+1)*(jcap+2),nlath)
C-CRA          dimension t(2*nlath+1,nlon+2,nsig)
C-CRA          dimension trigs(nlon*2),ifax(10),isatv(nsig)
C-CRA          dimension ts((jcap+1)*(jcap+2))
 
          dimension cshat((62+1)*(62+2))
          dimension pln((62+1)*(62+2),48)
          dimension t(2*48+1,192+2,28)
          dimension trigs(192*2),ifax(10),isatv(28)
          dimension ts((62+1)*(62+2))
c--------
c-------- internal scratch dynamic space follows:
c--------
c--------
      nbsig=1
      nesig=nsig
      do k=2,nsig
      if(isatv(k-1) .lt. 0.)then
      if(nbsig .eq. k-1)nbsig=k
      end if
      end do
      do k=nsig-1,1,-1
      if(isatv(k+1) .lt. 0.)then
      if(nesig .eq. k+1)nesig=k
      end if
      end do
      nsigl=nesig-nbsig+1
      if(nsigl .le. 0)then
      print *,' negative sigma levels in satc'
      stop
      end if
ccmic$  do all shared (t,jcap,nlon,nlath,pln,trigs,ifax)
ccmic$*        shared (nbsig,nesig,cshat)
ccmic$*        private (k,i,ts)
      do k=nbsig,nesig
       call ts2g0(ts,t(1,1,k),jcap,nlon,nlath,pln,trigs,ifax)
c-----------------------------multiply coefs by correlation
          do i=1,(jcap+1)*(jcap+2)
           ts(i)=ts(i)*cshat(i)
          end do
       call s2g0(ts,t(1,1,k),jcap,nlon,nlath,pln,trigs,ifax)
      end do
      return
      end
