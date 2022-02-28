      subroutine initsat(msat,nlath,nlon,nsig,rt,cshat,pln,
     *       trigs,ifax,jcap,isatv,sfile)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    initsat    setup initial rhs for sat. temps      
c   prgmmr: parrish          org: w/nmc22    date: 90-10-11
c
c abstract: setup initial rhs for sat. temps  include horizontally
c                    and vertically correlated error
c
c program history log:
c   90-10-11  parrish
c   92-07-21
c
c   input argument list:
c     msat     - number of satellite profiles               
c     nlath    - number of gaussian lats in one hemisphere
c     nlon     - number of longitudes
c     nsig     - number of sigma levels
c     cshat    - diagonal spectral error covariance matrix (inverse)
c     pln      - p(n,l)
c     trigs,ifax - used by fft
c     jcap     - spectral trunctation
c     isatv    - array flags for use of sat. covariance
c     iscra    - unit number for conventional scratch unit
c     ntrecs   - number of records of conventional temperatures
c
c   output argument list:
c     rt       - output vector after inclusion of sat. info.
c
c attributes:
c   language: cft77
c   machine:  cray ymp
c
c$$$
c--------
c-------- multiply input vector t by inverse of correlated obs error
c--------
c
C-CRA          dimension cshat((jcap+1)*(jcap+2)),isatv(nsig)
C-CRA          dimension trigs(nlon*2),ifax(10)
C-CRA          dimension sfile(*)
C-CRA          dimension tval(nsig)
C-CRA          real rt(2*nlath+1,nlon+2,nsig)
C-CRA          real pln((jcap+1)*(jcap+2),nlath)
C-CRA          real rr(2*nlath+1,nlon+2,nsig)
 
          dimension cshat((62+1)*(62+2)),isatv(28)
          dimension trigs(192*2),ifax(10)
          dimension sfile(*)
          dimension tval(28)
          real rt(2*48+1,192+2,28)
          real pln((62+1)*(62+2),48)
          real rr(2*48+1,192+2,28)
c--------
c-------- local space
c-------
c     dimension indl(nsig)
c--------
      ngrd=(2*nlath+1)*(nlon+2)
C-CRA                rt=0.
c       real rt(2*nlath+1,nlon+2,nsig)
          DO ITMP=1,(2*nlath+1)*(nlon+2)*nsig
          rt(ITMP,1,1)=0.
          ENDDO
C-CRA                rr=0.
c       real rr(2*nlath+1,nlon+2,nsig)
          DO ITMP=1,(2*nlath+1)*(nlon+2)*nsig
          rr(ITMP,1,1)=0.
          ENDDO
      if(msat .eq. 0)go to 440
      iacnt=0
      do 610 ll=1,msat
      iacnt=iacnt+1
      numt=sfile(iacnt)
      print *,numt
      do 602 lll=1,numt
c     nlevs=sfile(1+iacnt)
      jlat=sfile(1+iacnt)
      jlon=sfile(2+iacnt)
      ibeg=sfile(3+iacnt)
      iacnt=iacnt+3
      nxig=nsig-ibeg+1
      call sgemv('T',nxig,nxig,-1.,sfile(iacnt+1+nsig),nxig,
     *     sfile(iacnt+1),1,1.,rr(jlat,jlon,ibeg),ngrd)
c     do 644 n=1,nsig
c44   tval(n)=sfile(iacnt+n)
c     iacnt=iacnt+nsig
c     do 133 llm=1,nsig
c     do 133 llx=1,nsig
c     rr(jlat,jlon,llm)=rr(jlat,jlon,llm)-
c    *     tval(llx)
c    *     *sfile((llm-1)*nsig+llx+iacnt)
c133  continue
c     iacnt=iacnt+nsig*nsig
      iacnt=iacnt+nsig*nsig+nsig
 602  continue
 610  continue
 611  continue
      call satc(rr,nsig,jcap,nlon,nlath,pln,trigs,ifax,cshat,isatv)
      iacnt=0
      do 510 ll=1,msat
      iacnt=iacnt+1
      numt=sfile(iacnt)
      do 102 lll=1,numt
c     nlevs=sfile(1+iacnt)
      jlat=sfile(1+iacnt)
      jlon=sfile(2+iacnt)
      ibeg=sfile(3+iacnt)
      iacnt=iacnt+3+nsig
      nxig=nsig-ibeg+1
      call sgemv('N',nxig,nxig,1.,sfile(iacnt+1),nxig,
     *     rr(jlat,jlon,ibeg),ngrd,1.,rt(jlat,jlon,ibeg),ngrd)
c     do 1133 llm=1,nsig
c     do 1133 llx=1,nsig
c     rt(jlat,jlon,llm)=rt(jlat,jlon,llm)+
c    *     rr(jlat,jlon,llx)*sfile((llx-1)*nsig+llm+iacnt)
c133  continue
      iacnt=iacnt+nsig*nsig
 102  continue
 510  continue
 511  continue
440   continue
      return
      end
