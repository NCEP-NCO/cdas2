      subroutine fulldivt(u,v,t,vort,div,plon,plat,z0,
     *       nsig,jcap,nlon,nlath,pln,qln,rln,trigs,ifax,del2,
     *       wgts,a3,sigl,sigi,jiter,ds,rlats)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    fulldivt    compute full field divergence tendency
c   prgmmr: parrish          org: w/nmc22    date: 94-02-11
c
c abstract: get divergence tendency (only dynamic terms so far)
c
c program history log:
c   94-02-11  parrish
c
c   input argument list:
c     u,v,t    - u, v, t on gaussian grid
c     vort,div - vort, div on gaussian grid
c     plon,plat- lon, lat derivatives of psfc on gaussian grid
c     z0       - terrain height on gaussian grid
c     nsig     - number of sigma levels
c     jcap     - triangular truncation of spectral coefficients
c     nlon     - num of longitudes
c     nlath    - number of gaussian lats in one hemisphere
c     rlats    - gaussian latitudes
c     trigs,ifax  - used by fft
c     del2     - constants for application of del**2 operator
c     wgts     - gaussian integration weights.
c     a3       - hydrostatic matrix
c     sigl,sigi- sigma coordinates
c     jiter    - outer iteration counter
c
c   output argument list:
c     div      - on output, contains sigdot, the vertical velocity
c     ds       - div tendency (spectral coefficients)
c
c attributes:
c   language: cft77
c   machine:  cray ymp
c
c$$$
c-----------
C-CRA          dimension sigl(nsig),sigi(nsig+1)
C-CRA          dimension vort(2*nlath+1,nlon+2,nsig)
C-CRA          dimension div(2*nlath+1,nlon+2,nsig)
C-CRA          real t(2*nlath+1,nlon+2,nsig)
C-CRA          dimension a3(nsig,nsig)
C-CRA          dimension del2((jcap+1)*(jcap+2))
C-CRA          dimension wgts(2*nlath)
C-CRA          dimension u(2*nlath+1,nlon+2,nsig)
C-CRA          dimension v(2*nlath+1,nlon+2,nsig)
C-CRA          dimension plon(2*nlath+1,nlon+2),plat(2*nlath+1,nlon+2)
C-CRA          dimension z0(2*nlath+1,nlon+2)
C-CRA          dimension trigs(nlon*2),ifax(10)
C-CRA          dimension ds((jcap+1)*(jcap+2),nsig)
C-CRA          dimension pln((jcap+1)*(jcap+2),nlath)
C-CRA          dimension qln((jcap+1)*(jcap+2),nlath)
C-CRA          dimension rln((jcap+1)*(jcap+2),nlath)
C-CRA          dimension rlats(2*nlath)
c--------
c-------- internal scratch dynamic space follows:
c--------
C-CRA          dimension ts((jcap+1)*(jcap+2),nsig)
C-CRA          dimension bige(2*nlath+1,nlon+2,nsig)
C-CRA          dimension ud(2*nlath+1,nlon+2),vd(2*nlath+1,nlon+2)
C-CRA          dimension uw(2*nlath+1,nlon+2,nsig),vw(2*nlath+1,nlon+2,nsig)
C-CRA          dimension coriolis(2*nlath+1,nlon+2)
C-CRA          dimension p(2*nlath+1,nlon+2)
C-CRA          dimension dsms(nsig)
C-CRA          dimension factor((jcap+1)*(jcap+2))
c-----------
          dimension sigl(28),sigi(28+1)
          dimension vort(2*48+1,192+2,28)
          dimension div(2*48+1,192+2,28)
          real t(2*48+1,192+2,28)
          dimension a3(28,28)
          dimension del2((62+1)*(62+2))
          dimension wgts(2*48)
          dimension u(2*48+1,192+2,28)
          dimension v(2*48+1,192+2,28)
          dimension plon(2*48+1,192+2)
          dimension plat(2*48+1,192+2)
          dimension z0(2*48+1,192+2)
          dimension trigs(192*2),ifax(10)
          dimension ds((62+1)*(62+2),28)
          dimension pln((62+1)*(62+2),48)
          dimension qln((62+1)*(62+2),48)
          dimension rln((62+1)*(62+2),48)
          dimension rlats(2*48)
c--------
c-------- internal scratch dynamic space follows:
c--------
          dimension ts((62+1)*(62+2),28)
          dimension bige(2*48+1,192+2,28)
          dimension ud(2*48+1,192+2),vd(2*48+1,192+2)
          dimension uw(2*48+1,192+2,28)
          dimension vw(2*48+1,192+2,28)
          dimension coriolis(2*48+1,192+2)
          dimension p(2*48+1,192+2)
          dimension dsms(28)
          dimension factor((62+1)*(62+2))
c--------
c-------- compute ud,vd, bige        
c--------
         ii=-1
         do m=0,jcap
          ii=ii+2
          factor(ii)=.5
          factor(ii+1)=0.
          if(m.lt.jcap) then
           do l=1,jcap-m
            ii=ii+2
            factor(ii)=1.
            factor(ii+1)=1.
           end do
          end if
         end do
         ng=(2*nlath+1)*nlon
         nc=(jcap+1)*(jcap+2)
C-CRA             coriolis=0.
             DO i=1,2*nlath+1*nlon+2
                coriolis(i,1)=0.
             ENDDO 
         omega=conmc('omega$')
C        CALL EXIT
         do j=1,nlath
          jr=2*nlath+1-j
C-CRA              coriolis(j,1:nlon)=2.*omega*sin(rlats(j))
C-CRA              coriolis(jr,1:nlon)=-coriolis(j,1:nlon)
              DO k=1,nlon
              coriolis(j,k)=2.*omega*sin(rlats(j))
              coriolis(jr,k)=-coriolis(j,k)
              ENDDO
         end do
         gascon=conmc('rd$')
         eaccel=9.8
c---------------------get sigdot  (overwrites divergence)
C-CRA             p=0.
             DO i=1,(2*nlath+1)*(nlon+2)
             p(i,1)=0.
             ENDDO
         do k=1,nsig
          do i=1,ng
           ak=(u(i,1,k)*plon(i,1)+v(i,1,k)*plat(i,1)
     *         +div(i,1,k))*(sigi(k)-sigi(k+1))
           p(i,1)=p(i,1)+ak
           div(i,1,k)=p(i,1)
          end do
         end do
         sigsum=0.
         do k=1,nsig
          sigsum=sigsum+sigi(k)-sigi(k+1)
          do i=1,ng
           div(i,1,k)=div(i,1,k)-sigsum*p(i,1)
          end do
         end do
c--------------------now compute div tendency (dynamics only)
         do k=1,nsig
          do j=1,2*nlath
           do i=1,nlon
            term1=v(j,i,k)*(vort(j,i,k)+coriolis(j,i))
            term2=gascon*t(j,i,k)*plon(j,i)
            term3=0.
            if(k.gt.1) 
     *        term3=term3+.5*div(j,i,k-1)*(u(j,i,k)-u(j,i,k-1))
     *            /(sigl(k)-sigl(k-1))
            if(k.lt.nsig)
     *        term3=term3+.5*div(j,i,k)*(u(j,i,k)-u(j,i,k+1))
     *            /(sigl(k)-sigl(k+1))
            ud(j,i)=term1-term2-term3
            term1=-u(j,i,k)*(vort(j,i,k)+coriolis(j,i))
            term2=gascon*t(j,i,k)*plat(j,i)
            term3=0.
            if(k.gt.1) 
     *        term3=term3+.5*div(j,i,k-1)*(v(j,i,k)-v(j,i,k-1))
     *            /(sigl(k)-sigl(k-1))
            if(k.lt.nsig)
     *        term3=term3+.5*div(j,i,k)*(v(j,i,k)-v(j,i,k+1))
     *            /(sigl(k)-sigl(k+1))
            vd(j,i)=term1-term2-term3
            term1=.5*(u(j,i,k)**2+v(j,i,k)**2)
            bige(j,i,k)=term1
           end do
          end do
C-CRA              uw(1:ng,1,k)=ud(1:ng,1)
C-CRA              vw(1:ng,1,k)=vd(1:ng,1)
              DO i=1,ng
              uw(i,1,k)=ud(i,1)
              vw(i,1,k)=vd(i,1)
              ENDDO
C-CRA              ud=eaccel*z0
              DO i=1,(2*nlath+1)*(nlon+2)
              ud(i,1)=eaccel*z0(i,1)
              ENDDO
          do l=1,nsig
C-CRA               ud(1:ng,1)=ud(1:ng,1)+eaccel*a3(k,l)*t(1:ng,1,l)
               DO i=1,ng
               ud(i,1)=ud(i,1)+eaccel*a3(k,l)*t(i,1,l)
               ENDDO
          end do
          do j=1,2*nlath
           do i=1,nlon
            bige(j,i,k)=bige(j,i,k)+ud(j,i)
           end do
          end do
         end do
c---------
c--------- now get div (ud,vd) - del2 (bige)
c---------
         do k=1,nsig
          call g2s0(ts(1,k),bige(1,1,k),jcap,nlon,nlath,
     *         wgts,pln,trigs,ifax)
          call grad2s(ds(1,k),uw(1,1,k),vw(1,1,k),jcap,nlon,nlath,
     *                qln,rln,trigs,ifax,wgts,del2)
C-CRA              ts(1:nc,k)=del2(1:nc)*ts(1:nc,k)
C-CRA              ds(1:nc,k)=ds(1:nc,k)+ts(1:nc,k)
              DO i=1,nc
              ts(i,k)=del2(i)*ts(i,k)
              ds(i,k)=ds(i,k)+ts(i,k)
              ENDDO
         end do
C-CRA             dsms=0.
             DO i=1,nsig
             dsms(i)=0.
             ENDDO
         do k=1,nsig
          do i=1,nc
           dsms(k)=dsms(k)+factor(i)*ds(i,k)**2
          end do
         end do
c        print *,' for outer iteration = ',jiter,' divt stats follow'
c        write(6,50)dsms
c50       format(1h ,6e13.4)
       return
       end
