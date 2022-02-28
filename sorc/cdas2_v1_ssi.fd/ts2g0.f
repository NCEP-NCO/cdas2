       subroutine ts2g0(ts,t,jcap,nlon,nlath,pln,trigs,ifax)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    ts2g0       transpose of s2g0
c   prgmmr: parrish          org: w/nmc22    date: 90-09-21
c
c abstract: summation of scalar spherical harmonic series.
c
c program history log:
c   90-09-21  parrish
c
c   input argument list:
c     ts       - spectral coefs
c     jcap     - triangular truncation
c     nlon     - number of longitudes
c     nlath    - number of gaussian lats in one hemisphere
c     pln      - spherical harmonics
c     trigs,ifax - used by fft
c
c   output argument list:
c     t        - values of desired field on gaussian grid
c
c attributes:
c   language: cft77
c   machine:  cray ymp
c
c$$$
c
C-CRA             dimension ts((jcap+1)*(jcap+2))
C-CRA             dimension t(2*nlath+1,nlon+2)
C-CRA             dimension trigs(nlon*2),ifax(10)
C-CRA             dimension pln((jcap+1)*(jcap+2),nlath)
C-CRA             dimension work(2*(2*nlath+1)*(nlon+2))
C-CRA             dimension wgts(2*jcap+2),te(2*jcap+2),to(2*jcap+2)
 
             dimension ts((62+1)*(62+2))
             dimension t(2*48+1,192+2)
             dimension trigs(192*2),ifax(10)
             dimension pln((62+1)*(62+2),48)
             dimension work(2*(2*48+1)*(192+2))
             dimension wgts(2*62+2),te(2*62+2),to(2*62+2)
c--------
C-CRA                   wgts=2.*nlon
c          dimension wgts(2*jcap+2),te(2*jcap+2),to(2*jcap+2)
          DO ITMP=1,2*jcap+2
          wgts(ITMP)=2.*nlon
          ENDDO
         wgts(1)=nlon
         wgts(2)=0.
c--------
c-------- first do fourier analysis in longitude
c--------
         lot=nlath*2
         nlax=lot+1
C-CRA             call rfftmlt(t,work,trigs,ifax,nlax,1,nlon,lot,-1)
             call fft99m (t,work,trigs,ifax,nlax,1,nlon,lot,-1)
C-CRA                   ts=0.
c          dimension ts((jcap+1)*(jcap+2))
          DO ITMP=1,(jcap+1)*(jcap+2)
          ts(ITMP)=0.
          ENDDO
         do j=1,nlath
          jr=2*nlath+1-j
c---------- separate even and odd parts
          do ll=1,2*jcap+2
           te(ll)=(t(j,ll)+t(jr,ll))*wgts(ll)
           to(ll)=(t(j,ll)-t(jr,ll))*wgts(ll)
          end do
          ii0=0
          do m=0,jcap,2
           do ll=1,2*(jcap+1-m)
            ts(ii0+ll)=ts(ii0+ll)+pln(ii0+ll,j)*te(ll)
           end do
           if(m.lt.jcap) then
            ii0=ii0+2*(jcap+1-m)
            do ll=1,2*(jcap-m)
             ts(ii0+ll)=ts(ii0+ll)+pln(ii0+ll,j)*to(ll)
            end do
            ii0=ii0+2*(jcap-m)
           end if
          end do
         end do
       return
       end
