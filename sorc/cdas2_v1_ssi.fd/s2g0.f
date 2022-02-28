       subroutine s2g0(ts,t,jcap,nlon,nlath,pln,trigs,ifax)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    s2g0       inverse of g2s0.
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
C-CRA             dimension te(2*jcap+2),to(2*jcap+2)
 
             dimension ts((62+1)*(62+2))
             dimension t(2*48+1,192+2)
             dimension trigs(192*2),ifax(10)
             dimension pln((62+1)*(62+2),48)
             dimension work(2*(2*48+1)*(192+2))
             dimension te(2*62+2),to(2*62+2)
c--------
C-CRA                   t=0.
c          dimension t(2*nlath+1,nlon+2)
          DO ITMP=1,(2*nlath+1)*(nlon+2)
          t(ITMP,1)=0.
          ENDDO
         do j=1,nlath
          jr=2*nlath+1-j
          ii0=0
C-CRA                    te=0.
c          dimension te(2*jcap+2),to(2*jcap+2)
          DO ITMP=1,2*jcap+2
          te(ITMP)=0.
          ENDDO
C-CRA                    to=0.
c          dimension te(2*jcap+2),to(2*jcap+2)
          DO ITMP=1,2*jcap+2
          to(ITMP)=0.
          ENDDO
          do m=0,jcap,2
           do ll=1,2*(jcap+1-m)
            te(ll)=te(ll)+pln(ii0+ll,j)*ts(ii0+ll)
           end do
           if(m.lt.jcap) then
            ii0=ii0+2*(jcap+1-m)
            do ll=1,2*(jcap-m)
             to(ll)=to(ll)+pln(ii0+ll,j)*ts(ii0+ll)
            end do
            ii0=ii0+2*(jcap-m)
           end if
          end do
c----------
c---------- now combine even and odd parts
c----------
          do ll=1,2*(jcap+1)
           t(j,ll)=te(ll)+to(ll)
           t(jr,ll)=te(ll)-to(ll)
          end do
         end do
c--------
c-------- finally do fourier sums in longitude
c--------
         lot=nlath*2
         nlax=lot+1
C-CRA             call rfftmlt(t,work,trigs,ifax,nlax,1,nlon,lot,1)
             call fft99m (t,work,trigs,ifax,nlax,1,nlon,lot,1)
       return
       end
