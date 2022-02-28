       subroutine s2mg2x(ts,t,jcap,nlath,nlon,pln)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    s2mg2x     spectral variance to grid variance
c   prgmmr: parrish          org: w/nmc22    date: 90-12-07
c
c abstract: compute grid variance from spectral diag covariance.
c program history log:
c   90-12-07  parrish
c
c   input argument list:
c     ts       - spectral variance spectrum
c     jcap     - triangular truncation
c     nlath    - number of gaussian lats in one hemisphere
c     nlon     - number of longitudes
c     slat     - sin(gaussian latitudes)
c     pln      - spherical harmonics
c
c   output argument list:
c     t        - variance on grid
c
c attributes:
c   language: cft77
c   machine:  cray ymp
c
c$$$
c
C-CRA             dimension ts((jcap+1)*(jcap+2))
C-CRA             dimension t(2*nlath+1,nlon+2)
C-CRA             dimension pln((jcap+1)*(jcap+2),nlath)
C-CRA             dimension work(2*nlath+1,nlon+2)
C-CRA             dimension cos2(0:jcap,nlon),sin2(0:jcap,nlon)
C-CRA             dimension te(2*jcap+2),to(2*jcap+2)
 
             dimension ts((62+1)*(62+2))
             dimension t(2*48+1,192+2)
             dimension pln((62+1)*(62+2),48)
             dimension work(2*48+1,192+2)
             dimension cos2(0:62,192),sin2(0:62,192)
             dimension te(2*62+2),to(2*62+2)
c--------
c-------- internal scratch dynamic space follows:
c--------
c--------
         dlon=8.*atan(1.)/nlon
         do i=1,nlon
          angle=(i-1.)*dlon
          cos2(0,i)=1.
          sin2(0,i)=0.
          do l=1,jcap
           arg=l*angle
           cos2(l,i)=4.*cos(arg)**2
           sin2(l,i)=4.*sin(arg)**2
          end do
         end do
c--------
c-------- now sum in latitude (after zeroing output arrays)
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
            te(ll)=te(ll)+pln(ii0+ll,j)**2*ts(ii0+ll)
           end do
           if(m.lt.jcap) then
            ii0=ii0+2*(jcap+1-m)
            do ll=1,2*(jcap-m)
             to(ll)=to(ll)+pln(ii0+ll,j)**2*ts(ii0+ll)
            end do
            ii0=ii0+2*(jcap-m)
           end if
          end do
c----------
c---------- now combine even and odd parts
c----------
          do ll=1,2*(jcap+1)
           t(j,ll)=te(ll)+to(ll)
           t(jr,ll)=te(ll)+to(ll)
          end do
         end do
c--------
c------- finally do squared fourier sums in longitude
c-------
C-CRA                   work=t
c          dimension work(2*nlath+1,nlon+2)
          DO ITMP=1,(2*nlath+1)*(nlon+2)
          work(ITMP,1)=t(ITMP,1)
          ENDDO
C-CRA                   t=0.
c          dimension t(2*nlath+1,nlon+2)
          DO ITMP=1,(2*nlath+1)*(nlon+2)
          t(ITMP,1)=0.
          ENDDO
         do l=0,jcap
          lr=2*l+1
          li=2*l+2
          do i=1,nlon
           do j=1,2*nlath
            t(j,i)=t(j,i)+work(j,lr)*cos2(l,i)+work(j,li)*sin2(l,i)
           end do
          end do
         end do
       return
       end
