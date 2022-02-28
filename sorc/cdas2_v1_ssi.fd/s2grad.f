       subroutine s2grad(ds,u,v,jcap,nlon,nlath,qln,rln,trigs,ifax)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    s2grad     div to u,v
c   prgmmr: parrish          org: w/nmc22    date: 94-04-08
c
c abstract: spectral divergence coefs to grid u,v.
c
c program history log:
c   94-04-08  parrish
c
c   input argument list:
c     ds       - divergence coefficients
c     jcap     - triangular truncation
c     nlon     - number of longitudes
c     nlath    - number of gaussian lats in one hemisphere
c     qln      - q(n,l)
c     rln      - r(n,l)
c
c   output argument list:
c     u        - longitude component of winds
c     v        - latitude component of winds
c
c attributes:
c   language: cft77
c   machine:  cray ymp
c
c$$$
c
C-CRA             dimension ds((jcap+1)*(jcap+2))
C-CRA             dimension u(2*nlath+1,nlon+2)
C-CRA             dimension v(2*nlath+1,nlon+2)
C-CRA             dimension qln((jcap+1)*(jcap+2),nlath)
C-CRA             dimension rln((jcap+1)*(jcap+2),nlath)
C-CRA             dimension trigs(2*nlon),ifax(10)
C-CRA             dimension work(2*(2*nlath+1)*(nlon+2))
C-CRA             dimension ue(2*(jcap+1)),uo(2*(jcap+1))
C-CRA             dimension ve(2*(jcap+1)),vo(2*(jcap+1))
 
             dimension ds((62+1)*(62+2))
             dimension u(2*48+1,192+2)
             dimension v(2*48+1,192+2)
             dimension qln((62+1)*(62+2),48)
             dimension rln((62+1)*(62+2),48)
             dimension trigs(2*192),ifax(10)
             dimension work(2*(2*48+1)*(192+2))
             dimension ue(2*(62+1)),uo(2*(62+1))
             dimension ve(2*(62+1)),vo(2*(62+1))
c--------
c-------- internal scratch dynamic space follows:
c--------
c---------------
C-CRA                   v=0.
c          dimension v(2*nlath+1,nlon+2)
          DO ITMP=1,(2*nlath+1)*(nlon+2)
          v(ITMP,1)=0.
          ENDDO
         do j=1,nlath
          jr=2*nlath+1-j
          ii0=0
C-CRA                    ve=0.
c          dimension ve(2*(jcap+1)),vo(2*(jcap+1))
          DO ITMP=1,2*(jcap+1)
          ve(ITMP)=0.
          ENDDO
C-CRA                    vo=0.
c          dimension ve(2*(jcap+1)),vo(2*(jcap+1))
          DO ITMP=1,2*(jcap+1)
          vo(ITMP)=0.
          ENDDO
          do m=0,jcap,2
           do ll=1,2*(jcap+1-m)
            vo(ll)=vo(ll)-rln(ii0+ll,j)*ds(ii0+ll)
           end do
           if(m.lt.jcap) then
            ii0=ii0+2*(jcap+1-m)
            do ll=1,2*(jcap-m)
             ve(ll)=ve(ll)-rln(ii0+ll,j)*ds(ii0+ll)
            end do
            ii0=ii0+2*(jcap-m)
           end if
          end do
c----------
c---------- now combine even and odd parts
c----------
          do ll=1,2*(jcap+1)
           v(j,ll)=ve(ll)+vo(ll)
           v(jr,ll)=ve(ll)-vo(ll)
          end do
         end do
c--------
c-------- finally do fourier sums in longitude
c--------
         lot=nlath*2
         nlax=lot+1
C-CRA             call rfftmlt(v,work,trigs,ifax,nlax,1,nlon,lot,1)
             call fft99m (v,work,trigs,ifax,nlax,1,nlon,lot,1)
c---------------
c-------------multiply div by i
         do i=1,(jcap+1)*(jcap+2),2
          divr=ds(i)
          divi=ds(i+1)
          ds(i)=divi
          ds(i+1)=-divr
         end do
C-CRA                   u=0.
c          dimension u(2*nlath+1,nlon+2)
          DO ITMP=1,(2*nlath+1)*(nlon+2)
          u(ITMP,1)=0.
          ENDDO
         do j=1,nlath
          jr=2*nlath+1-j
          ii0=0
C-CRA                    ue=0.
c          dimension ue(2*(jcap+1)),uo(2*(jcap+1))
          DO ITMP=1,2*(jcap+1)
          ue(ITMP)=0.
          ENDDO
C-CRA                    uo=0.
c          dimension ue(2*(jcap+1)),uo(2*(jcap+1))
          DO ITMP=1,2*(jcap+1)
          uo(ITMP)=0.
          ENDDO
          do m=0,jcap,2
           do ll=1,2*(jcap+1-m)
            ue(ll)=ue(ll)+qln(ii0+ll,j)*ds(ii0+ll)
           end do
           if(m.lt.jcap) then
            ii0=ii0+2*(jcap+1-m)
            do ll=1,2*(jcap-m)
             uo(ll)=uo(ll)+qln(ii0+ll,j)*ds(ii0+ll)
            end do
            ii0=ii0+2*(jcap-m)
           end if
          end do
c----------
c---------- now combine even and odd parts
c----------
          do ll=1,2*(jcap+1)
           u(j,ll)=ue(ll)+uo(ll)
           u(jr,ll)=ue(ll)-uo(ll)
          end do
         end do
c--------
c-------- finally do fourier sums in longitude
c--------
C-CRA             call rfftmlt(u,work,trigs,ifax,nlax,1,nlon,lot,1)
             call fft99m (u,work,trigs,ifax,nlax,1,nlon,lot,1)
       return
       end
