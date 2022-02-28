       subroutine ts2grad(ds,u,v,jcap,nlon,nlath,qln,rln,trigs,ifax)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    transpose of s2grad  
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
C-CRA             dimension wgts(2*(jcap+1))
 
             dimension ds((62+1)*(62+2))
             dimension u(2*48+1,192+2)
             dimension v(2*48+1,192+2)
             dimension qln((62+1)*(62+2),48)
             dimension rln((62+1)*(62+2),48)
             dimension trigs(2*192),ifax(10)
             dimension work(2*(2*48+1)*(192+2))
             dimension ue(2*(62+1)),uo(2*(62+1))
             dimension ve(2*(62+1)),vo(2*(62+1))
             dimension wgts(2*(62+1))
c--------
c-------- internal scratch dynamic space follows:
c--------
c--------
C-CRA                   wgts=2.*nlon
c          dimension wgts(2*(jcap+1))
          DO ITMP=1,2*(jcap+1)
          wgts(ITMP)=2.*nlon
          ENDDO
         wgts(1)=nlon
         wgts(2)=0.
C-CRA                   ds=0.
c          dimension ds((jcap+1)*(jcap+2))
          DO ITMP=1,(jcap+1)*(jcap+2)
          ds(ITMP)=0.
          ENDDO
c--------
c-------- first do fourier analysis in longitude
c--------
         lot=nlath*2
         nlax=lot+1
C-CRA             call rfftmlt(u,work,trigs,ifax,nlax,1,nlon,lot,-1)
             call fft99m (u,work,trigs,ifax,nlax,1,nlon,lot,-1)
         do j=1,nlath
          jr=2*nlath+1-j
c---------- separate even and odd parts
          do ll=1,2*jcap+2
           ue(ll)=(u(j,ll)+u(jr,ll))*wgts(ll)
           uo(ll)=(u(j,ll)-u(jr,ll))*wgts(ll)
          end do
          ii0=0
          do m=0,jcap,2
           do ll=1,2*(jcap+1-m)
            ds(ii0+ll)=ds(ii0+ll)+qln(ii0+ll,j)*ue(ll)
           end do
           if(m.lt.jcap) then
            ii0=ii0+2*(jcap+1-m)
            do ll=1,2*(jcap-m)
             ds(ii0+ll)=ds(ii0+ll)+qln(ii0+ll,j)*uo(ll)
            end do
            ii0=ii0+2*(jcap-m)
           end if
          end do
         end do
c---------------
c-------------multiply div by i
         do i=1,(jcap+1)*(jcap+2),2
          divr=ds(i)
          divi=ds(i+1)
          ds(i)=-divi
          ds(i+1)=divr
         end do
c--------
c-------- next v, do fourier sums in longitude
c--------
C-CRA             call rfftmlt(v,work,trigs,ifax,nlax,1,nlon,lot,-1)
             call fft99m (v,work,trigs,ifax,nlax,1,nlon,lot,-1)
c---------------
         do j=1,nlath
          jr=2*nlath+1-j
c---------- separate even and odd parts
          do ll=1,2*jcap+2
           ve(ll)=(v(j,ll)+v(jr,ll))*wgts(ll)
           vo(ll)=(v(j,ll)-v(jr,ll))*wgts(ll)
          end do
          ii0=0
          do m=0,jcap,2
           do ll=1,2*(jcap+1-m)
            ds(ii0+ll)=ds(ii0+ll)-rln(ii0+ll,j)*vo(ll)
           end do
           if(m.lt.jcap) then
            ii0=ii0+2*(jcap+1-m)
            do ll=1,2*(jcap-m)
             ds(ii0+ll)=ds(ii0+ll)-rln(ii0+ll,j)*ve(ll)
            end do
            ii0=ii0+2*(jcap-m)
           end if
          end do
         end do
       return
       end
