       subroutine tgrad2s(ds,u,v,jcap,nlon,nlath,qln,rln,trigs,ifax,
     *            wgts,del2)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    transpose of grad2s  
c   prgmmr: parrish          org: w/nmc22    date: 94-04-08
c
c abstract: spectral divergence coefs to grid u,v.
c
c program history log:
c   94-04-08  parrish
c   08-04-04  ebisuzaki, use f90 dynamic arrays
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
c   language: f90
c   machine:  AIX
c
c$$$
c
             dimension ds((jcap+1)*(jcap+2))
             dimension u(2*nlath+1,nlon+2)
             dimension v(2*nlath+1,nlon+2)
             dimension qln((jcap+1)*(jcap+2),nlath)
             dimension rln((jcap+1)*(jcap+2),nlath)
             dimension trigs(2*nlon),ifax(10)
             dimension wgts(2*nlath)
             dimension del2((jcap+1)*(jcap+2))
             dimension work(2*(2*nlath+1)*(nlon+2))
             dimension ue(2*(jcap+1)),uo(2*(jcap+1))
             dimension ve(2*(jcap+1)),vo(2*(jcap+1))
             dimension factor(2*jcap+2,nlath)
 
c--------
c-------- internal scratch dynamic space follows:
c--------
c--------
         do j=1,nlath
          factor(1,j)=wgts(j)/nlon
          factor(2,j)=0.
          do i=3,2*jcap+2
           factor(i,j)=.5*factor(1,j)
          end do
         end do
c--------
C-CRA                   u=0.
c          dimension u(2*nlath+1,nlon+2)
          DO ITMP=1,(2*nlath+1)*(nlon+2)
          u(ITMP,1)=0.
          ENDDO
C-CRA                   v=0.
c          dimension v(2*nlath+1,nlon+2)
          DO ITMP=1,(2*nlath+1)*(nlon+2)
          v(ITMP,1)=0.
          ENDDO
C-CRA                   ds=del2*ds
c          dimension ds((jcap+1)*(jcap+2))
          DO ITMP=1,(jcap+1)*(jcap+2)
          ds(ITMP)=del2(ITMP)*ds(ITMP)
          ENDDO
c---------------
         do j=1,nlath
C-CRA                    vo=0.
c          dimension ve(2*(jcap+1)),vo(2*(jcap+1))
          DO ITMP=1,2*(jcap+1)
          vo(ITMP)=0.
          ENDDO
C-CRA                    ve=0.
c          dimension ve(2*(jcap+1)),vo(2*(jcap+1))
          DO ITMP=1,2*(jcap+1)
          ve(ITMP)=0.
          ENDDO
          jr=2*nlath+1-j
          ii0=0
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
c---------- separate even and odd parts
          do ll=1,2*jcap+2
           v(j,ll)=(ve(ll)+vo(ll))*factor(ll,j)
           v(jr,ll)=(ve(ll)-vo(ll))*factor(ll,j)
          end do
         end do
c--------
c-------- next v, do fourier sums in longitude
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
c---------- separate even and odd parts
          do ll=1,2*jcap+2
           u(j,ll)=(ue(ll)+uo(ll))*factor(ll,j)
           u(jr,ll)=(ue(ll)-uo(ll))*factor(ll,j)
          end do
         end do
c--------
c-------- first do fourier analysis in longitude
c--------
C-CRA             call rfftmlt(u,work,trigs,ifax,nlax,1,nlon,lot,1)
             call fft99m (u,work,trigs,ifax,nlax,1,nlon,lot,1)
       return
       end
