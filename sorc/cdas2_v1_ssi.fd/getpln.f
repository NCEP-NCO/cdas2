       subroutine getpln(pln,qln,rln,jcap,nlath,ap,bp,slat,pe0,
     *          qe0,ro0,aqr,bqr,gr,clat,del2,del2out)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    getpln     generate legendre polynomials
c   prgmmr: parrish          org: w/nmc22    date: 90-09-21
c
c abstract: summation of scalar spherical harmonic series.
c
c program history log:
c   90-09-21  parrish
c   08-04-04  ebisuzaki, change some loops to f90 
c
c   input argument list:
c     jcap     - triangular truncation
c     nlath    - number of gaussian lats in one hemisphere
c     ap,bp    - recursion constants for spherical harmonics
c     slat     - sin(gaussian latitudes)
c     pe0      - starting functions for spherical harmonics
c
c   output argument list:
c     pln      - legendre polynomials
c
c attributes:
c   language: f90
c   machine:  AIX
c
c$$$
c
C-CRA             dimension ap(0:jcap,0:jcap),aqr(0:jcap,0:jcap)
C-CRA             dimension bp(0:jcap,0:jcap),bqr(0:jcap,0:jcap)
C-CRA             dimension gr(0:jcap,0:jcap)
C-CRA             dimension del2(0:jcap,0:jcap)
C-CRA             dimension slat(nlath),clat(nlath)
C-CRA             dimension pe0(nlath,0:jcap)
C-CRA             dimension qe0(nlath,0:jcap)
C-CRA             dimension ro0(nlath,0:jcap)
C-CRA             dimension pln((jcap+1)*(jcap+2),nlath)
C-CRA             dimension qln((jcap+1)*(jcap+2),nlath)
C-CRA             dimension rln((jcap+1)*(jcap+2),nlath)
C-CRA             dimension del2out((jcap+1)*(jcap+2))
C-CRA             dimension iadr(0:jcap,0:jcap)
C-CRA             real pe(nlath,0:jcap),po(nlath,0:jcap)
C-CRA             real qe(nlath,0:jcap),qo(nlath,0:jcap)
C-CRA             real re(nlath,0:jcap),ro(nlath,0:jcap)
 
             dimension ap(0:62,0:62),aqr(0:62,0:62)
             dimension bp(0:62,0:62),bqr(0:62,0:62)
             dimension gr(0:62,0:62)
             dimension del2(0:62,0:62)
             dimension slat(48),clat(48)
             dimension pe0(48,0:62)
             dimension qe0(48,0:62)
             dimension ro0(48,0:62)
             dimension pln((62+1)*(62+2),48)
             dimension qln((62+1)*(62+2),48)
             dimension rln((62+1)*(62+2),48)
             dimension del2out((62+1)*(62+2))
             dimension iadr(0:62,0:62)
             real pe(48,0:62),po(48,0:62)
             real qe(48,0:62),qo(48,0:62)
             real re(48,0:62),ro(48,0:62)
c--------
c-------- internal scratch dynamic space follows:
c--------
c--------
         ii=-1
         do m=0,jcap
          do l=0,jcap-m
           ii=ii+2
           iadr(l,m)=ii
          end do
         end do
         do m=0,jcap
          del2out(iadr(0,m))=del2(m,0)
          del2out(iadr(0,m)+1)=0.
          if(m.lt.jcap) then
           do l=1,jcap-m
            del2out(iadr(l,m))=del2(m,l)
            del2out(iadr(l,m)+1)=del2(m,l)
           end do
          end if
         end do
         pln=0.
         qln=0.
         rln=0.
         po=0.
         pe=pe0
         qo=0.
         qe=qe0
         re=0.
         ro=ro0
         do l=0,jcap
          do m=0,jcap-l,2
c------------ first even terms (m=0,2,...)
           do j=1,nlath
            pln(iadr(l,m),j)=pe(j,l)
            qln(iadr(l,m),j)=qe(j,l)
            rln(iadr(l,m),j)=ro(j,l)
           end do
c------------ now do odd  (m=1,3,...)
           if(m+1.le.jcap-l) then
            mp=m+1
            do j=1,nlath
              po(j,l)=ap(m,l)*slat(j)*pe(j,l)+bp(m,l)*po(j,l)
              qo(j,l)=aqr(m,l)*slat(j)*qe(j,l)
     *                     +bqr(m,l)*qo(j,l)
              re(j,l)=aqr(m,l)*slat(j)*ro(j,l)
     *                 +bqr(m,l)*re(j,l)+gr(m,l)*pe(j,l)*clat(j)
            end do
            do j=1,nlath
             pln(iadr(l,mp),j)=po(j,l)
             qln(iadr(l,mp),j)=qo(j,l)
             rln(iadr(l,mp),j)=re(j,l)
            end do
c-------------- get next pe
            do j=1,nlath
             pe(j,l)=ap(mp,l)*slat(j)*po(j,l)+bp(mp,l)*pe(j,l)
             qe(j,l)=aqr(mp,l)*slat(j)*qo(j,l)
     *                 +bqr(mp,l)*qe(j,l)
             ro(j,l)=aqr(mp,l)*slat(j)*re(j,l)
     *                 +bqr(mp,l)*ro(j,l)+gr(mp,l)*po(j,l)*clat(j)
            end do
           end if
          end do
         end do
         do i=1,(jcap+1)*(jcap+2)*nlath,2
          pln(i+1,1)=pln(i,1)
          qln(i+1,1)=qln(i,1)
          rln(i+1,1)=rln(i,1)
         end do
       return
       end
