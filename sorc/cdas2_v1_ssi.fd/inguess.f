      subroutine inguess(gu,gv,gt,gp,gq,gmtns,sigi,sigl,
     *  inges,jcap,nsig,nlath,nlon,hourg,idateg,
     *   pln,qln,rln,trigs,ifax,
     *     ml2lm,factslm,factvlm)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    inguessv   same as inguess, but add vort,div, del(ps)
c   prgmmr: parrish          org: w/nmc22    date: 94-02-11
c
c abstract: augment inguess with vort, div, grad (log(psfc))
c
c program history log:
c   94-02-11  parrish
c
c   input argument list:
c     inges    - unit number of guess coefs
c     jcap     - triangular truncation
c     nsig     - number of sigma levels
c     nlath    - number of gaussian lats in one hemisphere
c     nlon     - number of longitudes
c     ap,bp,aqr,bqr,gr - recursion constants for spherical harmonics
c     slat,clat - sin and cos of gaussian latitudes
c     pe0,qe0,ro0 - starting functionf for spherical harmonic recursions
c     trigs,ifax - used by
c     del2     - n*(n+1)/a**2
c
c   output argument list:
c     gu       - guess u on grid
c     gv       - guess v on grid
c     gt       - guess t on grid
c     gp       - guess log(sfcp) on grid
c     gq       - guess specific humidity on grid
c     gmtns    - guess mountains
c     vortb,divb, plonb,platb - guess vort,div, grad(log(psfc))
c     sigi     - sigma values at interfaces of  sigma layers
c     sigl     - sigma values at mid-point of each sigma layer
c     hourg    - hour of guess field
c     idateg   - date of guess field
c
c attributes:
c   language: cft77
c   machine:  cray ymp
c
c$$$
C-CRA          dimension gu(2*nlath+1,nlon+2,nsig)
C-CRA          dimension gv(2*nlath+1,nlon+2,nsig)
C-CRA          dimension gt(2*nlath+1,nlon+2,nsig)
C-CRA          dimension gp(2*nlath+1,nlon+2)
C-CRA          dimension gq(2*nlath+1,nlon+2,nsig)
C-CRA          dimension gmtns(2*nlath+1,nlon+2)
C-CRA          dimension sigl(nsig),sigi(nsig+1)
C-CRA          integer idateg(4)
C-CRA          dimension del2((jcap+1)*(jcap+2))
C-CRA          dimension trigs(nlon*2),ifax(10)
C-CRA          dimension ml2lm((jcap+1)*(jcap+2))
C-CRA          dimension factslm((jcap+1)*(jcap+2))
C-CRA          dimension factvlm((jcap+1)*(jcap+2))
c--------
c-------- local space
c--------
C-CRA          real zc((jcap+1)*(jcap+2),nsig)
C-CRA          real dc((jcap+1)*(jcap+2),nsig)
C-CRA          real tc((jcap+1)*(jcap+2),nsig)
C-CRA          real qc((jcap+1)*(jcap+2),nsig)
C-CRA          real pc((jcap+1)*(jcap+2))
C-CRA          real rc((jcap+1)*(jcap+2))
C-CRA          character*4 on85(8)
C-CRA          real pln((jcap+1)*(jcap+2),nlath)
C-CRA          real qln((jcap+1)*(jcap+2),nlath)
C-CRA          real rln((jcap+1)*(jcap+2),nlath)

          dimension gu(2*48+1,192+2,28)
          dimension gv(2*48+1,192+2,28)
          dimension gt(2*48+1,192+2,28)
          dimension gp(2*48+1,192+2)
          dimension gq(2*48+1,192+2,28)
          dimension gmtns(2*48+1,192+2)
          dimension sigl(28),sigi(28+1)
          integer idateg(4)
          dimension del2((62+1)*(62+2))
          dimension trigs(192*2),ifax(10)
          dimension ml2lm((62+1)*(62+2))
          dimension factslm((62+1)*(62+2))
          dimension factvlm((62+1)*(62+2))
c--------
c-------- local space
c--------
          real zc((62+1)*(62+2),28)
          real dc((62+1)*(62+2),28)
          real tc((62+1)*(62+2),28)
          real qc((62+1)*(62+2),28)
          real pc((62+1)*(62+2))
          real rc((62+1)*(62+2))
          character*4 on85(8)
          real pln((62+1)*(62+2),48)
          real qln((62+1)*(62+2),48)
          real rln((62+1)*(62+2),48)
c--------
c-------- read in guess, putting into internal format.
c--------
      call rdgesc(zc,dc,tc,qc,pc,rc,hourg,idateg,sigi,sigl,
     *  inges,jcap,nsig,on85,ml2lm,factslm,factvlm)
c--------
c-------- reconstruct variables on grid
c--------
ccmic$ do all shared (nsig,pc,gp,jcap,nlon,nlath,pln,trigs,ifax)
ccmic$*       shared (rc,gmtns,tc,gt,qc,gq,zc,dc,gu,gv,qln,rln)
ccmic$*       private(kk,k)
         do kk=1,nsig*3+2
          if(kk.eq.3*nsig+1) 
     *      call s2g0(pc,gp,jcap,nlon,nlath,pln,trigs,ifax)
          if(kk.eq.3*nsig+2) 
     *      call s2g0(rc,gmtns,jcap,nlon,nlath,pln,trigs,ifax)
          k=mod(kk-1,nsig)+1
          if(kk.ge.1.and.kk.le.nsig) 
     *      call s2g0(tc(1,k),gt(1,1,k),jcap,nlon,nlath,pln,
     *                 trigs,ifax)
          if(kk.ge.nsig+1.and.kk.le.2*nsig)
     *      call s2g0(qc(1,k),gq(1,1,k),jcap,nlon,nlath,pln,
     *                 trigs,ifax)
          if(kk.ge.2*nsig+1.and.kk.le.3*nsig)
     *      call s2gvec(zc(1,k),dc(1,k),gu(1,1,k),gv(1,1,k),
     *       jcap,nlon,nlath,qln,rln,trigs,ifax)
         end do
      return
      end
