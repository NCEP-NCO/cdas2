      subroutine qoper(u,v,vort,t,plon,plat,nsig,jcap,nlon,nlath,
     *  pln,qln,rln,trigs,ifax,del2,wgts,a3,sigl,sigi,ds,iback,rlats,
     *   rus,rvs,rts,rvorts,rplons,rplats)
c-------------------
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    qoper   tangent linear model for divergency tendency
c   prgmmr: parrish          org: w/nmc22    date: 94-02-12
c
c abstract: tangent linear model (tlm) for divergence tendency.
c         special note: vertical advection terms not included yet
c-----                       (they are included in fulldivt computation)
c
c program history log:
c   94-02-12  parrish
c
c   input argument list:
c     u,v,vort,t,plon,plat - perturbation u,v, etc. on gaussian grid
c     nsig     - number of sigma layers
c     jcap     - triangular truncation
c     nlon     - number of longitudes
c     nlath    - number of gaussian lats in one hemisphere
c     pln,qln,rln - spherical harmonics
c     trigs,ifax - used by fft
c     del2     - n*(n+1)/a**2
c     wgts     - gaussian integration weights
c     a3       - hydrostatic matrix
c     sigl,sigi - vertical coordinate stuff
c     iback    - unit number where reference fields are stored
c
c   output argument list:
c     ds       - perturbation divergence tendency coefficients
c
c attributes:
c   language: cft77
c   machine:  cray ymp
c
c$$$
c
C-CRA          dimension sigl(nsig),sigi(nsig+1)
C-CRA          dimension a3(nsig,nsig)
C-CRA          dimension del2((jcap+1)*(jcap+2))
C-CRA          dimension wgts(2*nlath)
C-CRA          dimension u(2*nlath+1,nlon+2,nsig)
C-CRA          dimension v(2*nlath+1,nlon+2,nsig)
C-CRA          dimension vort(2*nlath+1,nlon+2,nsig)
C-CRA          dimension plon(2*nlath+1,nlon+2)
C-CRA          dimension plat(2*nlath+1,nlon+2)
C-CRA          dimension trigs(nlon*2),ifax(10)
C-CRA          dimension pln((jcap+1)*(jcap+2),nlath)
C-CRA          dimension qln((jcap+1)*(jcap+2),nlath)
C-CRA          dimension rln((jcap+1)*(jcap+2),nlath)
C-CRA          dimension ds((jcap+1)*(jcap+2),nsig)
C-CRA          dimension rlats(2*nlath)
C-CRA          dimension rus(2*nlath+1,nlon+2,nsig)
C-CRA          dimension rvs(2*nlath+1,nlon+2,nsig)
C-CRA          dimension rts(2*nlath+1,nlon+2,nsig)
C-CRA          dimension rvorts(2*nlath+1,nlon+2,nsig)
C-CRA          dimension rplons(2*nlath+1,nlon+2)
C-CRA          dimension rplats(2*nlath+1,nlon+2)
C-CRA          dimension ts((jcap+1)*(jcap+2),nsig)
C-CRA          dimension ps((jcap+1)*(jcap+2))
C-CRA          dimension uw(2*nlath+1,nlon+2,nsig),vw(2*nlath+1,nlon+2,nsig)
C-CRA          dimension pw(2*nlath+1,nlon+2)
C-CRA          dimension coriolis(2*nlath+1,nlon+2)
C-CRA          real t(2*nlath+1,nlon+2,nsig)
 
          dimension sigl(28),sigi(28+1)
          dimension a3(28,28)
          dimension del2((62+1)*(62+2))
          dimension wgts(2*48)
          dimension u(2*48+1,192+2,28)
          dimension v(2*48+1,192+2,28)
          dimension vort(2*48+1,192+2,28)
          dimension plon(2*48+1,192+2)
          dimension plat(2*48+1,192+2)
          dimension trigs(192*2),ifax(10)
          dimension pln((62+1)*(62+2),48)
          dimension qln((62+1)*(62+2),48)
          dimension rln((62+1)*(62+2),48)
          dimension ds((62+1)*(62+2),28)
          dimension rlats(2*48)
          dimension rus(2*48+1,192+2,28)
          dimension rvs(2*48+1,192+2,28)
          dimension rts(2*48+1,192+2,28)
          dimension rvorts(2*48+1,192+2,28)
          dimension rplons(2*48+1,192+2)
          dimension rplats(2*48+1,192+2)
          dimension ts((62+1)*(62+2),28)
          dimension ps((62+1)*(62+2))
          dimension uw(2*48+1,192+2,28),vw(2*48+1,192+2,28)
          dimension pw(2*48+1,192+2)
          dimension coriolis(2*48+1,192+2)
          real t(2*48+1,192+2,28)
c--------
c-------- internal scratch dynamic space follows:
c--------
c--------
         ng=(2*nlath+1)*nlon
         nc=(jcap+1)*(jcap+2)
         omega=conmc('omega$')
         gascon=conmc('rd$')
         eaccel=9.8
ccmic$  do all shared (nlath,coriolis,nlon,omega,rlats,uw,gascon,t)
ccmic$*        shared (rplons,rts,plon,rvs,vort,v,rvorts,vw,rus,rplats)
ccmic$*        shared (plat,u,eaccel,a3,t,nsig)
ccmic$*        private (j,k,i,l)
         do j=1,2*nlath
C-CRA                    coriolis(j,1:nlon)=2.*omega*sin(rlats(j))
          DO ITMP=1,nlon
                    coriolis(j,ITMP)=2.*omega*sin(rlats(j))
          END DO
c--------
c-------- compute ud,vd, bige        
c------------------------compute full non-lin bal eq  
          do k=1,nsig
           do i=1,nlon
            uw(j,i,k)=-gascon*t(j,i,k)*rplons(j,i)
     *           -gascon*rts(j,i,k)*plon(j,i)
     *           +rvs(j,i,k)*vort(j,i,k)
     *           +v(j,i,k)*(rvorts(j,i,k)+coriolis(j,i))
            vw(j,i,k)=-rus(j,i,k)*vort(j,i,k)
     *           -gascon*t(j,i,k)*rplats(j,i)
     *           -gascon*rts(j,i,k)*plat(j,i)
     *            -u(j,i,k)*(rvorts(j,i,k)+coriolis(j,i))
            vort(j,i,k)=u(j,i,k)*rus(j,i,k)
     *           +v(j,i,k)*rvs(j,i,k)
           end do
          end do
          do k=1,nsig
           do l=1,nsig
            do i=1,nlon
             vort(j,i,k)=vort(j,i,k)+eaccel*a3(k,l)*t(j,i,l)
            end do
           end do
          end do
         end do
c---------
c--------- now get div (ud,vd) - del2 (bige)  (bige stored in vort)
c---------
ccmic$ do all shared (nsig,ts,vort,jcap,nlon,nlath,wgts,pln,trigs,ifax)
ccmic$*       shared (ds,uw,vw,qln,rln,del2,nc)
ccmic$*       private (k)
         do k=1,nsig
          call g2s0(ts(1,k),vort(1,1,k),jcap,nlon,nlath,
     *         wgts,pln,trigs,ifax)
          call grad2s(ds(1,k),uw(1,1,k),vw(1,1,k),jcap,nlon,nlath,
     *      qln,rln,trigs,ifax,wgts,del2)
C-CRA                    ts(1:nc,k)=del2(1:nc)*ts(1:nc,k)
          DO ITMP=1,nc
                    ts(ITMP,k)=del2(ITMP)*ts(ITMP,k)
          END DO
C-CRA                    ds(1:nc,k)=ds(1:nc,k)+ts(1:nc,k)
          DO ITMP=1,nc
                    ds(ITMP,k)=ds(ITMP,k)+ts(ITMP,k)
          END DO
         end do
       return
       end
