      subroutine qtoper(u,v,vort,t,plon,plat,nsig,jcap,nlon,nlath,
     *  pln,qln,rln,trigs,ifax,del2,wgts,a3,sigl,sigi,ds,iback,rlats,
     *   rus,rvs,rts,rvorts,rplons,rplats)
c-----------------
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    qtoper   adjoint of qoper, tlm for divergency tendency
c   prgmmr: parrish          org: w/nmc22    date: 94-02-12
c
c abstract: adjoint of tangent linear model (tlm) for div tendency.
c         special note: vertical advection terms not included here 
c                          (they are included in fulldivt)
c
c program history log:
c   94-02-12  parrish
c
c   input argument list:
c     ds       - perturbation divergence tendency coefficients
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
c     u,v,vort,t,plon,plat - perturbation u,v, etc. on gaussian grid
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
C-CRA          dimension plon(2*nlath+1,nlon+2),plat(2*nlath+1,nlon+2)
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
C-CRA          dimension pw(2*nlath+1,nlon+2)
C-CRA          dimension ts((jcap+1)*(jcap+2),nsig)
C-CRA          dimension ps((jcap+1)*(jcap+2))
C-CRA          dimension uw(2*nlath+1,nlon+2,nsig),vw(2*nlath+1,nlon+2,nsig)
C-CRA          dimension coriolis(2*nlath+1,nlon+2)
C-CRA          real t(2*nlath+1,nlon+2,nsig)
 
          dimension sigl(28),sigi(28+1)
          dimension a3(28,28)
          dimension del2((62+1)*(62+2))
          dimension wgts(2*48)
          dimension u(2*48+1,192+2,28)
          dimension v(2*48+1,192+2,28)
          dimension vort(2*48+1,192+2,28)
          dimension plon(2*48+1,192+2),plat(2*48+1,192+2)
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
          dimension pw(2*48+1,192+2)
          dimension ts((62+1)*(62+2),28)
          dimension ps((62+1)*(62+2))
          dimension uw(2*48+1,192+2,28)
          dimension vw(2*48+1,192+2,28)
          dimension coriolis(2*48+1,192+2)
          real t(2*48+1,192+2,28)
c--------
c-------- internal scratch dynamic space follows:
c--------
c--------
c-------- compute ud,vd, bige  (stored in vort)       
c--------
         ng=(2*nlath+1)*nlon
         nc=(jcap+1)*(jcap+2)
         omega=conmc('omega$')
         gascon=conmc('rd$')
         eaccel=9.8
c---------
ccmic$ do all shared (nsig,ts,vort,jcap,nlon,nlath,wgts,pln,trigs,ifax)
ccmic$*       shared (ds,uw,vw,qln,rln,del2,nc)
ccmic$*       private (k)
         do k=1,nsig
C-CRA                    ts(1:nc,k)=del2(1:nc)*ds(1:nc,k)
          DO ITMP=1,nc
                    ts(ITMP,k)=del2(ITMP)*ds(ITMP,k)
          END DO
          call tg2s0(ts(1,k),vort(1,1,k),jcap,nlon,nlath,
     *         wgts,pln,trigs,ifax)
          call tgrad2s(ds(1,k),uw(1,1,k),vw(1,1,k),jcap,nlon,nlath,
     *        qln,rln,trigs,ifax,wgts,del2)
         end do
c        print *,' in qtoper at 1, ',uw(1,1,1),vw(1,1,1),vort(1,1,1)
c---------
ccmic$  do all shared (nlath,coriolis,nlon,omega,rlats,uw,gascon,t)
ccmic$*        shared (rplons,rts,plon,rvs,vort,v,rvorts,vw,rus,rplats)
ccmic$*        shared (plat,u,eaccel,a3,t,nsig)
ccmic$*        private (j,k,i,l)
         do j=1,2*nlath
C-CRA                    coriolis(j,1:nlon)=2.*omega*sin(rlats(j))
          DO ITMP=1,nlon
                    coriolis(j,ITMP)=2.*omega*sin(rlats(j))
          END DO
          do k=1,nsig
           do l=1,nsig
            do i=1,nlon
             t(j,i,l)=t(j,i,l)+eaccel*a3(k,l)*vort(j,i,k)
            end do
           end do
          end do
c---------------------
          do i=1,nlon
           plon(j,i)=0.
           plat(j,i)=0.
          end do
          do k=1,nsig
           do i=1,nlon
            t(j,i,k)=t(j,i,k)-gascon*uw(j,i,k)*rplons(j,i)
     *           -gascon*vw(j,i,k)*rplats(j,i)
            u(j,i,k)=u(j,i,k)+vort(j,i,k)*rus(j,i,k)
     *            -vw(j,i,k)*(rvorts(j,i,k)+coriolis(j,i))
            v(j,i,k)=v(j,i,k)
     *           +vort(j,i,k)*rvs(j,i,k)
     *           +uw(j,i,k)*(rvorts(j,i,k)+coriolis(j,i))
            plon(j,i)=plon(j,i)
     *           -gascon*rts(j,i,k)*uw(j,i,k)
            plat(j,i)=plat(j,i)
     *           -gascon*rts(j,i,k)*vw(j,i,k)
            vort(j,i,k)=-rus(j,i,k)*vw(j,i,k)
     *           +rvs(j,i,k)*uw(j,i,k)
           end do
          end do
         end do
c------------------------compute full non-lin bal eq 
       return
       end
