      subroutine htoper(zs,ds,hs,qs,ps,u,v,vort,t,p,plon,plat,q,
     *     bhalf,bhalfp,nsig,jcap,nlon,nlath,del2,
     *     pln,qln,rln,trigs,ifax,
     *     agvz,wgvz,bvz,nmdszh,vz,vd,vh,vq,in,baln)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    htoper      transpose of hoper
c   prgmmr: parrish          org: w/nmc22    date: 90-10-06
c
c abstract: apply transpose of hoper, going from grid to spectral.
c
c program history log:
c   90-10-06  parrish
c
c   input argument list:
c     u,v,vort,t,p,plon,plat,q - u,v,etc on gaussian grid
c     bhalf    - background error stats
c     bhalfp   - background error stats surface pressure
c     nsig     - number of sigma levels
c     jcap     - triangular truncation
c     nlon     - number of longitudes
c     nlath    - number of gaussian lats in one hemisphere
c     del2     - n*(n+1)/a**2
c     trigs,ifax - used by fft
c     agvz     - mass-variable modes to temperature conversion
c     wgvz     - mass-variable modes to log(psfc) conversion
c     bvz      - mass-variable modes to divergence conversion
c     nmdszh   - number of modes used in balance eqn.
c     vz       - vertical mode matrix - z    
c     vd       - vertical mode matrix - d    
c     vh       - vertical mode matrix - temps
c     vq       - vertical mode matrix - q   
c     in       - total wavenumber index array
c     baln     - spectral balance operator constants
c
c
c   output argument list:
c     zs,ds,hs,qs,ps - coefs of vort, div, unbal t, unbal log(ps), q
c
c attributes:
c   language: cft77
c   machine:  cray ymp
c
c$$$
c
C-CRA          dimension agvz(0:jcap,nsig,nmdszh)
C-CRA          dimension wgvz(0:jcap,nmdszh)
C-CRA          dimension bvz(0:jcap,nsig,nmdszh)
C-CRA          dimension bhalf((jcap+1)*(jcap+2),nsig,4)
C-CRA          dimension bhalfp((jcap+1)*(jcap+2))
C-CRA          dimension zs((jcap+1)*(jcap+2),nsig)
C-CRA          dimension ds((jcap+1)*(jcap+2),nsig)
C-CRA          dimension hs((jcap+1)*(jcap+2),nsig)
C-CRA          dimension qs((jcap+1)*(jcap+2),nsig)
C-CRA          dimension ps((jcap+1)*(jcap+2))
C-CRA          dimension u(2*nlath+1,nlon+2,nsig)
C-CRA          dimension v(2*nlath+1,nlon+2,nsig)
C-CRA          dimension vort(2*nlath+1,nlon+2,nsig)
C-CRA          dimension q(2*nlath+1,nlon+2,nsig)
C-CRA          dimension vz(nsig,nsig),vd(nsig,nsig)
C-CRA          dimension vq(nsig,nsig),vh(nsig,nsig)
C-CRA          dimension del2((jcap+1)*(jcap+2))
C-CRA          dimension trigs(nlon*2),ifax(10)
C-CRA          dimension pln((jcap+1)*(jcap+2),nlath)
C-CRA          dimension qln((jcap+1)*(jcap+2),nlath)
C-CRA          dimension rln((jcap+1)*(jcap+2),nlath)
C-CRA          dimension in((jcap+1)*(jcap+2))
C-CRA          dimension baln((jcap+1)*(jcap+2))
C-CRA          dimension psd((jcap+1)*(jcap+2))
C-CRA          dimension zsf((jcap+1)*(jcap+2),nmdszh)
C-CRA          dimension work((jcap+1)*(jcap+2),nsig)
C-CRA          real t(2*nlath+1,nlon+2,nsig),p(2*nlath+1,nlon+2)
C-CRA          real plon(2*nlath+1,nlon+2),plat(2*nlath+1,nlon+2)
 
          dimension agvz(0:62,28,28)
          dimension wgvz(0:62,28)
          dimension bvz(0:62,28,28)
          dimension bhalf((62+1)*(62+2),28,4)
          dimension bhalfp((62+1)*(62+2))
          dimension zs((62+1)*(62+2),28)
          dimension ds((62+1)*(62+2),28)
          dimension hs((62+1)*(62+2),28)
          dimension qs((62+1)*(62+2),28)
          dimension ps((62+1)*(62+2))
          dimension u(2*48+1,192+2,28)
          dimension v(2*48+1,192+2,28)
          dimension vort(2*48+1,192+2,28)
          dimension q(2*48+1,192+2,28)
          dimension vz(28,28),vd(28,28)
          dimension vq(28,28),vh(28,28)
          dimension del2((62+1)*(62+2))
          dimension trigs(192*2),ifax(10)
          dimension pln((62+1)*(62+2),48)
          dimension qln((62+1)*(62+2),48)
          dimension rln((62+1)*(62+2),48)
          dimension in((62+1)*(62+2))
          dimension baln((62+1)*(62+2))
          dimension psd((62+1)*(62+2))
          dimension zsf((62+1)*(62+2),28)
          dimension work((62+1)*(62+2),28)
          real t(2*48+1,192+2,28),p(2*48+1,192+2)
          real plon(2*48+1,192+2),plat(2*48+1,192+2)
c--------
c-------- internal scratch dynamic space follows:
c--------
c--------
         nc=(jcap+1)*(jcap+2)
         ng=(2*nlath+1)*(nlon+2)
ccmic$ do all shared (nsig,psd,plon,plat,jcap,nlon,nlath,qln,rln,work)
ccmic$*     shared (trigs,ifax,ps,p,zs,vort,hs,t,qs,q,zs,ds,u,v,nc,pln)
ccmic$*       private(kk,k,i)
         do kk=1,nsig*3+2
          if(kk.eq.3*nsig+1)
     *      call ts2grad(psd,plon,plat,jcap,nlon,nlath,qln,rln,
     *            trigs,ifax)
          if(kk.eq.3*nsig+2)
     *      call ts2g0(ps,p,jcap,nlon,nlath,pln,trigs,ifax)
          k=mod(kk-1,nsig)+1
          if(kk.ge.1.and.kk.le.nsig) then
           call ts2g0(zs(1,k),vort(1,1,k),jcap,nlon,nlath,pln,
     *             trigs,ifax)
           call ts2gvec(work(1,k),ds(1,k),u(1,1,k),v(1,1,k),
     *          jcap,nlon,nlath,qln,rln,trigs,ifax)
           do i=1,nc
            zs(i,k)=zs(i,k)+work(i,k)
           end do
          end if
          if(kk.ge.nsig+1.and.kk.le.2*nsig)
     *      call ts2g0(hs(1,k),t(1,1,k),jcap,nlon,nlath,pln,
     *              trigs,ifax)
          if(kk.ge.2*nsig+1.and.kk.le.3*nsig)
     *      call ts2g0(qs(1,k),q(1,1,k),jcap,nlon,nlath,pln,
     *              trigs,ifax)
         end do
C-CRA                   ps=ps-del2*psd
c       dimension ps((jcap+1)*(jcap+2))
          DO ITMP=1,(jcap+1)*(jcap+2)
          ps(ITMP)=ps(ITMP)-del2(ITMP)*psd(ITMP)
          ENDDO
c--------
c-------- next do vertical transforms
c--------
c-------- tsum in vertical zs                       
      do j=1,nsig
         do i=1,nc
          work(i,j)=vz(1,j)*zs(i,1)
         end do
        do k=2,nsig
         do i=1,nc
          work(i,j)=vz(k,j)*zs(i,k)
     *             +work(i,j)
         end do
        end do
      end do
c--------
c-------- tsum in vertical qs                       
      do j=1,nsig
         do i=1,nc
          zs(i,j)=work(i,j)
         end do
         do i=1,nc
          work(i,j)=vq(1,j)*qs(i,1)
         end do
        do k=2,nsig
         do i=1,nc
          work(i,j)=vq(k,j)*qs(i,k)
     *             +work(i,j)
         end do
        end do
      end do
C-CRA                qs=work
c       dimension qs((jcap+1)*(jcap+2),nsig)
          DO ITMP=1,((jcap+1)*(jcap+2))*nsig
          qs(ITMP,1)=work(ITMP,1)
          ENDDO
C-CRA                work=ds
c       dimension work((jcap+1)*(jcap+2),nsig)
          DO ITMP=1,((jcap+1)*(jcap+2))*nsig
          work(ITMP,1)=ds(ITMP,1)
          ENDDO
C-CRA                ds=0.
c       dimension ds((jcap+1)*(jcap+2),nsig)
          DO ITMP=1,((jcap+1)*(jcap+2))*nsig
          ds(ITMP,1)=0.
          ENDDO
C-CRA                zsf=0.
c       dimension zsf((jcap+1)*(jcap+2),nmdszh)
          DO ITMP=1,((jcap+1)*(jcap+2))*nmdszh
          zsf(ITMP,1)=0.
          ENDDO
c--------------- tsum in vertical ds                       
         do j=1,nsig
          do k=1,nsig
           do i=1,nc
            ds(i,j)=ds(i,j)+work(i,k)*vd(k,j)
           end do
          end do
          if(j.le.nmdszh) then
           do k=1,nsig
            do i=1,nc
             zsf(i,j)=zsf(i,j)+bvz(in(i),k,j)*work(i,k)
            end do
           end do
          end if
         end do
c---------------do ttemp and tpsfc
C-CRA                   work=hs
c       dimension work((jcap+1)*(jcap+2),nsig)
          DO ITMP=1,((jcap+1)*(jcap+2))*nsig
          work(ITMP,1)=hs(ITMP,1)
          ENDDO
C-CRA                   hs=0.
c       dimension hs((jcap+1)*(jcap+2),nsig)
          DO ITMP=1,((jcap+1)*(jcap+2))*nsig
          hs(ITMP,1)=0.
          ENDDO
         do j=1,nsig
          if(j.le.nmdszh) then
           do i=1,nc
            zsf(i,j)=zsf(i,j)+wgvz(in(i),j)*ps(i)
           end do
           do k=1,nsig
            do i=1,nc
             zsf(i,j)=zsf(i,j)+agvz(in(i),k,j)*work(i,k)
            end do
           end do
          end if
          do k=1,nsig
           do i=1,nc
            hs(i,j)=hs(i,j)+work(i,k)*vh(k,j)
           end do
          end do
         end do
c------------------------tapply spectral balance operator
c------------------------to zs
        do k=1,nmdszh
         ii0=2*(jcap+1)
         im0=0
         do m=1,jcap
          do ll=1,2*(jcap+1-m)
           zs(im0+ll,k)=zs(im0+ll,k)+baln(ii0+ll)*zsf(ii0+ll,k)
          end do
          ii0=ii0+2*(jcap+1-m)
          im0=im0+2*(jcap+2-m)
         end do
         ii0=0
         ip0=2*(jcap+1)
         do m=0,jcap-1
          do ll=1,2*(jcap-m)
           zs(ip0+ll,k)=zs(ip0+ll,k)+baln(ip0+ll)*zsf(ii0+ll,k)
          end do
          ii0=ii0+2*(jcap+1-m)
          ip0=ip0+2*(jcap-m)
         end do
        end do
      do j=1,nsig
       if(j .eq. 1)then
         do i=1,nc
           ps(i)=ps(i)*bhalfp(i)
         end do
       end if
       do i=1,nc
        zs(i,j)=zs(i,j)*bhalf(i,j,1)
        ds(i,j)=ds(i,j)*bhalf(i,j,2)
        hs(i,j)=hs(i,j)*bhalf(i,j,3)
        qs(i,j)=qs(i,j)*bhalf(i,j,4)
       end do
      end do
      return
      end
