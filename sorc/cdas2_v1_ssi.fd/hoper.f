       subroutine hoper(zs,ds,hs,qs,ps,u,v,vort,t,p,plon,plat,q,
     *     bhalf,bhalfp,nsig,jcap,nlon,nlath,del2,
     *     pln,qln,rln,trigs,ifax,
     *     agvz,wgvz,bvz,nmdszh,vz,vd,vh,vq,in,baln)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    hoper      analysis variables to grid variables
c   prgmmr: parrish          org: w/nmc22    date: 90-10-06
c
c abstract: convert analysis variables to grid variables
c
c program history log:
c   90-10-06  parrish
c   94-02-02  parrish
c
c   input argument list:
c     zs,ds,hs,qs,ps - coefs of vort, div, unbal t, unbal log(ps), q
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
c   output argument list:
c     u,v,vort,t,p,plon,plat,q - u,v,etc on gaussian grid
c
c attributes:
c   language: cft77
c   machine:  cray ymp
c
c$$$
C-CRA             real t(2*nlath+1,nlon+2,nsig),p(2*nlath+1,nlon+2)
C-CRA             real plon(2*nlath+1,nlon+2),plat(2*nlath+1,nlon+2)
C-CRA             dimension agvz(0:jcap,nsig,nmdszh)
C-CRA             dimension wgvz(0:jcap,nmdszh)
C-CRA             dimension bvz(0:jcap,nsig,nmdszh)
C-CRA             dimension u(2*nlath+1,nlon+2,nsig)
C-CRA             dimension v(2*nlath+1,nlon+2,nsig)
C-CRA             dimension vort(2*nlath+1,nlon+2,nsig)
C-CRA             dimension q(2*nlath+1,nlon+2,nsig)
C-CRA             dimension vz(nsig,nsig),vd(nsig,nsig)
C-CRA             dimension vq(nsig,nsig),vh(nsig,nsig)
C-CRA             dimension bhalf((jcap+1)*(jcap+2),nsig,4)
C-CRA             dimension bhalfp((jcap+1)*(jcap+2))
C-CRA             dimension zs((jcap+1)*(jcap+2),nsig)
C-CRA             dimension ds((jcap+1)*(jcap+2),nsig)
C-CRA             dimension hs((jcap+1)*(jcap+2),nsig)
C-CRA             dimension qs((jcap+1)*(jcap+2),nsig)
C-CRA             dimension ps((jcap+1)*(jcap+2))
C-CRA             dimension del2((jcap+1)*(jcap+2))
C-CRA             dimension trigs(nlon*2),ifax(10)
C-CRA             dimension pln((jcap+1)*(jcap+2),nlath)
C-CRA             dimension qln((jcap+1)*(jcap+2),nlath)
C-CRA             dimension rln((jcap+1)*(jcap+2),nlath)
C-CRA             dimension in((jcap+1)*(jcap+2))
C-CRA             dimension baln((jcap+1)*(jcap+2))
c--------
c-------- internal scratch dynamic space follows:
c--------
C-CRA             dimension psd((jcap+1)*(jcap+2))
C-CRA             dimension zsf((jcap+1)*(jcap+2),nmdszh)
C-CRA             dimension work((jcap+1)*(jcap+2),nsig)

             real t(2*48+1,192+2,28),p(2*48+1,192+2)
             real plon(2*48+1,192+2),plat(2*48+1,192+2)
             dimension agvz(0:62,28,28)
             dimension wgvz(0:62,28)
             dimension bvz(0:62,28,28)
             dimension u(2*48+1,192+2,28)
             dimension v(2*48+1,192+2,28)
             dimension vort(2*48+1,192+2,28)
             dimension q(2*48+1,192+2,28)
             dimension vz(28,28),vd(28,28)
             dimension vq(28,28),vh(28,28)
             dimension bhalf((62+1)*(62+2),28,4)
             dimension bhalfp((62+1)*(62+2))
             dimension zs((62+1)*(62+2),28)
             dimension ds((62+1)*(62+2),28)
             dimension hs((62+1)*(62+2),28)
             dimension qs((62+1)*(62+2),28)
             dimension ps((62+1)*(62+2))
             dimension del2((62+1)*(62+2))
             dimension trigs(192*2),ifax(10)
             dimension pln((62+1)*(62+2),48)
             dimension qln((62+1)*(62+2),48)
             dimension rln((62+1)*(62+2),48)
             dimension in((62+1)*(62+2))
             dimension baln((62+1)*(62+2))
c--------
c-------- internal scratch dynamic space follows:
c--------
             dimension psd((62+1)*(62+2))
             dimension zsf((62+1)*(62+2),28)
             dimension work((62+1)*(62+2),28)
c--------
         nc=(jcap+1)*(jcap+2)
         ng=(2*nlath+1)*(nlon+2)
c--------
c-------- first sum in vertical, and zero various arrays)
         do k=1,nsig
          if(k .eq. 1)then
C-CRA               p=0.
C-CRA               plon=0.
C-CRA               plat=0.
               DO i=1,(2*nlath+1)*(nlon+2)
               p(i,1)=0.
               plon(i,1)=0.
               plat(i,1)=0.
               ENDDO
           do i=1,nc
            ps(i)=ps(i)*bhalfp(i)
           end do
          end if
          do i=1,nc
           zs(i,k)=zs(i,k)*bhalf(i,k,1)
           ds(i,k)=ds(i,k)*bhalf(i,k,2)
           hs(i,k)=hs(i,k)*bhalf(i,k,3)
           qs(i,k)=qs(i,k)*bhalf(i,k,4)
          end do
          do i=1,ng
           u(i,1,k)=0.
           v(i,1,k)=0.
           vort(i,1,k)=0.
           t(i,1,k)=0.
           q(i,1,k)=0.
          end do
         end do
c------------------------apply spectral balance operator
c------------------------to zs
C-CRA            zsf=0.
            DO i=1,(jcap+1)*(jcap+2)*nmdszh
            zsf(i,1)=0.
            ENDDO
        do k=1,nmdszh
         ii0=2*(jcap+1)
         im0=0
         do m=1,jcap
          do ll=1,2*(jcap+1-m)
           zsf(ii0+ll,k)=zsf(ii0+ll,k)+baln(ii0+ll)*zs(im0+ll,k)
          end do
          ii0=ii0+2*(jcap+1-m)
          im0=im0+2*(jcap+2-m)
         end do
         ii0=0
         ip0=2*(jcap+1)
         do m=0,jcap-1
          do ll=1,2*(jcap-m)
           zsf(ii0+ll,k)=zsf(ii0+ll,k)+baln(ip0+ll)*zs(ip0+ll,k)
          end do
          ii0=ii0+2*(jcap+1-m)
          ip0=ip0+2*(jcap-m)
         end do
        end do
c---------------do temp and psfc
C-CRA             work=0.
             DO i=1,(jcap+1)*(jcap+2)*nsig
             work(i,1)=0.
             ENDDO
         do k=1,nsig
          if(k .eq. 1)then
           do j=1,nmdszh
            do i=1,nc
             ps(i)=ps(i)
     *             +wgvz(in(i),j)*zsf(i,j)
            end do
           end do
          end if
          do j=1,nmdszh
           do i=1,nc
            work(i,k)=work(i,k)
     *             +agvz(in(i),k,j)*zsf(i,j)
           end do
          end do
          do j=1,nsig
           do i=1,nc
            work(i,k)=work(i,k)+vh(k,j)*hs(i,j)
           end do
          end do
         end do
         do i=1,nsig*nc
          hs(i,1)=work(i,1)
          work(i,1)=0.
         end do
c--------------- sum in vertical ds                       
         do k=1,nsig
          do j=1,nsig
           do i=1,nc
            work(i,k)=work(i,k)+vd(k,j)*ds(i,j)
           end do
          end do
          do j=1,nmdszh
           do i=1,nc
            work(i,k)=work(i,k)+bvz(in(i),k,j)*zsf(i,j)
           end do
          end do
         end do
         do i=1,nc*nsig
          ds(i,1)=work(i,1)
          work(i,1)=0.
         end do
c--------
c-------- sum in vertical qs                       
         do k=1,nsig
          do j=1,nsig
           do i=1,nc
            work(i,k)=work(i,k)+vq(k,j)*qs(i,j)
           end do
          end do
         end do
         do i=1,nsig*nc
          qs(i,1)=work(i,1)
          work(i,1)=0.
         end do
c-------
c-------- sum in vertical zs                       
         do k=1,nsig
          do j=1,nsig
           do i=1,nc
            work(i,k)=work(i,k)+vz(k,j)*zs(i,j)
           end do
          end do
         end do
         do i=1,nsig*nc
          zs(i,1)=work(i,1)
         end do
         do i=1,nc
          psd(i)=-del2(i)*ps(i)
         end do
ccmic$ do all shared (nsig,psd,plon,plat,jcap,nlon,nlath,pln,qln,rln)
ccmic$*       shared (trigs,ifax,ps,p,zs,vort,hs,t,qs,q,zs,ds,u,v)
ccmic$*       private(kk,k)
         do kk=1,nsig*3+2
          if(kk.eq.3*nsig+1)
     *      call s2grad(psd,plon,plat,jcap,nlon,nlath,qln,rln,
     *           trigs,ifax)
          if(kk.eq.3*nsig+2)
     *      call s2g0(ps,p,jcap,nlon,nlath,pln,trigs,ifax)
          k=mod(kk-1,nsig)+1
          if(kk.ge.1.and.kk.le.nsig) then
           call s2g0(zs(1,k),vort(1,1,k),jcap,nlon,nlath,pln,
     *              trigs,ifax)
           call s2gvec(zs(1,k),ds(1,k),u(1,1,k),v(1,1,k),
     *        jcap,nlon,nlath,qln,rln,trigs,ifax)
          end if
          if(kk.ge.nsig+1.and.kk.le.2*nsig)
     *      call s2g0(hs(1,k),t(1,1,k),jcap,nlon,nlath,pln,
     *              trigs,ifax)
          if(kk.ge.2*nsig+1.and.kk.le.3*nsig)
     *      call s2g0(qs(1,k),q(1,1,k),jcap,nlon,nlath,pln,
     *              trigs,ifax)
         end do
       return
       end
