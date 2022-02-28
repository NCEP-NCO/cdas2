      subroutine gtbhalf(ineofs,bhalf,bhalfp,jcap,nsig,nlath,a,
     *  jcapstat,nsigstat,agvz,wgvz,bvz,nmdszh,vz,vd,vh,vq,sigl)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    getbhalf   obtain bhalf and vert. functs.
c   prgmmr: parrish          org: w/nmc22    date: 90-10-06
c
c abstract: obtain bhalf, (sqrt(bhat)) and vertical functions
c
c program history log:
c   90-10-06  parrish
c
c   input argument list:
c     ineofs   - input file for input statistics
c     jcap     - triangular truncation
c     nsig     - number of sigma levels
c     nlath    - number of gaussian lats in one hemisphere
c     a        - a(4): scaling factors for forecast error spectra
c     jcapstat - triangular truncation for statistics
c     nsigstat - number of sigma levels for statistics
c     nmdszh   - number of vertical modes used in balance stats
c     sigl     - sigma layer values for model
c
c   output argument list:
c     bhalf    - sqrt(bhat) where bhat is background error spectrum
c     bhalfp   - sqrt(bhatp) where bhat is sur. press. back.error spectrum
c     agvz,wgvz,bvz - arrays to convert mass variable to t,ln(ps),div
c     vz,vd,vh,vq - vertical functions for vort,div, temp, and spec hum
c
c attributes:
c   language: cft77
c   machine:  cray ymp
c
c$$$
C-CRA          dimension sigl(nsig)
C-CRA          dimension aibw(nsig,nsig),w(nsig),vz(nsig,nsig)
C-CRA          dimension vd(nsig,nsig),vh(nsig,nsig),vq(nsig,nsig)
C-CRA          dimension chalf(0:jcapstat,0:jcapstat,nsig,4)
C-CRA          dimension chalfp(0:jcapstat,0:jcapstat)
C-CRA          dimension bhalf((jcap+1)*(jcap+2),nsig,4),a(4)
C-CRA          dimension bhalfp((jcap+1)*(jcap+2))
C-CRA          dimension beta(0:jcap,nsig),gamma(0:jcap,nsig)
C-CRA          dimension gammac(0:jcapstat,nsig)
C-CRA          dimension gammap(0:jcap),gammapc(0:jcapstat)
C-CRA          dimension betac(0:jcapstat,nsig)
C-CRA          dimension rlsg(nsig),tbar(nsig),a3(nsig,nsig)
C-CRA          dimension ainv(nsig,nsig)
C-CRA          dimension alphar(nsig),lwork(nsig),mwork(nsig)
C-CRA          dimension agvz(0:jcap,nsig,nmdszh),wgvz(0:jcap,nmdszh)
C-CRA          dimension bvz(0:jcap,nsig,nmdszh)

          dimension sigl(28)
          dimension aibw(28,28),w(28),vz(28,28)
          dimension vd(28,28),vh(28,28),vq(28,28)
          dimension chalf(0:126,0:126,28,4)
          dimension chalfp(0:126,0:126)
          dimension bhalf((62+1)*(62+2),28,4),a(4)
          dimension bhalfp((62+1)*(62+2))
          dimension beta(0:62,28),gamma(0:62,28)
          dimension gammac(0:126,28)
          dimension gammap(0:62),gammapc(0:126)
          dimension betac(0:126,28)
          dimension rlsg(28),tbar(28),a3(28,28)
          dimension ainv(28,28)
          dimension alphar(28),lwork(28),mwork(28)
          dimension agvz(0:62,28,28),wgvz(0:62,28)
          dimension bvz(0:62,28,28)
c--------
      print *,' read in vert eofs and horiz error spectra ',
     *     'from unit ',ineofs
      if(nsig.ne.nsigstat) then
      print *,' new vertical resolution, interpolate eofs'
      print *,' stopping code '
      end if
      rewind ineofs
      read(ineofs)msig,mlath,mmdszh,kcap,rlsg,aibw,w,
     *   rogc,tbar,a3,vz,vd,vh,vq,chalf,chalfp
c     This place passed
      read(ineofs)betac,gammac,gammapc
      rewind ineofs
      print *,' for vert eofs, nsig=',msig,
     *     ', nmdszh=',mmdszh,', jcap=',kcap,', nlath=',mlath
c    This place passed
C-CRA          beta=0.
c          dimension beta(0:62,28),gamma(0:62,28)
          DO j=1,nsig
          DO i=0,jcap
          beta(i,j)=0.
          ENDDO
          ENDDO
      do 190 k=1,nsig
      do 190 n=0,min(jcap,jcapstat)
        beta(n,k)=betac(n,k)
190   continue
c--------
c-------- now take care of horizontal part of stats
c--------
C-CRA          gamma=0
c          dimension beta(0:62,28),gamma(0:62,28)
          DO j=1,nsig
          DO i=0,jcap
          gamma(i,j)=0
          ENDDO
          ENDDO

c    This palce passed

      do 200 m=1,nsig
      do 200 n=0,min(jcap,jcapstat)
        gamma(n,m)=gammac(n,m)
200   continue
      do n=0,min(jcap,jcapstat)
       gammap(n)=gammapc(n)
      end do
C-CRA          bhalf=0.
c          dimension bhalf((62+1)*(62+2),28,4),a(4)
          DO i=1,(jcap+1)*(jcap+2)*nsig*4
          bhalf(i,1,1)=0.
          ENDDO
      do ll=1,4
       do k=1,nsig
        ii=-1
        do m=0,min(jcap,jcapstat)
         do l=0,min(jcap,jcapstat)-m
          ii=ii+2
          bhalf(ii,k,ll)=chalf(m,l,k,ll)
          bhalf(ii+1,k,ll)=bhalf(ii,k,ll)
         end do
        end do
       end do
      end do
C-CRA          bhalfp=0.
c          dimension bhalfp((62+1)*(62+2))
          DO i=1,(jcap+1)*(jcap+2)
          bhalfp(i)=0.
          ENDDO
      ii=-1
      do m=0,min(jcap,jcapstat)
       do l=0,min(jcap,jcapstat)-m
        ii=ii+2
        bhalfp(ii)=chalfp(m,l)
        bhalfp(ii+1)=bhalfp(ii)
       end do
      end do
      do 300 l=1,4
      do 300 i=1,(jcap+1)*(jcap+2)*nsig
        if(bhalf(i,1,l) .lt. 0.)then
        print *,' warning '
        print *,i,l,bhalf(i,1,l)
        end if
        bhalf(i,1,l)=a(l)*bhalf(i,1,l)
300   continue
      do i=1,(jcap+1)*(jcap+2)
        if(bhalfp(i) .lt. 0.)then
        print *,' warning surface'
        print *,i,bhalfp(i)
        end if
        bhalfp(i)=a(3)*bhalfp(i)
      end do
      do 400 k=1,nsig
        bhalf(1,k,1)=0.
        bhalf(1,k,2)=0.
400   continue
      do 500 i=1,nsig*4*(jcap+1)*(jcap+2)
        bhalf(i,1,1)=sqrt(bhalf(i,1,1))
500   continue
C-CRA          bhalfp=sqrt(bhalfp)
c          dimension bhalfp((62+1)*(62+2))
          DO i=1,(jcap+1)*(jcap+2)
          bhalfp(i)=sqrt(bhalfp(i))
          ENDDO


      close(ineofs)
c-------------------
C-CRA          agvz=0.
c          dimension agvz(0:62,28,28),wgvz(0:62,28)
          DO k=1,nmdszh
          DO j=1,nsig
          DO i=0,jcap
          agvz(i,j,k)=0.
          ENDDO
          ENDDO
          ENDDO
C-CRA          wgvz=0.
          DO k=1,nmdszh
          DO i=0,jcap
          wgvz(i,k)=0.
          ENDDO
          ENDDO
C-CRA          bvz=0.
c          dimension bvz(0:62,28,28)
          DO k=1,nmdszh
          DO j=1,nsig
          DO i=0,jcap
          bvz(i,j,k)=0.
          ENDDO
          ENDDO
          ENDDO
      do j=1,nsig
       do k=1,nsig
C-CRA            if(j.le.nmdszh)
C-CRA         *    wgvz(1:jcap,j)=wgvz(1:jcap,j)
C-CRA         *      +w(k)*gammap(1:jcap)*vz(k,j)
            if(j.le.nmdszh) then
              DO i=1,jcap
              wgvz(i,j)=wgvz(i,j)+w(k)*gammap(i)*vz(k,j)
              ENDDO
            endif
        do i=1,nsig
C-CRA             if(k.le.nmdszh)
C-CRA         *     agvz(1:jcap,j,k)=agvz(1:jcap,j,k)
C-CRA         *        +aibw(j,i)*gamma(1:jcap,j)*vz(i,k)
             if(k.le.nmdszh) then
               DO l=1,jcap
               agvz(l,j,k)=agvz(l,j,k)+aibw(j,i)*gamma(l,j)*vz(i,k)
               ENDDO
             endif
        end do
C-CRA            if(k.le.nmdszh)
C-CRA         *    bvz(1:jcap,j,k)=beta(1:jcap,j)*vz(j,k)
            if(k.le.nmdszh) then
              DO i=1,jcap
              bvz(i,j,k)=beta(i,j)*vz(j,k)
              ENDDO
            endif
       end do
      end do
      return
      end
