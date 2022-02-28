      subroutine rdgesc(zc,dc,tc,qc,pc,rc,hourg,idateg,sigi,sigl,
     *  inges,jcap,nsig,on85,ml2lm,factslm,factvlm)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    rdgesc     read sigma coefs and reorder.
c   prgmmr: parrish          org: w/nmc22    date: 90-10-10
c
c abstract: read guess sigma coefs, and reorder to internal format.
c
c program history log:
c   90-10-10  parrish
c
c   input argument list:
c     inges    - unit number of guess coefs
c     jcap     - triangular truncation
c     nsig     - number of sigma levels
c
c   output argument list:
c     zc,dc,tc,qc,pc,rc - ges sig coefs of vort,div,t,q,ln(ps),z0
c     hourg    - guess forecast hour
c     idateg   - initial date of guess
c     sigi     - sigma values at interface of each sigma layer
c     sigl     - sigma values at mid-point of each sigma layer
c     on85     - on85 date record for guess coefs
c
c attributes:
c   language: cft77
c   machine:  cray ymp
c
c$$$
c
C-CRA          dimension zc((jcap+1)*(jcap+2),nsig)
C-CRA          dimension dc((jcap+1)*(jcap+2),nsig)
C-CRA          dimension tc((jcap+1)*(jcap+2),nsig)
C-CRA          dimension qc((jcap+1)*(jcap+2),nsig)
C-CRA          dimension pc((jcap+1)*(jcap+2))
C-CRA          dimension rc((jcap+1)*(jcap+2))
C-CRA          dimension idateg(4),sigi(nsig+1),sigl(nsig)
C-CRA          dimension ml2lm((jcap+1)*(jcap+2))
C-CRA          dimension factslm((jcap+1)*(jcap+2))
C-CRA          dimension factvlm((jcap+1)*(jcap+2))
C-CRA          dimension z((jcap+1)*(jcap+2))
C-CRA          character*4 on85(8)
 
          dimension zc((62+1)*(62+2),28)
          dimension dc((62+1)*(62+2),28)
          dimension tc((62+1)*(62+2),28)
          dimension qc((62+1)*(62+2),28)
          dimension pc((62+1)*(62+2))
          dimension rc((62+1)*(62+2))
          dimension idateg(4),sigi(28+1),sigl(28)
          dimension ml2lm((62+1)*(62+2))
          dimension factslm((62+1)*(62+2))
          dimension factvlm((62+1)*(62+2))
          dimension z((62+1)*(62+2))
          character*4 on85(8)
c--------
c-------- local space
c--------
c-------
      nc=(jcap+1)*(jcap+2)
      rewind inges
c-------- hour,idate, etc.
      read(inges)on85
      read(inges)hourg,idateg,sigi,sigl
c-------- terrain coefs
      read(inges)z
      do i=1,nc
       rc(i)=factslm(i)*z(ml2lm(i))
      end do
c-------- sfcp coefficients
      read(inges)z
      do i=1,nc
       pc(i)=factslm(i)*z(ml2lm(i))
      end do
c-------- temp coefficients
      do k=1,nsig
        read(inges)z
        do i=1,nc
         tc(i,k)=factslm(i)*z(ml2lm(i))
        end do
      end do
c-------- div and vort
      do k=1,nsig
        read(inges)z
        do i=1,nc
         dc(i,k)=factvlm(i)*z(ml2lm(i))
        end do
        read(inges)z
        do i=1,nc
         zc(i,k)=factvlm(i)*z(ml2lm(i))
        end do
      end do
c-------- q coefs
      do k=1,nsig
        read(inges)z
        do i=1,nc
         qc(i,k)=factslm(i)*z(ml2lm(i))
        end do
      end do
      write(6,700)jcap,nsig,hourg,idateg
700   format(' guess sigma coefficients read in, jcap,nsig=',
     *  2i6,/,' hour,idate=',f10.1,4i4)
      close (inges)
      return
      end
