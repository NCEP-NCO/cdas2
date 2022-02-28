       subroutine gtbdivt(idivt,bdivt,nsig,jcap,nsigdivt,
     *            jcapdivt,sigl)
c-------------
c-------------bring in estimates of divtend error variance,
c-----------  interpolate in vertical, and truncate in horizontal
c-----------
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    gtbdivt    read divtend error variances.
c   prgmmr: parrish          org: w/nmc22    date: 94-02-11
c
c abstract: read divtend error variance, interpolate in vert, truncate
c            in horizontal, as necessary.
c
c program history log:
c   94-02-11  parrish
c
c   input argument list:
c     idivt    - input unit number containing divtend error variances
c     nsig     - number of sigma levels
c     jcap     - triangular truncation
c     nsigdivt - number of sigma levels in input divtend errors
c     jcapdivt - triangular truncation for input divtend errors
c     sigl     - analysis sigma levels
c
c   output argument list:
c     bdivt    - inverse of divtend error variances.
c
c attributes:
c   language: cft77
c   machine:  cray ymp
c
c$$$
C-CRA             dimension bdivt(0:jcap,nsig),sigl(nsig)
C-CRA             dimension cdivt(0:jcapdivt,nsigdivt),sigldivt(nsigdivt)
C-CRA             dimension tdivt(0:jcapdivt,nsig)
C-CRA             dimension rlsg(nsig),rlsgdivt(nsigdivt)
C-CRA             dimension grid(nsig)

             dimension bdivt(0:62,28),sigl(28)
             dimension cdivt(0:126,28)
             dimension sigldivt(28)
             dimension tdivt(0:126,28)
             dimension rlsg(28),rlsgdivt(28)
             dimension grid(28)
c-----------
         rewind idivt
         read(idivt)msigdivt,mcapdivt,sigldivt,cdivt
         close(idivt)
C-CRA             rlsgdivt=log(sigldivt)
             DO i=1,nsigdivt
             rlsgdivt(i)=log(sigldivt(i))
             ENDDO
C-CRA             rlsg=log(sigl)
             DO i=1,nsig
             rlsg(i)=log(sigl(i))
             ENDDO
C-CRA             grid=rlsg
             DO i=1,nsig
             grid(i)=rlsg(i)
             ENDDO
         call gdcrdn(grid,nsig,rlsgdivt,nsigdivt)
         do k=1,nsig
          i0=grid(k)
          i0=max(1,min(nsigdivt-1,i0))
          i1=i0+1
          del=grid(k)-i0
          w0=1.-del
          w1=del
          do n=0,jcapdivt
           tdivt(n,k)=w0*cdivt(n,i0)+w1*cdivt(n,i1)
          end do
          do n=0,min(jcap,jcapdivt)
           bdivt(n,k)=tdivt(n,k)
          end do
          if(jcap.gt.jcapdivt) then
           do n=jcapdivt+1,jcap
            bdivt(n,k)=tdivt(jcapdivt,k)
           end do
          end if
          do n=1,jcap
           bdivt(n,k)=1./bdivt(n,k)
          end do
         end do
       return
       end
