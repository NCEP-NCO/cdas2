      subroutine wranlc(zc,dc,tc,qc,pc,rc,hourg,idateg,sigi,sigl,
     *  ioanl,jcap,nsig,on85,on85dt,lm2ml,factsml,factvml)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    wranlc     reorder and write sigma coefs.
c   prgmmr: parrish          org: w/nmc22    date: 90-10-10
c
c abstract: reorder from internal format and write sigma coefs.
c
c program history log:
c   90-10-10  parrish
c   08-04-04  ebisuzaki, f90 dynamic arrays, f90 loops
c
c   input argument list:
c     zc,dc,tc,qc,pr,rc - analysis sigma coefs for vort,div, etc
c     hourg    - analysis forecast hour
c     idateg   - initial date of analysis
c     sigi     - sigma values at interface of each sigma layer
c     sigl     - sigma values at mid-point of each sigma layer
c     ioanl    - unit number of analysis coefs
c     jcap     - triangular truncation
c     nsig     - number of sigma levels
c     on85     - on85 date record for guess coefs
c     on85dt   - on85 date record for data
c
c   output argument list:
c     no output arguments
c
c attributes:
c   language: f90
c   machine:  AIX
c
c$$$
c
          dimension zc((jcap+1)*(jcap+2),nsig)
          dimension dc((jcap+1)*(jcap+2),nsig)
          dimension tc((jcap+1)*(jcap+2),nsig)
          dimension qc((jcap+1)*(jcap+2),nsig)
          dimension pc((jcap+1)*(jcap+2))
          dimension rc((jcap+1)*(jcap+2))
          dimension idateg(4),sigi(nsig+1),sigl(nsig)
          dimension lm2ml((jcap+1)*(jcap+2))
          dimension factsml((jcap+1)*(jcap+2))
          dimension factvml((jcap+1)*(jcap+2))
          dimension dummy(201-nsig-1-nsig)
          dimension z((jcap+1)*(jcap+2))
 

c          integer iwash(2)
           character*4 on85(8),on85dt(8)
c          character*4 cwash(4)
c--------
c-------- local space
c--------
c     equivalence (cwash,iwash)
c     data iwash/X'00000000E6C1E2C8',X'C9D5C7E3D6D5C3E1'/

      integer iwash(4,4)
      data iwash/0,0,0,0, 230,193,226,200,
     1   201,213,199,227, 214,213,195,225/

c--------
c-------- set up index arrays for converting to output coefs
c--------
      nc=(jcap+1)*(jcap+2)
c--------
c-------- fix up stuff for record # 2
c--------
      dummy=0.
      waves=jcap
      xlayers=nsig
      trun=1.
      order=2.
      realform=1.
c-----------------------following 2 lines corrected on
c---------------------mark irdell request, 2-9/94  (dp)
      gencode=78.
      if(jcap.eq.62) gencode=80.
c--------
c-------- update on85 date word and idate, using date word from
c-------- data.
c--------
c     do i=1,4
c      on85(4+i)=cwash(i)
c     end do
      do i = 1, 4
          do j = 1, 4
             on85(4+i)(j:j) = char(iwash(j,i))
          enddo   
      enddo       

      print *,' on85 follows:'
      call prnon85(on85)
      print *,' on85dt follows:'
      call prnon85(on85dt)
	write(*,*) '>>w3fs03: ihour,iyear,month,iday=',
     1      ihour,iyear,month,iday

      call w3fs03(on85dt(3),ihour,iyear,month,iday,1)

	write(*,*) '<<w3fs03: ihour,iyear,month,iday=',
     1      ihour,iyear,month,iday
	write(*,*) 'idateg(4)=',idateg(4)

      idateg(1)=ihour
      idateg(2)=month
      idateg(3)=iday
c
c     iyear is 2 digit .. want 4 digit code
c        either same guess as initial guess or next year
      if (mod(idateg(4),100).ne.iyear) then
	  idateg(4) = idateg(4) + 1
      endif

      on85(3)=on85dt(3)
c--------
      rewind ioanl
c-------- hour,idate, etc.
      write(ioanl)on85
      print *,' on85 written to output coefs file follows:'
      call prnon85(on85)
      write(ioanl)hourg,idateg,sigi,sigl,dummy,waves,xlayers,trun,
     *   order,realform,gencode
c-------- terrain coefs
      do i=1,nc
       z(i)=factsml(i)*rc(lm2ml(i))
      end do
      write(ioanl)z
c-------- sfcp coefficients
      do i=1,nc
       z(i)=factsml(i)*pc(lm2ml(i))
      end do
      write(ioanl)z
c-------- temp coefficients
      do k=1,nsig
       do i=1,nc
        z(i)=factsml(i)*tc(lm2ml(i),k)
       end do
       write(ioanl)z
      end do
c------- div and vort
      do k=1,nsig
       do i=1,nc
        z(i)=factvml(i)*dc(lm2ml(i),k)
       end do
       write(ioanl)z
       do i=1,nc
        z(i)=factvml(i)*zc(lm2ml(i),k)
       end do
       write(ioanl)z
      end do
c-------- q coefs
      do k=1,nsig
       do i=1,nc
        z(i)=factsml(i)*qc(lm2ml(i),k)
       end do
       write(ioanl)z
      end do
      write(6,700)jcap,nsig,hourg,idateg
700   format(' some sigma coefficients written, jcap,nsig=',
     *   2i6,/,' hour,idate=',f10.1,4i5)
       rewind ioanl
       close (ioanl)
        write(*,*) '<<wrlanc'
      return
      end
