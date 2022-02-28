      subroutine rdfact(factor,idateg,hourg,nlath,nlon,isfc)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    rdfact  read and compute factor       
c   prgmmr: derber           org: w/nmc23    date: 92-09-08
c
c abstract: read in factor for use in reducing to 10m winds          
c
c program history log:
c   92-09-08  derber 
c
c   input argument list:
c     isfc     - unit number of bges file
c     idateg   - date time array for guess
c     hourg    - hour of guess
c     nlath    - number of latitudes on gaussian grid
c     nlon     - number of longitudes on gaussian grid
c
c   output argument list:
c     factor   - factor for reducing bottom sigma values to 10m
c
c attributes:
c   language: cft77
c   machine:  cray ymp
c
c$$$
c
C-CRA          dimension fldr (nlon,2*nlath-2)
C-CRA          real factor(2*nlath+1,nlon+2)
C-CRA          integer idateg(4)
C-CRA          integer lab85(8),idate(4)          
C-CRA          integer idate5(5)
 
          dimension fldr (192,2*48-2)
          real factor(2*48+1,192+2)
          integer idateg(4)
          integer lab85(8),idate(4)          
          integer idate5(5)
c--------
c-------- scratch space
c--------
c
c
c   read surface file to get 10m winds
c
C-CRA                factor=1.   
c       real factor(2*nlath+1,nlon+2)
          DO ITMP=1,(2*nlath+1)*(nlon+2)
          factor(ITMP,1)=1.   
          ENDDO
      print *, 'calculating factor'
      rewind isfc
      read (isfc,end=7781,err=7781) lab85
      read (isfc,end=7781,err=7781) fhour,idate
      idate5(1)=idateg(4)
      idate5(2)=idateg(2)
      idate5(3)=idateg(3)
      idate5(4)=idateg(1)
      idate5(5)=0
      call w3fs21(idate5,nming)
      nming=nming+60*hourg
      idate5(1)=idate(4)
      idate5(2)=idate(2)
      idate5(3)=idate(3)
      idate5(4)=idate(1)
      idate5(5)=0
      call w3fs21(idate5,nmins)
      nmins=nmins+60*fhour
      print 101,fhour,idate
  101 format(' fhour=',f5.0,' idate=',4i5)
      print *,' for bges file, nmins=',nmins
      print *,' for ges file, nming=',nming
      if(nmins.ne.nming) go to 7781
c
      do 102 j = 1,13
c
        read(isfc,end=7781,err=7781)
  102 continue
      read(isfc,end=7781,err=7781) fldr
      do 1101 j=1,nlon
      do 1101 i = 1, 2*nlath-2
        factor(i+1,j)=fldr(j,2*nlath-1-i)
 1101 continue
c
      sumn=0.
      sums=0.
      do 780 j=1,nlon
        sumn=factor(2*nlath-1,j) +sumn
        sums=factor(2,j) +sums
 780  continue
      sumn=sumn/nlon
      sums=sums/nlon
      do 781 j=1,nlon
        factor(2*nlath,j)=sumn
        factor(1,j)=sums
 781  continue
      xmax=-1.e20
      xmin=1.e20
      do 1102 j=1,nlon
      do 1102 i = 1, 2*nlath-2
        xmax=max(factor(i,j),xmax)
        xmin=min(factor(i,j),xmin)
 1102 continue
      print *, 'factor  max and min = ',xmax,xmin
7781  continue
      close(isfc)
      return    
      end
