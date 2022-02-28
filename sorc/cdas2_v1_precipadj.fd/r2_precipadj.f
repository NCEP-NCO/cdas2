	program adjbges
c
c	adjust bges
c
c	take observed pentad precip - analyzed precip 
c	and use it to correct surface moisture
c
c
c	unit 10 = input bges
c	unit 11 = input observed precip (same grid) (prev pentad)
c	unit 12 = input model precip (same grid) (prev pentad)
c		  (optional) land mask
c		  (optional) runoff
c	unit 51  = output bges
c	unit 52  = output grib file
c
c	assume that bges and observed precip are correct
c	use idate(4) from bges file to get date of month
c
c v1.3 .. used by R2 2001-2004 (over the counter)
c
c v1.4 5/2004 for NCO operations
c	read decoded grib files in real*4 rather than real*8
c		to eliminate double precision version of wgrib
c	do pentad averaging .. eliminate pentad averging program.
c
c v1.5 12/2012 for wcoss .. some f95 notation Wesley Ebisuzaki

	parameter (lonf= 192 )
	parameter (latg= 94 )
	parameter (lsoil= 2 )

c	parameter (lonf=192)
c	parameter (latg=94)
c	parameter (lsoil=2)

	parameter (ijdim= lonf * latg)
	parameter (lugi=10,lugi1=11,lugi2=12)
	parameter (lugo1=51)
	parameter (lugo2=52)
	parameter (nhour= 6 )	! analysis cycle time
	parameter (undef=9.999e20) ! wgrib undefined number
	parameter (undeflow=9.998e20) ! wgrib undefined number

c	rlev1 = 0 .. 10 cm below ground
c	rlev2 = 10 .. 200 cm below ground
	parameter(rlev1 = 112 * 256 * 256 + 0 * 256 + 10)
	parameter(rlev2 = 112 * 256 * 256 + 10 * 256 + 200)

	character*32 label
	real fhour ! forecast hour
	integer idate(4) ! initial date code
	real sfctmp(ijdim) ! sfc temp
	real soilm(ijdim,lsoil) ! soil moisture
	real snow(ijdim) ! sfc temp
	real temp(ijdim), tempv(ijdim,lsoil) ! temp array

	real obspre(ijdim) ! observed precip
	real mdlpre(ijdim) ! model precip
	real offfset, bias(ijdim,lsoil) ! correction term in cm
	real land(ijdim) ! land mask, 1=land 0=water
	real runoff(ijdim) ! runoff
	real*4 tmp1(ijdim), tmp2(ijdim), tmp3(ijdim)
	real soilcap(lsoil)   ! soil moisture capacity for each layer (cm)
	integer year, month, day, hour
	real factor, rtmp(ijdim)
	logical dorun, flag1(ijdim), flag2(ijdim), flag3(ijdim)

	data soilcap /100.,1900./  ! layer depth in mm

	dorun = .true.

c	initialize bias
	bias = 0.0

c	read observed precip

c       read real*4, convert to real*8
	read(lugi1,end=100,err=100) tmp1
	obspre = tmp1
	goto 110
100	continue
c	    no obs precip
            write(*,*) "WARNING!!!!"
            write(*,*) "Obs precip read error or end of record"
            write(*,*) "ZERO SOIL WETNESS CORRECTION"
	    do i = 1, ijdim
		obspre(i) = undef
	    enddo
110     continue

	call maxmin(obspre,ijdim,rmin,rmax)
	write(*,*) 'obspre min/max (mm/s) =', rmin, rmax
c
c	read model precip
c
	mdlpre = 0.0
	land = 0.0
	runoff = 0.0
	flag1 = .true.
	flag2 = .true.
	flag3 = .true.

	n = 0
120	continue
	    read(lugi2,end=140,err=130) tmp1
	    read(lugi2,end=130,err=130) tmp2
	    read(lugi2,end=130,err=130) tmp3

	    rtmp = tmp1	
	    call maxmin(rtmp,ijdim,rmin,rmax)
	    write(*,*) 'n mdlpre min/max=', n, rmin, rmax

	    do i = 1, ijdim
		if (tmp1(i).gt.undeflow) flag1(i) = .false.
		if (tmp2(i).gt.undeflow) flag2(i) = .false.
		if (tmp3(i).gt.undeflow) flag3(i) = .false.
		mdlpre(i) = mdlpre(i) + tmp1(i)
		land(i) = land(i) + tmp2(i)
		runoff(i) = runoff(i) + tmp3(i)
	    enddo
	    n = n + 1
	goto 120
130	continue
	    write(*,*) 'WARNING!!!!'
	    write(*,*) 'READ ERROR - MODEL precip - precipadj'
	    write(*,*) 'WARNING!!!!'
	    n = 0
140	continue
	if (n.ne.20.and.n.ne.24) then
            write(*,*) "WARNING!!!!"
            write(*,*) "Bad model precip file"
            write(*,*) "running without soil moisture adjustment"
	    mdlpre = undef
	    land = 1.0
	    runoff = 0.0
	    n = 0
	    dorun = .false.
	    stop 8
	else
c	    Average model parameters
	    factor = 1.0 / n
	    do i = 1, ijdim
		if (flag1(i)) then
		    mdlpre(i) = factor*mdlpre(i)
		else
		    mdlpre(i) = undef
		endif
		if (flag2(i)) then
		    land(i) = factor*land(i)
		else
		    land(i) = 1.0
		endif
		if (flag3(i)) then
		    runoff(i) = factor*runoff(i)
		else
		    runoff(i) = 0.0
		    land(i) = 0.0
		endif
	    enddo
	endif
	
c	fixup values, runoff in mm (per 6 hours)
	if (dorun) then
	   do i = 1, ijdim
	    if (abs((runoff(i)-undef)/undef).lt.1e-5) runoff(i) = 0.0
	   enddo
	endif

	close(lugi1)
	close(lugi2)
	write(*,*) 'obspre (pnt 1 100):',obspre(1), obspre(100)
	write(*,*) 'mdlpre:(pnt i 100)',mdlpre(1), mdlpre(100)
	write(*,*) 'runoff is ',dorun
	if (dorun) write(*,*) 'runoff:',runoff(1), runoff(100)

	call maxmin(obspre,ijdim,rmin,rmax)
	write(*,*) 'obspre min/max=', rmin, rmax
	call maxmin(mdlpre,ijdim,rmin,rmax)
	write(*,*) 'mdlpre min/max=', rmin, rmax
	call maxmin(runoff,ijdim,rmin,rmax)
	write(*,*) 'runoff min/max=', rmin, rmax
	call maxmin(land,ijdim,rmin,rmax)
	write(*,*) 'land min/max=', rmin, rmax

	read(lugi) label
	read(lugi) fhour, idate
	write(lugo1) label
	write(lugo1) fhour,idate
	write(6,*) 'bguess date ',fhour,' yr=',idate(4),' mo=',
     1		idate(2),' day=',idate(3),' hr=',idate(1)

	year = idate(4)
c	Y2K, y2k problem
	if (year.lt.200) year = year + 1900
	month = idate(2)
	day = idate(3)
	hour = idate(1)

	call pdays(year,month,day,factor)

c	PRATE = kg/m/m/s
c	x 86400 to get to mm/day
c	x factor which is a correction for leap year
c	x (nhour/24) to get correction for a n-hour cycle
c	net result .. correction in mm

	write(*,*) ' time factor= ',factor, ' should be approx 1'
	factor = 86400.0 * factor * (nhour/24.0)
	write(*,*) ' scaling factor= ',factor

	read(lugi) sfctmp
	read(lugi) soilm
	read(lugi) snow

	call maxmin(sfctmp,ijdim,rmin,rmax)
	write(*,*) 'sfctmp min/max=', rmin, rmax
	call maxmin(snow,ijdim,rmin,rmax)
	write(*,*) 'snow min/max=', rmin, rmax
	call maxmin(soilm(1,1),ijdim,rmin,rmax)
	write(*,*) 'soilm(:,1) min/max=', rmin, rmax
	call maxmin(soilm(1,2),ijdim,rmin,rmax)
	write(*,*) 'soilm(:,2) min/max=', rmin, rmax

c	now to do the correction

c
c	v1.4
c
c	case1:
c	   obspre .ge. mdlpre-runoff .and. runoff .gt. 0
c		=> assume extra precip would have gone into runoff 
c		=> do nothing
c	   obspre .lt. mdlpre-runoff .and. runoff .gt. 0
c		=> apply precip correction
c	   obspre .ne. mdlpre .and. runoff .eq. 0
c		=> assume no runoff either way
c		=> apply precip correction
c
c	   surface temp < 0C .. do nothing .. frozen ground
c
c	note:
c		runoff == 0 => v1.0 scheme
c		temperature correction
c		correction is applied to top layer first
c
c	write diagnostic file
c		soilm(ijdim,lsoil)
c
cif (lsoil.eq.2) then
c    write(*,*) '>>writing soilm'
c    call setctr(7, 1, 80)
c    call ezgbw(soilm(1,1),lonf,latg,.false.,.false.,
c    2          144,.false.,rlev1,year,month,day,hour,lugo2,
c    3		0,0,'hour','uninit anl')
c    call ezgbw(soilm(1,2),lonf,latg,.false.,.false.,
c    2          144,.false.,rlev2,year,month,day,hour,lugo2,
c    3		0,0,'hour','uninit anl')
c    write(*,*) '<<writing soilm'
cendif


	do i = 1, ijdim

c	    check for water point
	    if (land(i).eq.0) goto 200

c	    check for undefined values of precip
	    if (abs((obspre(i)-undef)/undef).lt.1e-5) goto 200
c	    should be defined but check anyways
	    if (abs((mdlpre(i)-undef)/undef).lt.1e-5) goto 200

c	    check for surface temperature
	    if (sfctmp(i).lt.273.16) goto 200

c	    figure out the precip surplus
	    if (runoff(i).gt.0.0) then
		offset = (obspre(i) - (mdlpre(i)-runoff(i)))*factor
		if (offset.gt.0.0) offset = 0.0
	    else
	        offset = (obspre(i) - mdlpre(i))*factor
	    endif

	    if (offset.eq.0.0) goto 200
c
c	    do k = lsoil, 1, -1  => bottom to top
c	    do k = 1, lsoil      => top to bottom

	    do k = 1, lsoil
		smoist = soilm(i,k)*soilcap(k) + offset
		if (smoist.lt.0.0) then
		    bias(i,k) = offset - smoist
		    offset = smoist
		    soilm(i,k) = 0.0
		else if (smoist.gt.soilcap(k)) then
		    bias(i,k) = offset - (smoist - soilcap(k))
		    offset = smoist - soilcap(k)
		    soilm(i,k) = 1.0
		else
		    bias(i,k) = offset
		    soilm(i,k) = smoist / soilcap(k)
		    goto 200
		endif
	    enddo
200	    continue
	enddo

c
c	write diagnostic file
c		soilm(ijdim,lsoil)
c		bias(ijdim,lsoil)
c
	if (lsoil.eq.2) then
c    write(*,*) '>>new diag'
c    call ezgbw(soilm(1,1),lonf,latg,.false.,.false.,
c    2          144,.false.,rlev1,year,month,day,hour,lugo2,
c    3		0,0,'hour','init anl')
c    call ezgbw(soilm(1,2),lonf,latg,.false.,.false.,
c    2          144,.false.,rlev2,year,month,day,hour,lugo2,
c    3		0,0,'hour','init anl')
c
	    call ezgbw1(bias(1,1),lonf,latg,234,.false.,rlev1,
     2		year,month,day,hour,lugo2)
	    call ezgbw1(bias(1,2),lonf,latg,234,.false.,rlev2,
     2		year,month,day,hour,lugo2)

	    call maxmin(bias(1,1),ijdim,rmin,rmax)
  	    write(*,*) 'bias(:,1) min/max=', rmin, rmax
	    call maxmin(bias(1,2),ijdim,rmin,rmax)
  	    write(*,*) 'bias(:,2) min/max=', rmin, rmax

	    call maxmin(soilm(1,1),ijdim,rmin,rmax)
            write(*,*) 'new soilm(:,1) min/max=', rmin, rmax

           call maxmin(soilm(1,2),ijdim,rmin,rmax)
            write(*,*) 'new soilm(:,2) min/max=', rmin, rmax

	endif

	write(lugo1) sfctmp
	write(lugo1) soilm
	write(lugo1) snow

c	copy rest of the file

	read(lugi) tempv
	write(lugo1) tempv

1000	read(lugi,end=1200) temp
	    write(lugo1) temp
	    goto 1000
1200	continue
	stop
	end

	subroutine maxmin(a,idim,rmin,rmax)

	parameter (undef=9.999e20) ! wgrib undefined number
	parameter (undeflow=9.998e20) ! wgrib undefined number

	real a(idim), rmax, rmin
	rmax=-1e20
	rmin= 1e20
	do i = 1, idim
	    if (a(i).lt.undeflow) then
		if (rmin.gt.a(i)) rmin = a(i)
		if (rmax.lt.a(i)) rmax = a(i)
	    endif
	enddo
	return
	end
	subroutine pdays(year,month,day,factor)
c
c	returns the scaling factor for 6 vs 5 days
c	factor = (no. of days in preceeding pentad) /
c		no of days in current pentad
c
	integer year,month,day
	real factor

	factor = 1.0
	if (month.eq.1 .or. month.ge.4) return
	
	if (mod(year,4).ne.0) return
	if (mod(year,100).eq.0 .and. mod(year,400).ne.0) return

c	must be leap year
	if (month.eq.2) then
	   if (day.ge.25) then
		factor = 5.0 / 6.0
		return
	   endif
	   return
	endif
	if (day.eq.1) then
	    factor = 5.0 / 6.0
	    return
	endif
	if (day.le.6) then
	    factor = 6.0 / 5.0
	    return
	endif
	return
	end

c
c	added variance as a type
c
c=================================================================
c
c	change center/model id (in common block)
c
	subroutine setctr(ictr, isctr, imodel)
c
	integer ictr, isctr, imodel
	integer center, subcen, model
	common /gribc/center, subcen, model
	center = ictr
	subcen = isctr
	model = imodel
	return
	end
c=================================================================
c
	block data gribitc
	integer center, subcen, model
	common /gribc/center, subcen, model
c	nmc
	data center/7/
c	reanalysis
	data subcen/3/
c	who knows
	data model/195/
	end
c
c=================================================================
	subroutine ezgbw1(data,nx,ny,iparm,islola,level,
     2	iyear,imonth,iday,ihour,iunit)
c
c
c	ez grib write 1			W. Ebisuzaki
c
c	sets many defaults - so you don't have to 
c		lola - grid,  gaussian (future)
c	writes only 1 field that is for an instantaneous time
c
c
c	real data(nx,ny)	- data to written as grib record
c	integer iparm		- grib parameter number ex. 33 -> U wind
c	logical islola		- .true. -> lola  .false. -> gaussian grid
c	real level		- 0.0 -> 1.0 (sigma), 1.+ ->2000 (prs)
c				  -1 -> surface
c
c
	real data(nx,ny), level
	integer iparm,iyear,imonth,iday,ihour,iunit
	logical islola, ldummy(1)

	call ezgbw(data,nx,ny,ldummy,.false.,iparm,islola,level,
     2	iyear,imonth,iday,ihour,iunit,0,0,'hour','inst')

	return
	end

c=================================================================
	subroutine ezgbw2(data,nx,ny,iparm,islola,level,
     2	iyear,imonth,iday,ihour,iunit,itime1,itime2,timeun,type)
c
c	ez grib write 2			W. Ebisuzaki
c
c	sets many defaults - so you don't have to 
c		lola - grid,  gaussian (future)
c	writes only 1 field that is for a time average/accumulation/inst
c
c
	real data(nx,ny), level
	integer iparm,iyear,imonth,iday,ihour,iunit
	logical islola, ldummy(1)
	integer itime1,itime2
	character*(*) timeun,type

	call ezgbw(data,nx,ny,ldummy,.false., iparm,islola,level,
     2	iyear,imonth,iday,ihour,iunit,itime1,itime2,timeun,type)
	return
	end
c=================================================================
	subroutine ezgbw(data,nx,ny,bit,hasbit,iparm,islola,level,
     2	iyear,imonth,iday,ihour,iunit,itime1,itime2,timeun,type)
c
c
c	ez grib write 				W. Ebisuzaki
c
c	sets many defaults - so you don't have to 
c		lola - grid,  gaussian (future)
c
c	real data(nx,ny)	- data to written as grib record
c	logical bit(*), hasbit  - bitmap of defined values
c	integer iparm		- grib parameter number ex. 33 -> U wind
c	logical islola		- .true. -> lola  .false. -> gaussian grid
c	real level		- 0.0 -> 1.0 (sigma), 1.+ ->2000 (prs)
c				  -1 -> surface
c	integer itime1, itime2	- time average from time1 to time2
c	character*(*) timeun	- 'hour', 'day', 'month'
c	character*(*) type	- 'ave', 'acc', 'inst', 'climo', 'var'
c				- 'init anl', 'uninit anl'
c note: type='climo'
c	itime1 = length of mean
c	itime2 = number of years of data
c
	parameter (ilpds=28)
	parameter (mxbit=16)
	parameter (iptv=132)
	parameter (mxnxny=384*190)
	parameter (ngrib = 100 + ilpds + mxnxny*(mxbit+1)/8+10)

	real data(nx,ny), grib(ngrib), level
	integer iyear, imonth, iday, ihour
	character*(*) timeun, type
	logical hasbit, bit(*)

	integer levelt, level1, level2, itime1, itime2
	integer iftu, time1, time2, timerg, nave, nmiss
	integer dscale(255)

	integer gridty, lengrb
	logical init, islola

        integer center, subcen, model
        common /gribc/center, subcen, model

	data init/.false./

	save init, dscale

	if (nx*ny.gt.mxnxny) then
	    write(*,*) '**************'
	    write(*,*) '*** ezgbwr ***'
	    write(*,*) 'array too small: edit, recompile'
	    write(*,*) '*** Fatal Error ***'
	    call exit(99)
	    stop
	endif

	if (.not.init) then
	    call idsdef(iptv,dscale)
	    init = .true.
	endif

	nave = 0
	nmiss = 0

c	gridtyp 0 -> lola   4 -> gaussian grid
	gridty = 0
	if (.not.islola) gridty = 4

c	set colatitude
	if (.not.islola) then
	    if (ny.eq.94) then
c		T62 Gaussian lat
	        colat1 = (90.0-88.542)*3.141592654/180.0
	    elseif (ny.eq.190) then
c		T126 Gaussian lat
	        colat1 = (90.0-89.227)*3.141592654/180.0
	    else
		write(*,*) '*** Missing Gaussian colat in gribit ***'
		write(*,*) '*** please add ***'
		write(*,*) '*** 0 used ***'
	        colat1 = 0.0
	    endif
	endif

c	set level
	level1 = 0
	if (level.gt.256*256) then
	    levelt = level / (256*256)
	    level1 = mod(int(level)/256, 256)
	    level2 = mod(int(level), 256)
	elseif (level.ge.2000.0) then
c	    meters above ground - 2000
	    level2 = level - 2000.0
	    levelt = 105
	elseif (level.gt.1.0) then
c	    must be pressure level
	    level2 = level
	    levelt = 100
	elseif (level.ge.0.0 .and. level.le.1.0) then
c	    must be sigma level
	    level2 = 10000.0 * level
	    levelt = 107
	elseif (level.eq.-1.0) then
c	    must be the surface
	    level2 = 0
	    levelt = 1
	elseif (level.gt.-256.0.and.level.lt.-1.0) then
	    rlevel = -level
	    levelt = int(rlevel)

            if (levelt.le.100.or.levelt.eq.102.or.levelt.eq.103 .or.
     1            levelt.eq.105.or.levelt.eq.107.or.levelt.eq.109 .or.
     1            levelt.eq.111.or.levelt.eq.113.or.levelt.eq.125 .or.
     1            levelt.eq.160.or.levelt.eq.200.or.levelt.eq.201) then
                level1 = 0
	        level2 = int(1000*(rlevel - int(rlevel)) + 0.5)
	    else
	        level2 = int(1000*(rlevel - int(rlevel)) + 0.5)
	        level1 = level2 / 256
	        level2 = mod(level2,256)
            endif

	else if (level.eq.-257.0) then
c	    sigma 0..1
	    levelt = 108
	    level2 = 100
	else
c	    ????
	    write(*,*) 'ezgbw2: unknown level?', level
	    call exit(99)
	    stop
	endif

c	for an average/accum.
	if (timeun.eq.'hour'.or.timeun.eq.'HOUR') then
	    iftu = 1
	elseif (timeun.eq.'day'.or.timeun.eq.'DAY') then
	    iftu = 2
	elseif (timeun.eq.'month'.or.timeun.eq.'MONTH') then
	    iftu = 3
	elseif (timeun.eq.'year'.or.timeun.eq.'YEAR') then
	    iftu = 4
	    write(*,*) 'ezwgb: check year parameter!!'
	    stop
	else
	    write(*,*) 'ezwgb bad timeun:',timeun
	    call exit(99)
	endif

	write(*,*) 'type = ',type
	if (type.eq.'init anl' .or. type.eq.'INIT ANL') then
	    timerg = 1
	    time1 = 0
	    time2 = 0
	else if (type.eq.'uninit anl' .or. type.eq.'UNINIT ANL') then
	    timerg = 0
	    time1 = 0
	    time2 = 0
	else if (type.eq.'inst'.or.type.eq.'INST'.or.itime1.eq.
     1			itime2) then
	    timerg = 10
	    time1 = 0
	    time2 = min0(itime1, itime2)
	else if (type.eq.'acc'.or.type.eq.'ACC') then
	    timerg = 4
	    time1 = min0(itime1, itime2)
	    time2 = max0(itime1, itime2)
	else if (type.eq.'ave'.or.type.eq.'AVE') then
	    timerg = 3
	    time1 = min0(itime1, itime2)
	    time2 = max0(itime1, itime2)
	else if (type.eq.'climo'.or.type.eq.'CLIMO') then
c	    time1 = length of average
c	    time2 = number of years of data
	    timerg = 51
	    time1 = 0
	    time2 = itime1
	    nave = itime2
	else if (type.eq.'var'.or.type.eq.'VAR') then
	    timerg = 118
	    time1 = min0(itime1, itime2)
	    time2 = max0(itime1, itime2)
	else
	     write(*,*) 'ezwgb2 type:',type
	     call exit(99)
	endif

	ids = dscale(iparm)
	if (type.eq.'var'.or.type.eq.'VAR') ids = 2*ids

c	bitmap?
	ibms = 0
	if (hasbit) ibms = 1

	call gribit(data,bit,gridty,nx,ny,mxbit,colat1,ilpds,
     1		iptv,center,subcen,model,ibms,iparm,
     2		levelt, level1, level2,
     3		iyear, imonth, iday, ihour, 
     4		iftu, time1, time2, timerg,
     5		nave, nmiss, ids, grib, lengrb, ierr)

	if (ierr.eq.0) call binwrt(iunit, grib, lengrb)

	return
	end
c
c=====================================================================
c
	subroutine ezgbwk(data,lbm,kpds,kgds,iunit)
c
c	write grib record, kpds, kgds must be defined by caller
c
	parameter (mxbit=16)
	parameter (ilpds = 28)
	parameter (mxnxny=384*190)
	parameter (ngrib = 100 + ilpds + mxnxny*(mxbit+1)/8+10)

	real data(*)
	logical lbm(*)
	integer kpds(*), kgds(*)
	integer center, subcen, model
	real grib(ngrib)

c	data(nx,ny) -- data to write

c	lbm(nx,ny) -- optional bitmap

c	lola/gaussian
	idrt = kgds(1)

c	dim of data(nx,ny)
	nx = kgds(2)
	ny = kgds(3)

c	maximum number of bits

c	colatitude
c	colat1 = (90e3 - (kpds(6)))*acos(-1.0)/(180e3)
	colat1 = (90e3 - (kgds(4)))*acos(-1.0)/(180e3)

c	length of pds (usually 28)
c	ilpds

c	parameter table (usually 1)
	iptv = kpds(19)
c	forecast center
	center = kpds(1)
c	subcenter
	subcen = kpds(23)
c	generating model
	model = kpds(2)
c	has bitmap?
	ibms = mod(kpds(4)/64,2)
	write(*,*) 'ezgbwk1 ibms=',ibms

c	input parameter
	ipu = kpds(5)

c	level indicator
	itl=kpds(6)
c	old code
c        if(itl.LE.100.OR.itl.EQ.102.OR.itl.EQ.103.OR.
c     &    itl.EQ.105.OR.itl.EQ.107.OR.itl.EQ.111.OR.
c     &    itl.EQ.160)  then
c	mod 2/28/95
	if (itl.le.100 .or. itl.eq.102 .or. itl.eq.103 .or.
     1		itl.eq.105 .or. itl.eq.107 .or. itl.eq.109 .or.
     1		itl.eq.111 .or. itl.eq.113 .or. itl.eq.125 .or.
     1		itl.eq.160 .or. itl.eq.200 .or. itl.eq.201) then

            il1=0
            il2=kpds(7)
        else
            il1=mod(kpds(7)/256,256)
            il2=mod(kpds(7),256)
        endif

c	date
c	iyr = kpds(8)
	iyr = kpds(8) + (kpds(21)-1)*100

	imo = kpds(9)
	idy = kpds(10)
	ihr = kpds(11)
c	forecast unit, and times, and units
	iftu = kpds(13)
	ip1 = kpds(14)
	ip2 = kpds(15)
	itr = kpds(16)

c	number in/missing from average
	ina = kpds(17)
	inm = kpds(20)
c	integer decimal scaling
	ids = kpds(22)

c	write(*,*) 'idrt=',idrt,' nx/ny=',nx,ny,' mxbit=',mxbit
c	write(*,*) 'colat1=',colat1,' ilpds=',ilpds,' iptv=',iptv
c	write(*,*) 'icen=',icen,' igen=',igen,' ibms=',ibms,' param=',ipu
c	write(*,*) 'times:',itl,il1,il2
c	write(*,*) 'date:',iyr,imo,idy,ihr
c	write(*,*) 'forcast times:',iftu,ip1,ip2,itr
c	write(*,*) 'no. for ave=',ina,' missing=',inm,' dec scale=',ids
	call gribit(data,lbm,idrt,nx,ny,mxbit,colat1,ilpds,iptv,
     a      center,subcen,model,
     1	    ibms,ipu,itl,il1,il2,iyr,imo,idy,ihr,iftu,ip1,ip2,itr,
     2	    ina,inm,ids,grib,lengrb,ierr)

	if (ierr.eq.0) call binwrt(iunit, grib, lengrb)

	return
	end
c=====================================================================
	subroutine binwrt(iunit,string,n)
c
c	binary write
c
	character*1 string(n)
	write(iunit) string
	return
	end
c=====================================================================
CFPP$ NOCONCUR R
      SUBROUTINE GRIBIT(F,LBM,IDRT,IM,JM,MXBIT,COLAT1,
     &                  ILPDS,IPTV,ICEN,ISUBCN,IGEN,IBMS,IPU,ITL,
     &                  IL1,IL2,
     &                  IYR,IMO,IDY,IHR,IFTU,IP1,IP2,ITR,INA,INM,IDS,
     &                  GRIB,LGRIB,IERR)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    GRIBIT      CREATE GRIB MESSAGE
C   PRGMMR: IREDELL          ORG: W/NMC23    DATE: 92-10-31
C
C ABSTRACT: CREATE A GRIB MESSAGE FROM A FULL FIELD.
C   AT PRESENT, ONLY GLOBAL LATLON GRIDS AND GAUSSIAN GRIDS ARE ALLOWED.
C
C PROGRAM HISTORY LOG:
C   92-10-31  IREDELL
C
C USAGE:    CALL GRIBIT(F,LBM,IDRT,IM,JM,MXBIT,COLAT1,
C    &                  ILPDS,IPTV,ICEN,IGEN,IBMS,IPU,ITL,IL1,IL2,
C    &                  IYR,IMO,IDY,IHR,IFTU,IP1,IP2,ITR,INA,INM,IDS,
C    &                  GRIB,LGRIB,IERR)
C   INPUT ARGUMENT LIST:
C     F        - REAL (IM*JM) FIELD DATA TO PACK INTO GRIB MESSAGE
C     LBM      - LOGICAL (IM*JM) BITMAP TO USE IF IBMS=1
C     IDRT     - INTEGER DATA REPRESENTATION TYPE
C                (0 FOR LATLON OR 4 FOR GAUSSIAN)
C     IM       - INTEGER LONGITUDINAL DIMENSION
C     JM       - INTEGER LATITUDINAL DIMENSION
C     MXBIT    - INTEGER MAXIMUM NUMBER OF BITS TO USE (0 FOR NO LIMIT)
C     COLAT1   - REAL FIRST COLATITUDE OF GAUSSIAN GRID IF IDRT=4
C     ILPDS    - INTEGER LENGTH OF THE PDS (USUALLY 28)
C     IPTV     - INTEGER PARAMETER TABLE VERSION (USUALLY 1)
C     ICEN     - INTEGER FORECAST CENTER (USUALLY 7)
C     ISUBCN   - INTEGER SUBCENTER
C     IGEN     - INTEGER MODEL GENERATING CODE
C     IBMS     - INTEGER BITMAP FLAG (0 FOR NO BITMAP)
C     IPU      - INTEGER PARAMETER AND UNIT INDICATOR
C     ITL      - INTEGER TYPE OF LEVEL INDICATOR
C     IL1      - INTEGER FIRST LEVEL VALUE (0 FOR SINGLE LEVEL)
C     IL2      - INTEGER SECOND LEVEL VALUE
C    &                  IYR,IMO,IDY,IHR,IFTU,IP1,IP2,ITR,INA,INM,IDS,
C    &                  GRIB,LGRIB,IERR)
C     IYR      - INTEGER YEAR (YYYY)
C     IMO      - INTEGER MONTH
C     IDY      - INTEGER DAY
C     IHR      - INTEGER HOUR
C     IFTU     - INTEGER FORECAST TIME UNIT (1 FOR HOUR)
C     IP1      - INTEGER FIRST TIME PERIOD
C     IP2      - INTEGER SECOND TIME PERIOD (0 FOR SINGLE PERIOD)
C     ITR      - INTEGER TIME RANGE INDICATOR (10 FOR SINGLE PERIOD)
C     INA      - INTEGER NUMBER INCLUDED IN AVERAGE
C     INM      - INTEGER NUMBER MISSING FROM AVERAGE
C     IDS      - INTEGER DECIMAL SCALING
C
C   OUTPUT ARGUMENT LIST:
C     GRIB     - CHARACTER (LGRIB) GRIB MESSAGE
C     LGRIB    - INTEGER LENGTH OF GRIB MESSAGE
C                (NO MORE THAN 100+ILPDS+IM*JM*(MXBIT+1)/8)
C     IERR     - INTEGER ERROR CODE (0 FOR SUCCESS)
C
C SUBPROGRAMS CALLED:
C   GTBITS     - COMPUTE NUMBER OF BITS AND ROUND DATA APPROPRIATELY
C   W3FI72     - ENGRIB DATA INTO A GRIB1 MESSAGE
C
C ATTRIBUTES:
C   LANGUAGE: CRAY FORTRAN
C
C$$$
      REAL F(IM*JM)
      LOGICAL LBM(IM*JM)
      CHARACTER GRIB(*)
      INTEGER IPDS(25),IGDS(18),IBDS(9)
C
      PARAMETER(JLPDS=28)
C     INTEGER IBM(IM*JM*IBMS+1-IBMS)
C     REAL FR(IM*JM)
      INTEGER IBM(192 * 94)
      REAL FR(192 * 94)
      CHARACTER PDS(JLPDS)
C
c	write(*,*) '>>gribit'
      NF=IM*JM
      LONI=NINT(360.E3/IM)
      IF(IDRT.EQ.0) THEN
        LAT1=NINT(90.E3)
        LATI=NINT(180.E3/(JM-1))
        IF(IM.EQ.144.AND.JM.EQ.73) THEN
          IGRID=2
        ELSEIF(IM.EQ.360.AND.JM.EQ.181) THEN
          IGRID=3
        ELSE
          IGRID=255
        ENDIF
      ELSEIF(IDRT.EQ.4) THEN
        LAT1=NINT(90.E3-180.E3/ACOS(-1.)*COLAT1)
        LATI=JM/2
        IF(IM.EQ.192.AND.JM.EQ.94) THEN
          IGRID=98
        ELSE
          IGRID=255
        ENDIF
      ELSE
        IERR=40
        RETURN
      ENDIF
      IPDS(01)=ILPDS    ! LENGTH OF PDS
      IPDS(02)=IPTV     ! PARAMETER TABLE VERSION ID
      IPDS(03)=ICEN     ! CENTER ID
      IPDS(04)=IGEN     ! GENERATING MODEL ID
      IPDS(05)=IGRID    ! GRID ID
      IPDS(06)=1        ! GDS FLAG
      IPDS(07)=IBMS     ! BMS FLAG
      IPDS(08)=IPU      ! PARAMETER UNIT ID
      IPDS(09)=ITL      ! TYPE OF LEVEL ID
      IPDS(10)=IL1      ! LEVEL 1 OR 0
      IPDS(11)=IL2      ! LEVEL 2
c     IPDS(12)=IYR      ! YEAR
      IPDS(12)=mod((IYR-1),100) + 1
      IPDS(13)=IMO      ! MONTH
      IPDS(14)=IDY      ! DAY
      IPDS(15)=IHR      ! HOUR
      IPDS(16)=0        ! MINUTE
      IPDS(17)=IFTU     ! FORECAST TIME UNIT ID
      IPDS(18)=IP1      ! TIME PERIOD 1
      IPDS(19)=IP2      ! TIME PERIOD 2 OR 0
      IPDS(20)=ITR      ! TIME RANGE INDICATOR
      IPDS(21)=INA      ! NUMBER IN AVERAGE
      IPDS(22)=INM      ! NUMBER MISSING
c     IPDS(23)=20       ! CENTURY
      IPDS(23)=((IYR-IPDS(12))/100) + 1
      IPDS(24)=ISUBCN   ! SUBCENTER
      IPDS(25)=IDS      ! DECIMAL SCALING
      IGDS(01)=0        ! NUMBER OF VERTICAL COORDS
      IGDS(02)=255      ! VERTICAL COORD FLAG
      IGDS(03)=IDRT     ! DATA REPRESENTATION TYPE
      IGDS(04)=IM       ! EAST-WEST POINTS
      IGDS(05)=JM       ! NORTH-SOUTH POINTS
      IGDS(06)=LAT1     ! LATITUDE OF ORIGIN
      IGDS(07)=0        ! LONGITUDE OF ORIGIN
      IGDS(08)=128      ! RESOLUTION FLAG
      IGDS(09)=-LAT1    ! LATITUDE OF END
      IGDS(10)=-LONI    ! LONGITUDE OF END
      IGDS(11)=LATI     ! LAT INCREMENT OR GAUSSIAN LATS
      IGDS(12)=LONI     ! LONGITUDE INCREMENT
      IGDS(13)=0        ! SCANNING MODE FLAGS
      DO I=14,18
        IGDS(I)=0     ! NOT USED
      ENDDO
      DO I=1,9
      IBDS(I)=0       ! BDS FLAGS
      ENDDO
      NBM=NF
c	write(*,*) 'gribit:ibms=',ibms
      IF(IBMS.NE.0) THEN
        NBM=0
        DO I=1,NF
          IF(LBM(I)) THEN
            IBM(I)=1
            NBM=NBM+1
          ELSE
            IBM(I)=0
          ENDIF
        ENDDO
        IF(NBM.EQ.NF) IPDS(7)=0
      ENDIF
c	write(*,*) 'gribit2 nbm=',nbm
      IF(NBM.EQ.0) THEN
        DO I=1,NF
          FR(I)=0.
        ENDDO
        NBIT=0
      ELSE
	fmin = f(1)
	fmax = fmin
	do ii = 1, nf
	    fmin = amin1(fmin,f(ii))
	    fmax = amax1(fmax,f(ii))
	enddo
c	write(*,*) 'fmin/fmax=',fmin,fmax

        CALL GTBITS(IPDS(7),IDS,NF,IBM,F,FR,FMIN,FMAX,NBIT)
c	write(*,*) 'fmin/max=',fmin,fmax,' nbit=',nbit
        IF(MXBIT.GT.0) NBIT=MIN(NBIT,MXBIT)
      ENDIF
c	write(*,*)'gribit4'
      CALL W3FI72(0,FR,0,NBIT,0,IPDS,PDS,
     &            1,255,IGDS,0,0,IBM,NF,IBDS,
     &            NFO,GRIB,LGRIB,IERR)
c	write(*,*)'gribit5'
      RETURN
      END

CFPP$ NOCONCUR R
      SUBROUTINE GTBITS(IBM,IDS,LEN,MG,G,GROUND,GMIN,GMAX,NBIT)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    GTBITS      COMPUTE NUMBER OF BITS AND ROUND FIELD.
C   PRGMMR: IREDELL          ORG: W/NMC23    DATE: 92-10-31
C
C ABSTRACT: THE NUMBER OF BITS REQUIRED TO PACK A GIVEN FIELD
C   AT A PARTICULAR DECIMAL SCALING IS COMPUTED USING THE FIELD RANGE.
C   THE FIELD IS ROUNDED OFF TO THE DECIMAL SCALING FOR PACKING.
C   THE MINIMUM AND MAXIMUM ROUNDED FIELD VALUES ARE ALSO RETURNED.
C   GRIB BITMAP MASKING FOR VALID DATA IS OPTIONALLY USED.
C
C PROGRAM HISTORY LOG:
C   92-10-31  IREDELL
C
C USAGE:    CALL GTBITS(IBM,IDS,LEN,MG,G,GMIN,GMAX,NBIT)
C   INPUT ARGUMENT LIST:
C     IBM      - INTEGER BITMAP FLAG (=0 FOR NO BITMAP)
C     IDS      - INTEGER DECIMAL SCALING
C                (E.G. IDS=3 TO ROUND FIELD TO NEAREST MILLI-VALUE)
C     LEN      - INTEGER LENGTH OF THE FIELD AND BITMAP
C     MG       - INTEGER (LEN) BITMAP IF IBM=1 (0 TO SKIP, 1 TO KEEP)
C     G        - REAL (LEN) FIELD
C
C   OUTPUT ARGUMENT LIST:
C     GROUND   - REAL (LEN) FIELD ROUNDED TO DECIMAL SCALING
C     GMIN     - REAL MINIMUM VALID ROUNDED FIELD VALUE
C     GMAX     - REAL MAXIMUM VALID ROUNDED FIELD VALUE
C     NBIT     - INTEGER NUMBER OF BITS TO PACK
C
C SUBPROGRAMS CALLED:
C   ISRCHNE  - FIND FIRST VALUE IN AN ARRAY NOT EQUAL TO TARGET VALUE
C
C ATTRIBUTES:
C   LANGUAGE: CRAY FORTRAN
C
C$$$
      DIMENSION MG(LEN),G(LEN),GROUND(LEN)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  ROUND FIELD AND DETERMINE EXTREMES WHERE BITMAP IS ON
      DS=10.**IDS
      IF(IBM.EQ.0) THEN
        GROUND(1)=NINT(G(1)*DS)/DS
        GMAX=GROUND(1)
        GMIN=GROUND(1)
        DO I=2,LEN
          GROUND(I)=NINT(G(I)*DS)/DS
          GMAX=MAX(GMAX,GROUND(I))
          GMIN=MIN(GMIN,GROUND(I))
        ENDDO
      ELSE
        I1=ISRCHNE(LEN,MG,1,0)
        IF(I1.GT.0.AND.I1.LE.LEN) THEN
          GROUND(I1)=NINT(G(I1)*DS)/DS
          GMAX=GROUND(I1)
          GMIN=GROUND(I1)
          DO I=I1+1,LEN
            IF(MG(I).NE.0) THEN
              GROUND(I)=NINT(G(I)*DS)/DS
              GMAX=MAX(GMAX,GROUND(I))
              GMIN=MIN(GMIN,GROUND(I))
            ENDIF
          ENDDO
        ELSE
          GMAX=0.
          GMIN=0.
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  COMPUTE NUMBER OF BITS
      NBIT=LOG((GMAX-GMIN)*DS+0.9)/LOG(2.)+1.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
      SUBROUTINE IDSDEF(IPTV,IDS)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM: IDSDEF         SETS DEFAULT DECIMAL SCALINGS
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 92-10-31
C
C ABSTRACT: SETS DECIMAL SCALINGS DEFAULTS FOR VARIOUS PARAMETERS.
C   A DECIMAL SCALING OF -3 MEANS DATA IS PACKED IN KILO-SI UNITS.
C
C PROGRAM HISTORY LOG:
C   92-10-31  IREDELL
C
C USAGE:    CALL IDSDEF(IPTV,IDS)
C   INPUT ARGUMENTS:
C     IPTV         PARAMTER TABLE VERSION (ONLY 1 OR 2 IS RECOGNIZED)
C   OUTPUT ARGUMENTS:
C     IDS          INTEGER (255) DECIMAL SCALINGS
C                  (UNKNOWN DECIMAL SCALINGS WILL NOT BE SET)
C
C ATTRIBUTES:
C   LANGUAGE: CRAY FORTRAN
C
C$$$
      DIMENSION IDS(255)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IPTV.EQ.1.OR.IPTV.EQ.2.or.IPTV.eq.132) THEN
        IDS(001)=-1     ! PRESSURE (PA)
        IDS(002)=-1     ! SEA-LEVEL PRESSURE (PA)
        IDS(003)=5      ! PRESSURE TENDENCY (PA/S)
                        !
                        !
        IDS(006)=-1     ! GEOPOTENTIAL (M2/S2)
        IDS(007)=0      ! GEOPOTENTIAL HEIGHT (M)
        IDS(008)=0      ! GEOMETRIC HEIGHT (M)
        IDS(009)=0      ! STANDARD DEVIATION OF HEIGHT (M)
                        !
        IDS(011)=1      ! TEMPERATURE (K)
        IDS(012)=1      ! VIRTUAL TEMPERATURE (K)
        IDS(013)=1      ! POTENTIAL TEMPERATURE (K)
        IDS(014)=1      ! PSEUDO-ADIABATIC POTENTIAL TEMPERATURE (K)
        IDS(015)=1      ! MAXIMUM TEMPERATURE (K)
        IDS(016)=1      ! MINIMUM TEMPERATURE (K)
        IDS(017)=1      ! DEWPOINT TEMPERATURE (K)
        IDS(018)=1      ! DEWPOINT DEPRESSION (K)
        IDS(019)=4      ! TEMPERATURE LAPSE RATE (K/M)
        IDS(020)=0      ! VISIBILITY (M)
                        ! RADAR SPECTRA 1 ()
                        ! RADAR SPECTRA 2 ()
                        ! RADAR SPECTRA 3 ()
                        !
        IDS(025)=1      ! TEMPERATURE ANOMALY (K)
        IDS(026)=-1     ! PRESSURE ANOMALY (PA)
        IDS(027)=0      ! GEOPOTENTIAL HEIGHT ANOMALY (M)
                        ! WAVE SPECTRA 1 ()
                        ! WAVE SPECTRA 2 ()
                        ! WAVE SPECTRA 3 ()
        IDS(031)=0      ! WIND DIRECTION (DEGREES)
        IDS(032)=1      ! WIND SPEED (M/S)
        IDS(033)=1      ! ZONAL WIND (M/S)
        IDS(034)=1      ! MERIDIONAL WIND (M/S)
        IDS(035)=-4     ! STREAMFUNCTION (M2/S) WNE
        IDS(036)=-4     ! VELOCITY POTENTIAL (M2/S) WNE
        IDS(037)=-1     ! MONTGOMERY STREAM FUNCTION (M2/S2)
        IDS(038)=8      ! SIGMA VERTICAL VELOCITY (1/S)
        IDS(039)=3      ! PRESSURE VERTICAL VELOCITY (PA/S)
        IDS(040)=4      ! GEOMETRIC VERTICAL VELOCITY (M/S)
        IDS(041)=6      ! ABSOLUTE VORTICITY (1/S)
        IDS(042)=6      ! ABSOLUTE DIVERGENCE (1/S)
        IDS(043)=6      ! RELATIVE VORTICITY (1/S)
        IDS(044)=6      ! RELATIVE DIVERGENCE (1/S)
        IDS(045)=4      ! VERTICAL U SHEAR (1/S)
        IDS(046)=4      ! VERTICAL V SHEAR (1/S)
        IDS(047)=0      ! DIRECTION OF CURRENT (DEGREES)
                        ! SPEED OF CURRENT (M/S)
                        ! U OF CURRENT (M/S)
                        ! V OF CURRENT (M/S)
        IDS(051)=4      ! SPECIFIC HUMIDITY (KG/KG)
        IDS(052)=0      ! RELATIVE HUMIDITY (PERCENT)
        IDS(053)=4      ! HUMIDITY MIXING RATIO (KG/KG)
        IDS(054)=1      ! PRECIPITABLE WATER (KG/M2)
        IDS(055)=-1     ! VAPOR PRESSURE (PA)
        IDS(056)=-1     ! SATURATION DEFICIT (PA)
        IDS(057)=1      ! EVAPORATION (KG/M2)
        IDS(058)=1      ! CLOUD ICE (KG/M2)
        IDS(059)=6      ! PRECIPITATION RATE (KG/M2/S)
        IDS(060)=0      ! THUNDERSTORM PROBABILITY (PERCENT)
        IDS(061)=1      ! TOTAL PRECIPITATION (KG/M2)
        IDS(062)=1      ! LARGE-SCALE PRECIPITATION (KG/M2)
        IDS(063)=1      ! CONVECTIVE PRECIPITATION (KG/M2)
        IDS(064)=6      ! WATER EQUIVALENT SNOWFALL RATE (KG/M2/S)
        IDS(065)=0      ! WATER EQUIVALENT OF SNOW DEPTH (KG/M2)
        IDS(066)=2      ! SNOW DEPTH (M)
                        ! MIXED-LAYER DEPTH (M)
                        ! TRANSIENT THERMOCLINE DEPTH (M)
                        ! MAIN THERMOCLINE DEPTH (M)
                        ! MAIN THERMOCLINE ANOMALY (M)
        IDS(071)=0      ! TOTAL CLOUD COVER (PERCENT)
        IDS(072)=0      ! CONVECTIVE CLOUD COVER (PERCENT)
        IDS(073)=0      ! LOW CLOUD COVER (PERCENT)
        IDS(074)=0      ! MIDDLE CLOUD COVER (PERCENT)
        IDS(075)=0      ! HIGH CLOUD COVER (PERCENT)
        IDS(076)=1      ! CLOUD WATER (KG/M2)
                        !
        IDS(078)=1      ! CONVECTIVE SNOW (KG/M2)
        IDS(079)=1      ! LARGE SCALE SNOW (KG/M2)
        IDS(080)=1      ! WATER TEMPERATURE (K)
        IDS(081)=0      ! SEA-LAND MASK ()
                        ! DEVIATION OF SEA LEVEL FROM MEAN (M)
        IDS(083)=5      ! ROUGHNESS (M)
        IDS(084)=0      ! ALBEDO (PERCENT)
        IDS(085)=1      ! SOIL TEMPERATURE (K)
        IDS(086)=0      ! SOIL WETNESS (KG/M2)
        IDS(087)=0      ! VEGETATION (PERCENT)
                        ! SALINITY (KG/KG)
        IDS(089)=4      ! DENSITY (KG/M3)
        IDS(090)=1      ! RUNOFF (KG/M2)
        IDS(091)=0      ! ICE CONCENTRATION ()
                        ! ICE THICKNESS (M)
        IDS(093)=0      ! DIRECTION OF ICE DRIFT (DEGREES)
                        ! SPEED OF ICE DRIFT (M/S)
                        ! U OF ICE DRIFT (M/S)
                        ! V OF ICE DRIFT (M/S)
                        ! ICE GROWTH (M)
                        ! ICE DIVERGENCE (1/S)
        IDS(099)=1      ! SNOW MELT (KG/M2)
                        ! SIG HEIGHT OF WAVES AND SWELL (M)
        IDS(101)=0      ! DIRECTION OF WIND WAVES (DEGREES)
                        ! SIG HEIGHT OF WIND WAVES (M)
                        ! MEAN PERIOD OF WIND WAVES (S)
        IDS(104)=0      ! DIRECTION OF SWELL WAVES (DEGREES)
                        ! SIG HEIGHT OF SWELL WAVES (M)
                        ! MEAN PERIOD OF SWELL WAVES (S)
        IDS(107)=0      ! PRIMARY WAVE DIRECTION (DEGREES)
                        ! PRIMARY WAVE MEAN PERIOD (S)
        IDS(109)=0      ! SECONDARY WAVE DIRECTION (DEGREES)
                        ! SECONDARY WAVE MEAN PERIOD (S)
        IDS(111)=0      ! NET SOLAR RADIATIVE FLUX AT SURFACE (W/M2)
        IDS(112)=0      ! NET LONGWAVE RADIATIVE FLUX AT SURFACE (W/M2)
        IDS(113)=0      ! NET SOLAR RADIATIVE FLUX AT TOP (W/M2)
        IDS(114)=0      ! NET LONGWAVE RADIATIVE FLUX AT TOP (W/M2)
        IDS(115)=0      ! NET LONGWAVE RADIATIVE FLUX (W/M2)
        IDS(116)=0      ! NET SOLAR RADIATIVE FLUX (W/M2)
        IDS(117)=0      ! TOTAL RADIATIVE FLUX (W/M2)
                        !
                        !
                        !
        IDS(121)=0      ! LATENT HEAT FLUX (W/M2)
        IDS(122)=0      ! SENSIBLE HEAT FLUX (W/M2)
        IDS(123)=0      ! BOUNDARY LAYER DISSIPATION (W/M2)
        IDS(124)=3      ! U WIND STRESS (N/M2)
        IDS(125)=3      ! V WIND STRESS (N/M2)
                        ! WIND MIXING ENERGY (J)
                        ! IMAGE DATA ()
        IDS(128)=-1     ! MEAN SEA-LEVEL PRESSURE (STDATM) (PA)
        IDS(129)=-1     ! MEAN SEA-LEVEL PRESSURE (MAPS) (PA)
        IDS(130)=-1     ! MEAN SEA-LEVEL PRESSURE (ETA) (PA)
        IDS(131)=1      ! SURFACE LIFTED INDEX (K)
        IDS(132)=1      ! BEST LIFTED INDEX (K)
        IDS(133)=1      ! K INDEX (K)
        IDS(134)=1      ! SWEAT INDEX (K)
        IDS(135)=10     ! HORIZONTAL MOISTURE DIVERGENCE (KG/KG/S)
        IDS(136)=4      ! SPEED SHEAR (1/S)
        IDS(137)=5      ! 3-HR PRESSURE TENDENCY (PA/S)
        IDS(138)=6      ! BRUNT-VAISALA FREQUENCY SQUARED (1/S2)
        IDS(139)=11     ! POTENTIAL VORTICITY (MASS-WEIGHTED) (1/S/M)
        IDS(140)=0      ! RAIN MASK ()
        IDS(141)=0      ! FREEZING RAIN MASK ()
        IDS(142)=0      ! ICE PELLETS MASK ()
        IDS(143)=0      ! SNOW MASK ()
        IDS(144)=4      ! SOILW
                        !
                        !
                        !
                        !
                        !
                        !
                        ! COVARIANCE BETWEEN V AND U (M2/S2)
                        ! COVARIANCE BETWEEN U AND T (K*M/S)
                        ! COVARIANCE BETWEEN V AND T (K*M/S)
                        !
                        !
        IDS(155)=0      ! GROUND HEAT FLUX (W/M2)
        IDS(156)=0      ! CONVECTIVE INHIBITION (W/M2)
                        ! CONVECTIVE APE (J/KG)
                        ! TURBULENT KE (J/KG)
                        ! CONDENSATION PRESSURE OF LIFTED PARCEL (PA)
        IDS(160)=0      ! CLEAR SKY UPWARD SOLAR FLUX (W/M2)
        IDS(161)=0      ! CLEAR SKY DOWNWARD SOLAR FLUX (W/M2)
        IDS(162)=0      ! CLEAR SKY UPWARD LONGWAVE FLUX (W/M2)
        IDS(163)=0      ! CLEAR SKY DOWNWARD LONGWAVE FLUX (W/M2)
        IDS(164)=0      ! CLOUD FORCING NET SOLAR FLUX (W/M2)
        IDS(165)=0      ! CLOUD FORCING NET LONGWAVE FLUX (W/M2)
        IDS(166)=0      ! VISIBLE BEAM DOWNWARD SOLAR FLUX (W/M2)
        IDS(167)=0      ! VISIBLE DIFFUSE DOWNWARD SOLAR FLUX (W/M2)
        IDS(168)=0      ! NEAR IR BEAM DOWNWARD SOLAR FLUX (W/M2)
        IDS(169)=0      ! NEAR IR DIFFUSE DOWNWARD SOLAR FLUX (W/M2)
			!
        IDS(170)=3      ! old (but current as of 1/94) u-stress (n/m**2)
        IDS(171)=3      ! old (but current as of 1/94) v-stress (n/m**2)
                        !
                        !
        IDS(172)=3      ! MOMENTUM FLUX (N/M2)
        IDS(173)=0      ! MASS POINT MODEL SURFACE ()
        IDS(174)=0      ! VELOCITY POINT MODEL SURFACE ()
        IDS(175)=0      ! SIGMA LAYER NUMBER ()
        IDS(176)=2      ! LATITUDE (DEGREES)
        IDS(177)=2      ! EAST LONGITUDE (DEGREES)

        IDS(178)=2      ! UMAS wne
        IDS(179)=2      ! VMAS wne
                        !
                        !
                        !
        IDS(181)=9      ! X-GRADIENT LOG PRESSURE (1/M)
        IDS(182)=9      ! Y-GRADIENT LOG PRESSURE (1/M)
        IDS(183)=5      ! X-GRADIENT HEIGHT (M/M)
        IDS(184)=5      ! Y-GRADIENT HEIGHT (M/M)
                        !
                        !
                        !
                        !
                        !
                        !
                        !
                        !
                        !
                        !
                        !
                        !
                        !
                        !
                        !
                        !
        IDS(201)=0      ! ICE-FREE WATER SURCACE (PERCENT)
                        !
                        !
        IDS(204)=0      ! DOWNWARD SOLAR RADIATIVE FLUX (W/M2)
        IDS(205)=0      ! DOWNWARD LONGWAVE RADIATIVE FLUX (W/M2)
                        !
        IDS(207)=0      ! MOISTURE AVAILABILITY (PERCENT)
                        ! EXCHANGE COEFFICIENT (KG/M2/S)
        IDS(209)=0      ! NUMBER OF MIXED LAYER NEXT TO SFC ()
                        !
        IDS(211)=0      ! UPWARD SOLAR RADIATIVE FLUX (W/M2)
        IDS(212)=0      ! UPWARD LONGWAVE RADIATIVE FLUX (W/M2)
        IDS(213)=0      ! NON-CONVECTIVE CLOUD COVER (PERCENT)
        IDS(214)=6      ! CONVECTIVE PRECIPITATION RATE (KG/M2/S)
        IDS(215)=7      ! TOTAL DIABATIC HEATING RATE (K/S)
        IDS(216)=7      ! TOTAL RADIATIVE HEATING RATE (K/S)
        IDS(217)=7      ! TOTAL DIABATIC NONRADIATIVE HEATING RATE (K/S)
        IDS(218)=2      ! PRECIPITATION INDEX (FRACTION)
        IDS(219)=1      ! STD DEV OF IR T OVER 1X1 DEG AREA (K)
        IDS(220)=4      ! NATURAL LOG OF SURFACE PRESSURE OVER 1 KPA ()
                        !
        IDS(222)=0      ! 5-WAVE GEOPOTENTIAL HEIGHT (M)
        IDS(223)=1      ! PLANT CANOPY SURFACE WATER (KG/M2)
                        !
                        !
                        ! BLACKADARS MIXING LENGTH (M)
                        ! ASYMPTOTIC MIXING LENGTH (M)
        IDS(228)=1      ! POTENTIAL EVAPORATION (KG/M2)
        IDS(229)=0      ! SNOW PHASE-CHANGE HEAT FLUX (W/M2)
                        !
        IDS(231)=3      ! CONVECTIVE CLOUD MASS FLUX (PA/S)
        IDS(232)=0      ! DOWNWARD TOTAL RADIATION FLUX (W/M2)
        IDS(233)=0      ! UPWARD TOTAL RADIATION FLUX (W/M2)
        IDS(224)=1      ! BASEFLOW-GROUNDWATER RUNOFF (KG/M2)
        IDS(225)=1      ! STORM SURFACE RUNOFF (KG/M2)
        IDS(234)=1      ! BASEFLOW-GROUNDWATER RUNOFF (KG/M2)
                        !
                        !
        IDS(238)=0      ! SNOW COVER (PERCENT)
        IDS(239)=1      ! SNOW TEMPERATURE (K)
                        !
        IDS(241)=7      ! LARGE SCALE CONDENSATION HEATING RATE (K/S)
        IDS(242)=7      ! DEEP CONVECTIVE HEATING RATE (K/S)
        IDS(243)=10     ! DEEP CONVECTIVE MOISTENING RATE (KG/KG/S)
        IDS(244)=7      ! SHALLOW CONVECTIVE HEATING RATE (K/S)
        IDS(245)=10     ! SHALLOW CONVECTIVE MOISTENING RATE (KG/KG/S)
        IDS(246)=7      ! VERTICAL DIFFUSION HEATING RATE (KG/KG/S)
        IDS(247)=7      ! VERTICAL DIFFUSION ZONAL ACCELERATION (M/S/S)
        IDS(248)=7      ! VERTICAL DIFFUSION MERID ACCELERATION (M/S/S)
        IDS(249)=10     ! VERTICAL DIFFUSION MOISTENING RATE (KG/KG/S)
        IDS(250)=7      ! SOLAR RADIATIVE HEATING RATE (K/S)
        IDS(251)=7      ! LONGWAVE RADIATIVE HEATING RATE (K/S)
                        ! DRAG COEFFICIENT ()
                        ! FRICTION VELOCITY (M/S)
                        ! RICHARDSON NUMBER ()
                        !
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
          FUNCTION ISRCHNE(N,IX,INCX,ITARGET)
          INTEGER IX(*),ITARGET
          J=1
          ISRCHNE=0
          IF(N.LE.0) RETURN
          IF(INCX.LT.0) J=1-(N-1)*INCX
          DO I=1,N
            IF(IX(J).NE.ITARGET) THEN
              ISRCHNE=I
              RETURN
            ENDIF
            J=J+INCX
          ENDDO
          RETURN
          END

