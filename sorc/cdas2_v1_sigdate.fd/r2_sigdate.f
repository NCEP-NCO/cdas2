	program sigdate
c
c	reads the date code from a sigma file
c
	integer hour, month, day, year
c
	read(11)
c
	read(11) fhour, hour, month, day, year
c
	ihour=nint(fhour)
c
	write(6,200) year, month, day, hour, ihour
200	format(i4,3(i2.2),2x,i10)
c
	stop
	end

