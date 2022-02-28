       subroutine w3fs11(idate,iyear,month,iday,ihour,nn)
c$$$   subprogram  documentation  block
c
c subprogram: w3fs11         nmc date word unpacker and packer
c   author: jones,r.e.       org: w342       date: 87-03-24
c
c abstract: obtains the components of the nmc date word (nmc office
c   notes 84 and 85), or given its components, forms an nmc type date
c   word. w3fs11 is the same as w3fs03 except for the order of the
c   parameters.
c
c program history log:
c   87-03-24  r.e.jones   convert to cyber 205 fortran 200
c   89-10-13  r.e.jones   convert to cray cft77 fortran
c
c usage:  call w3fs11 (idate, iyear, month, iday, ihour, nn)
c
c   input variables:
c     names  interface description of variables and types
c     ------ --------- -----------------------------------------------
c     idate  arg list  left 4 bytes of integer 64 bit word, or can be
c                      character*1 idate(4) or character*4 idate.
c                      if office note 85 label used as 4 64 bit words,
c                      idate is in the left 32 bits of the 2nd word.
c                      if office note 84 12 ids. (6 64 bit words on
c                      cray, date word in left 32 bits of 4th cray id
c                      word, or 7th 32 bit id word on nas.
c     iyear  arg list  integer   year (2 digits)
c     month  arg list  integer   month
c     iday   arg list  integer   day
c     ihour  arg list  integer   hour
c     nn     arg list  integer   code:
c                     .eq. 0 pack iyear, month, iday, ihour into idate
c                     .ne. 0 unpack idate into iyear, month, iday, ihour
c
c   output variables:
c     names  interface description of variables and types
c     ------ --------- -----------------------------------------------
c     idate  arg list  left 4 bytes of integer 64 bit word, or can be
c                      character*1 idate(4) or character*4 idate.
c                      if office note 85 label used as 4 64 bit words,
c                      idate is in the left 32 bits of the 2nd word.
c                      if office note 84 12 ids. (6 64 bit words on
c                      cray, date word in left 32 bits of 4th cray id
c                      word, or 7th 32 bit id word on nas.
c     iyear  arg list  integer   year (2 digits)
c     month  arg list  integer   month
c     iday   arg list  integer   day
c     ihour  arg list  integer   hour
c
c   subrograms called:
c     names                                                   library
c     ------------------------------------------------------- --------
c     char   mova2i                                            system
c
c attributes:
c   language: cray cft77 fortran
c   machine:  cray y-mp8/832
c
c$$$
c
      character idate(4)
c
      if (nn.ne.0) then
c
        iyear = mova2i(idate(1))
        month = mova2i(idate(2))
        iday  = mova2i(idate(3))
        ihour = mova2i(idate(4))
c
      else
c
        idate(1) = char(iyear)
        idate(2) = char(month)
        idate(3) = char(iday)
        idate(4) = char(ihour)
      endif
c
      return
      end
