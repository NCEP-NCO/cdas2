      subroutine w3fs03(idate,ihour,iyear,month,iday,nn)
c$$$   subprogram  documentation  block
c
c subprogram: w3fs03         nmc date word packer and unpacker
c   author: jones,r.e.       org: w342       date: 87-03-24
c
c abstract: obtains the components of the nmc date word (see nmc
c   o.n. 84 and 85) or given its components, forms an nmc type
c   date word. w3fs03 is the same as w3fs11 except for the order of
c   the parameters.
c
c program history log:
c 87-03-24  r.e.jones   convert to cyber 205 fortran 200
c 89-10-13  r.e.jones   convert to cray cft77 fortran
c
c usage:  call w3fs03 (idate,ihour,iyear,month,iday,nn)
c
c   input variables:
c     names  interface description of variables and types
c     ------ --------- -----------------------------------------------
c     idate  arg list  left 4 bytes of integer word
c                      word. (7th word of data field or 3rd word of
c                      the id table of a binary file if in half
c                      precision array or the 4th word or 2nd word if
c                      in integer array).
c     ihour  arg list  hour
c     iyear  arg list  year (2 digits)
c     month  arg list  month
c     iday   arg list  day
c     nn     arg list  code:
c                       = 0  pack ihour, iyear, month, iday, into idate
c                      <> 0  unpack idate into ihour, iyear, month, iday
c
c   output variables:
c     names  interface description of variables and types
c     ------ --------- -----------------------------------------------
c     idate  arg list  left 4 bytes of integer word
c                      word. (7th word of data field or 3rd word of
c                      the id table of a binary file if in half
c                      precision array or the 4th word or 2nd word if
c                      in integer array).
c     ihour  arg list  hour
c     iyear  arg list  year (2 digits)
c     month  arg list  month
c     iday   arg list  day
c
c   subprgrams called:
c     names   library
c     ------------------------------------------------------- --------
c     char mova2i                                              system
c
c remarks:    when nn.ne.0, the information in idate must be
c     formatted as diagrammed in appendix c of nmc o.n. 84.
c
c attributes:
c   language: cray cft77 fortran
c   machine:  cyay y-mp8/832
c
c$$$
c
      character*1 idate(4)
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
