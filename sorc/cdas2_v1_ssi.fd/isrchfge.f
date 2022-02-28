      function isrchfge(n,a,is,val)
      dimension a(n)
      isrchfge=n
      do nn=is,n
      if( a(nn).ge.val ) then
        isrchfge=nn
        return
      endif
      enddo
      return
      end
