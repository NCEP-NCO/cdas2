      function isrchfle(n,a,is,val)
      dimension a(n)
      isrchfle=n
      do nn=is,n
      if( a(nn).le.val ) then
        isrchfle=nn
        return
      endif
      enddo
      return
      end
