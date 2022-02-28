      function sdot(n,x,incx,y,incy)
      dimension x(n),y(n)
      sdot=0.
      ix=1
      iy=1
      do i=1,n
        sdot=sdot+x(ix)*y(iy)
        ix=ix+incx
        iy=iy+incy
      enddo
      return
      end
