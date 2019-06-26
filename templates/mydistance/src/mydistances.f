      SUBROUTINE aniso_euclidean_s(D,x,y,nx,ny,ndx,ndy,cmin,cmax,sr,
     *nsr,symm)

cf2py intent(inplace) D
cf2py integer intent(optional) :: cmin=0
cf2py integer intent(optional) :: cmax=-1
cf2py logical intent(optional) :: symm=0
cf2py intent(hide) nx, ny,ndx,ndy,nsr
cf2py threadsafe

      DOUBLE PRECISION D(nx,ny), x(nx,ndx), y(ny,ndy), sr(nsr)
      INTEGER nx,ny,i,j,k,cmin,cmax,ndx,ndy,nsr
      LOGICAL symm
      DOUBLE PRECISION dist, dev

      if (cmax.EQ.-1) then
          cmax = ny
      end if
      ! print *,nsr,sr(1)

      if(symm) then

        do j=cmin+1,cmax
          D(j,j) = 0.0D0
          do i=1,j-1
            dist = 0.0D0
            do k=1,ndx
              if (k.GT.1) then
                 dev=(x(i,k) - y(j,k))*sr(k-1)
              else
                 dev=(x(i,k) - y(j,k))
              end if
              dist = dist + dev*dev
            enddo
            D(i,j) = dsqrt(dist)
          enddo
        enddo
      else
        do j=cmin+1,cmax
          do i=1,nx
            dist = 0.0D0
            do k=1,ndx
              dev=(x(i,k) - y(j,k))*sr(k)
              dist = dist + dev*dev
            enddo
            D(i,j) = dsqrt(dist)
          enddo    
        enddo  
      endif
      RETURN
      END
