!=================================
!==  vel-smooth
!==             Rt
!==      nx_cdp:输出道对比
!==      nzs：纵向开始平滑深度
!==      nsmooth：平滑次数
!=================================
program smooth
   parameter(nx=560,nz=360,nx_cdp=80,nzs=40,nsmooth=200)
   dimension vp1(1:nx,1:nz),vp2(1:nx,1:nz)
   real:: vp1,vp2,vmin,vmax
   integer :: i,j,k

   character*256 fn1,fn2,fn3,fn4


   fn1='fault_vel_560_360.dat'
   fn2='fault_vel_560_360_initial.dat'
   fn3='CDP_line_compare.txt'

   
   open(10,file=fn1,access='direct',recl=4*nz)
        do i=1,nx
	   read(10,rec=i)(vp2(i,j),j=1,nz)
	enddo
    close(10)
   vmin=vp2(1,1)
   vmax=vp2(1,1)
   do i=1,nx
      do j=1,nz
        if (vmin>vp2(i,j))then
           vmin=vp2(i,j)
        elseif (vmax<vp2(i,j))then
           vmax=vp2(i,j)
        endif
      enddo
   enddo
   print *,vmin,vmax 
    do i=1,nx
        do j=1,nz
            vp1(i,j)=vp2(i,j)
        enddo
    enddo
  

  do k=1,nsmooth

    vp1(1,nzs)=(2*vp1(2,nzs)+2*vp1(1,nzs+1)+vp1(2,nzs+1)+4*vp1(1,nzs))/9
    vp1(1,nz)=(2*vp1(1,nz-1)+2*vp1(2,nz)+vp1(2,nz-1)+4*vp1(1,nz))/9
    vp1(nx,nzs)=(2*vp1(nx-1,nzs)+vp1(nx-1,nzs+1)+2*vp1(nx,nzs+1)+4*vp1(nx,nzs))/9
    vp1(nx,nz)=(vp1(nx-1,nz-1)+2*vp1(nx,nz-1)+2*vp1(nx-1,nz)+4*vp1(nx,nz))/9

    do i=2,nx-1
        vp1(i,nzs)=(2*vp1(i-1,nzs)+vp1(i-1,nzs+1)+2*vp1(i,nzs+1)+vp1(i+1,nzs+1)+2*vp1(i+1,nzs)+4*vp1(i,nzs))/12
        vp1(i,nz)=(2*vp1(i-1,nz)+vp1(i-1,nz-1)+2*vp1(i,nz-1)+vp1(i+1,nz-1)+2*vp1(i+1,nz)+4*vp1(i,nz))/12
    enddo

    do j=nzs+1,nz-1
        vp1(1,j)=(2*vp1(1,j-1)+vp1(2,j-1)+2*vp1(2,j)+vp1(2,j+1)+2*vp1(1,j+1)+4*vp1(1,j))/12
        vp1(nx,j)=(2*vp1(nx,j-1)+vp1(nx-1,j-1)+2*vp1(nx-1,j)+vp1(nx-1,j+1)+2*vp1(nx,j+1)+4*vp1(nx,j))/12
    enddo

    do i=2,nx-1
       do j=nzs+1,nz-1
vp1(i,j)=(vp1(i-1,j-1)+vp1(i,j-1)*2+vp1(i+1,j-1)+vp1(i-1,j)*2+vp1(i,j)*4+vp1(i+1,j)*2+vp1(i-1,j+1)+vp1(i,j+1)*2+vp1(i+1,j+1))/16
       enddo
    enddo

enddo
    open(11,file=fn2,access='direct',recl=4*nz)
        do i=1,nx
	   write(11,rec=i)(vp1(i,j),j=1,nz)
	   enddo
    close(11)

    open(14,file=fn3,status='unknown')
    do j=1,nz
    write(14,*)(vp2(nx_cdp,j)),(vp1(nx_cdp,j))
   ! write(14,*)(vp1(nx_cdp,j))
    enddo
    close(14)


end

