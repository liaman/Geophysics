# smooth for binary file
* smooth3D.f95 for 3D medium
* smooth.f95 for 2D medium
* smooth method as follows:
```f95
    do i=2,nx-1
        vp1(i,nzs)=(2*vp1(i-1,nzs)+vp1(i-1,nzs+1)+2*vp1(i,nzs+1)+vp1(i+1,nzs+1)+2*vp1(i+1,nzs)+4*vp1(i,nzs))/12
        vp1(i,nz)=(2*vp1(i-1,nz)+vp1(i-1,nz-1)+2*vp1(i,nz-1)+vp1(i+1,nz-1)+2*vp1(i+1,nz)+4*vp1(i,nz))/12
    enddo
```
