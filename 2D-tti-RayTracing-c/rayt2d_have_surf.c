/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* RAYT2D $Revision: 1.18 $; Date: 94/10/11 $	*/
#include "/home/rongtao/gpfs03/hc/cjbsegy.h"
#include "/home/rongtao/gpfs03/hc/fft.c"
#include "/home/rongtao/gpfs03/hc/alloc.c"
#include "/home/rongtao/gpfs03/hc/complex.c"

/*//#include "ealloc.c"
//#include "par.h"
//#include "getpars.c"
//#include "intcub.c"
//#include "cubicspline.c"
//#include "xindex.c"
//#include "cwp.h"*/

/* define initial value for traveltime as an arbitrary large number 
       定义初始值为走时任意大值            */
#define INITIAL_T 999999

/* Structures used internally */
/* one  ray */
typedef struct RayStruct {
	float t;		/* time 时间 */
	float px,pz;		/* slowness componets 慢度分量 */
	float x,z;		/* x,z coordinates x，z坐标*/
	float q2,p2;	        /* Cerveny's dynamic ray tracing solution 
                                       切尔韦尼动态射线追踪方法    */
	float v,dvdx,dvdz;	/* velocity and its derivatives 
                                       速度及其导数               */
	float dtv;	        /* traveltime variation 
                                       走时的变化                 */
} Ray;

/* ray tracing geometry parameters 射线追踪的几何参数 */
typedef struct RaytrStruct {
	float xs;		/* source lateral location 震源横向位置 */	 
	float zs;		/* source depth location   震源深度位置 */	 
	int nt;			/* number of time samples in ray tracing 
                                                           时间采样  */
	float dt;		/* time sample interval in ray tracing 
                                                           时间采样间隔 */	 
	int na;			/* number of rays          射线条数  */
	float fa;		/* first take-off angle of rays */
	float da;		/* increment of take-off angles */
	float amin;		/* minimum emergence angle  */	 
	float amax;		/* maximum emergence angle  */	 
	float fac;		/* factor to determine extrapolation radius
                                      确定外推半径的因子  */ 
} Raytr;

/* 2D section geometry parameters 二维剖面几何参数 */
typedef struct Geo2DStruct {
	int nx;			/* number of lateral samples  */
	float fx;		/* first lateral sample  */	 
	float dx;		/* lateral interval   */	 
	float odx;		/* 1/dx   */	 
	int nz;			/* number of depth samples */
	float fz;		/* first depth sample  */	 
	float dz;		/* depth interval  */	
	float odz;		/* 1/dz   */	 
}Geo2d;

/* Geometry of the recording surface 记录表面的几何形状 */
typedef struct SurfaceSegmentStruct {
        float x;        /* x coordinate of segment midpoint 分割中点的x坐标*/
        float z;        /* z coordinate of segment midpoint 分割中点的z坐标*/
        float s;        /* x component of unit-normal-vector x分量的单位法向量*/
        float c;        /* z component of unit-normal-vector z分量的单位法向量*/
} SurfaceSegment;

typedef struct SurfaceStruct {
        int ns;               /* number of surface segments 表面数 */
        float ds;             /* segment length 段长度 */
        SurfaceSegment *ss;   /* array[ns] of surface segments */
} Surface;

/* Prototypes of functions used internally */
void makeRay (Geo2d geo2dv,float *v,float *vxx,float *vxz,float *vzz,
	Raytr raytr,float a0,int *nrs,Ray *ray,int npv,float *pv);

void raytime2d(Raytr raytr,Geo2d geo2dv,float *vt,Geo2d geo2dt,float *t,
	int npv,float *vo,float *pv,float *tv,float *cs,FILE *fp);

void voint(Geo2d geo2dv,float *v,Geo2d geo2dt,float *ov2,int npv,float *vo);

void dv2(Geo2d geo2dv,float *v,float *vxx,float *vxz,float *vzz);

void eiknl(Geo2d geo2dt,float *t,float *ov2,float tmax);

void trans(int nx,int nz,int nxt,int nx0,float *v,float *vt);

void vel2Interp(Geo2d geo2dv,float *v,float *vxx,float *vxz,float *vzz,
     float x,float z,
     float *u,float *ux,float *uz,float *uxx,float *uxz,float *uzz);

void vintp(Geo2d geo2dv,float *v,float x,float z,float *u);

void decodeSurfaces(int *nrPtr,int **nxzPtr,float ***xPtr,float ***zPtr);

int decodeSurface(char *string,int *nxzPtr,float **xPtr,float **zPtr);

void makesurf(float dsmax,int nr,int *nu,float **xu,float **zu,
             Surface **r);
void zcoorTopog(float fxs,float dxs,int nxs,Surface *srf,float *sz,
                float *nangl);


void main(int argc, char **argv)
{
	int	na,nat,nt,nxs,nxo,nzo,nx,nz,nxt,nx0,mx,npv,nsrf,*nxzsrf;
	float	dt,xs,fxs,dxs,exs,fxo,fzo,dxo,dzo,exo,ezo,fa,ea,amin,eat,
		amax,da,fat,fac,tmax,aperx,temp,fx,fz,dx,dz,ex,ez,fxt,ext;
	float	*v,*vt,*t,*ov2,*pv=NULL,*pvt=NULL,*tv=NULL,*cs=NULL,*vo=NULL,
		**xsrf,**zsrf,*szi,*nangl;
	int	ixs,ixs0,nsize,isize,ek,ms;
        int     i,j,k;
    /*struct */
	Geo2d geo2dv, geo2dt;
	Raytr raytr;
        Surface *srf;


        char *pvfile,*tvfile,*csfile;
	FILE *vfp, *tfp, *jpfp, *pvfp=NULL, *tvfp=NULL, *csfp=NULL;


        char vfile[250]={"vel200200.dat"};
        char tfile[250]={"time.dat"};
        char rayfile[250]={"Raypath.txt"};

	vfp = fopen(vfile,"r");

        nx=200;    fx=0;     dx=5.0;
        nz=200;    fz=0;     dz=5.0;
 
        nt = 200001;
        dt = 0.000008;

        nxs=2;
        fxs = nx*dx/3;
	dxs = nx*dx/3;

	fa = -90.;  //first take-off angle of rays (degrees)	
	da = 1;    //increment of take-off angle 
	na = 180;    //number of rays
	amin = 0.;  //minimum emergence angle
	amax = 90.; //maximum emergence angle

	fac = 0.02;//factor to determine radius for extrapolation 确定外推半径

	ek = 1;  //flag of implementing eikonal in shadow zones在阴影区实现程函方程
	ms = 10; //print verbal information at every ms sources
	npv = 0; //flag of computing quantities for velocity analysis标记速度分析的计算量

	ixs0 = 0;
	isize = 0;
	nsize = nzo*nxo;


        tmax = (nt-1)*dt;
   
        nxo=nx;  //number of lateral samples in traveltime table	
        nzo=nz;  //number of depth samples in traveltime table	
      
        fxo=fx;  //first lateral sample in traveltime table	
        fzo=fz;  //first depth sample in traveltime table	
        dxo=dx;  //lateral interval in traveltime table	
        dzo=dz;  //depth interval in traveltime table	
         

	
	if(nx<3 || nz<3) err("nx and nz must be not less than 3!\n");
	if(nxo<3 || nzo<3) err("nxo and nzo must be not less than 3!\n");

	ex = fx+(nx-1)*dx;
	ez = fz+(nz-1)*dz;

      /* 2D section geometry parameters -velocity*/
	geo2dv.nx = nx;		geo2dv.fx = fx;	
	geo2dv.dx = dx;		geo2dv.odx = 1.0/dx;
	geo2dv.nz = nz;		geo2dv.fz = fz;	
	geo2dv.dz = dz;		geo2dv.odz = 1.0/dz;


	exo = fxo+(nxo-1)*dxo;
	ezo = fzo+(nzo-1)*dzo;

      /* 2D section geometry parameters -time*/
	geo2dt.nx = nxo;	geo2dt.fx = fxo;	
	geo2dt.dx = dxo;	geo2dt.odx = 1.0/dxo;
	geo2dt.nz = nzo;	geo2dt.fz = fzo;	
	geo2dt.dz = dzo;	geo2dt.odz = 1.0/dzo;




	exs = fxs+(nxs-1)*dxs;

        aperx = 0.5*(ex-fx);
	
	if(nxs>1) aperx += dxs;

	mx = 2.0*aperx/dx+2.99;

	if(mx>nx) mx = nx;


	

	if(amax>180 || amin<0)  
	err("amin and amx must be within 0 to 180 degrees!\n");

	fa *= PI/180;
	da *= PI/180;
	ea = fa+(na-1)*da;
	amin *= PI/180;
	amax *= PI/180;



      /* ray tracing geometry parameters */
	raytr.nt = nt;		raytr.dt = dt;
	raytr.da = da;		raytr.fac = fac;
	raytr.amin = amin;	raytr.amax = amax;


	

	
        tfp = fopen(tfile,"w");

	jpfp = fopen("Check_document.txt","w");
	
	
	if(npv){ 
		pvfile = "pvfile";
	    pvfp = fopen(pvfile,"r"); 
		tvfile = "tvfile";
	    tvfp = fopen(tvfile,"w");
		csfile = "csfile";
	    csfp = fopen(csfile,"w");
	}


	fprintf(jpfp,"\n");
	fprintf(jpfp," RAYT2D parameters\n");
	fprintf(jpfp," ================\n");
	fprintf(jpfp," vfile=%s \n",vfile);
	fprintf(jpfp," tfile=%s \n",tfile);
	fprintf(jpfp," one-way time: dt=%g nt=%d \n",dt,nt);
	fprintf(jpfp," nz=%d fz=%g dz=%g\n",nz,fz,dz);
	fprintf(jpfp," nx=%d fx=%g dx=%g\n",nx,fx,dx);
	fprintf(jpfp," nzo=%d fzo=%g dzo=%g\n",nzo,fzo,dzo);
	fprintf(jpfp," nxo=%d fxo=%g dxo=%g\n",nxo,fxo,dxo);
 	fprintf(jpfp," nxs=%d fxs=%g dxs=%g\n",nxs,fxs,dxs);
 	fprintf(jpfp," mx=%d aperx=%g \n",mx,aperx);
 	fprintf(jpfp," na=%d fa=%g da=%g\n",na,fa*180/PI,da*180/PI);
 	fprintf(jpfp," amin=%g amax=%g \n",amin*180/PI,amax*180/PI);
 	fprintf(jpfp," ek=%d ms=%d npv=%d \n",ek,ms,npv);
	if(npv)
 	  	fprintf(jpfp," pvfile=%s tvfile=%s csfile=%s\n",
			pvfile,tvfile,csfile);
	fflush(jpfp);

	/* ensure sources and output are in velocity grid*/
	if(fx>fxs || ex<exs || fz>0) {
            printf("This condition must NOT be satisfied: fx>fxs || ex<exs || fz>0\n");
            printf("fx=%g fxs=%g ex=%g exs=%g fz=%g\n",
                    fx,fxo,ex,exo,fz);
	    printf("source lies outside of specified x grid  \n");exit(0);
        }
	if(fx>fxo || ex<exo || fz>fzo || ez<ezo) {
            printf("NOT: fx>fxo||ex<exo||fz>fzo||ez<ezo\n");
            printf("fx=%g fxo=%g ex=%g exo=%g fz=%g fzo=%g,ez=%g ezo=%g\n",
                    fx,fxo,ex,exo,fz,fzo,ez,ezo);
	    printf("output lies outside of specified x grid  \n");exit(0);
        }

	/* allocate space */
	v = alloc1float(nx*nz);
	ov2 = alloc1float(nxo*nzo);
	t = alloc1float(nxo*nzo);
	if(npv) {
		pv = alloc1float(nx*nz);
		tv = alloc1float(nxo*nzo);
		vo = alloc1float(nxo*nzo);
		cs = alloc1float(nxo*nzo);
	}

	/* read velocity */
	if(fread(v,sizeof(float),nx*nz,vfp)!=nx*nz)
         {
	    printf("cannot read %d velocities from file %s",nx*nz,vfile);
            exit(0);
         }

	if(npv) 
	    if(fread(pv,sizeof(float),nx*nz,pvfp)!=nx*nz)
	    	printf("cannot read %d values from file %s",nx*nz,pvfile);

	fprintf(jpfp," finish velocity input\n");


	/*interpolate velocity and compute slowness squares 
           插入速度和计算慢度场*/
	voint(geo2dv,v,geo2dt,ov2,npv,vo);

	fprintf(jpfp," begin traveltime calculation\n");
	fflush(jpfp);
       
	/* Estimation of the recording surface 	*/
	//decodeSurfaces(&nsrf,&nxzsrf,&xsrf,&zsrf);//解码表面
        //makesurf(dxs*0.025,nsrf,nxzsrf,xsrf,zsrf,&srf);//创建表面
    
        szi = alloc1float(nxs);
	nangl = alloc1float(nxs);
        //zcoorTopog(fxs,dxs,nxs,srf,szi,nangl);
       
     FILE *fpraypath;
     fpraypath=fopen(rayfile,"w");


	/* loop over sources */
	for(ixs=ixs0,xs=fxs+ixs0*dxs;ixs<nxs;ixs++,xs+=dxs){
                printf("ixs=%4d/%d\n",ixs,nxs);

		/* redefine the velocity model due to aperture	*/
		temp = fxo;
		if(xs<temp) temp = xs;
		if(xs-aperx>temp) temp = xs-aperx;
		nx0 = (temp-fx)/dx;
		fxt = fx+nx0*dx;
		temp = exo;
		if(xs>temp) temp = xs;
		if(xs+aperx<temp) temp = xs+aperx;
		nxt = 1+(int)((temp-fx)/dx+0.99)-nx0;
		ext = fxt+(nxt-1)*dx;

		vt = alloc1float(nxt*nz);
		if(npv) pvt = alloc1float(nxt*nz);

	        /* reducing velocity volume due to aperture	*/
		trans(nx,nz,nxt,nx0,v,vt);
		if(npv) trans(nx,nz,nxt,nx0,pv,pvt);

		/* determine range of take-off angles	*/
		fat = fa + nangl[ixs];
		eat = ea + nangl[ixs];
		nat = na;
		if (fat<-(PI/2)) fat = -PI/2;
		if (eat>(PI/2)) eat = PI/2;
		if(xs==fxt && fat<0) {
			fat = 0.;
			nat = eat/da+1.5;
		} else if(xs==ext && eat>0)
			nat = 1.5-fa/da;


		/* update geometry information	*/
		raytr.xs = xs;		raytr.zs = szi[ixs];
		raytr.na = nat;		raytr.fa = fat;	
	
		geo2dv.nx = nxt;	geo2dv.fx = fxt;
		
	        /* compute traveltime by paraxial ray tracing	*/
		raytime2d(raytr,geo2dv,vt,geo2dt,t,npv,vo,pvt,tv,cs,fpraypath);

		/*make up in shadow zones by eikonal equation	
                     在阴影区实现程函方程     */
		if(ek) eiknl(geo2dt,t,ov2,tmax);
	
		/*write traveltime	*/
                for(i=0;i<nxo*nzo;i++)
                      if((t[i]==999999)||(t[i]>99999))
                             t[i]=0.0;
		fwrite(t,sizeof(float),nxo*nzo,tfp);

		/*write quantities for velocity analysis	*/
		if(npv) {
			fwrite(tv,sizeof(float),nxo*nzo,tvfp);
			fwrite(cs,sizeof(float),nxo*nzo,csfp);
		}

		if(ixs%ms==0) {
		  fprintf(jpfp," traveltime computed at source xs=%g,zs=%g\n"
                                  ,xs,szi[ixs]);
		  fflush(jpfp);
		}
		free1float(vt);
		if(npv) free1float(pvt);
	}
 	fprintf(jpfp," finish program rayt2d\n\n");
	
	fclose(tfp);
	fclose(vfp);
	fclose(jpfp);
	fclose(fpraypath);

	free1float(v);
	free1float(t);
	if(npv){
		free1float(pv);
		free1float(tv);
		free1float(cs);
	}


}

/**********************************************************************
        Subroutines adapted from Dave Hale's modeling library

        decodeReflectors
        decodeReflector
        makeref                 Autor: Dave Hale, CSM, 09/17/91
**********************************************************************/

void decodeSurfaces(int *nrPtr,int **nxzPtr,float ***xPtr,float ***zPtr)
/*************************************************************************
decodeSurfaces - parse surfaces parameter string 解析表面参数字符串
**************************************************************************
Output:
nrPtr           pointer to nr an int specifying number of surfaces = 2
nxzPtr          pointer to nxz specifying number of (x,z) pairs defining the
                surfaces
xPtr            pointer to array[x][nr] of x values for each surface
zPtr            array[z][nr] of z values for each surface
**************************************************************************/
{

        int nr,*nxz,ir;
        float **x,**z;
	char t[6144],*s;

        /* count surfaces */
        //nr = countparname("surf");
        if (nr==0) nr = 1;

        /* allocate space */
        nxz = alloc1(nr,sizeof(int));
        x = alloc1(nr,sizeof(float*));
        z = alloc1(nr,sizeof(float*));

        /* get surfaces */
        for (ir=0; ir<nr; ++ir) {
                //if (!getnparstring(ir+1,"surf",&s)) s = "0,0;99999,0";
                strcpy(t,s);
                if (!decodeSurface(t,&nxz[ir],&x[ir],&z[ir]))
                        err("Surface number %d specified "
                                "incorrectly!\n",ir+1);
        }

        /* set output parameters before returning */
        *nrPtr = nr;
        *nxzPtr = nxz;
        *xPtr = x;
        *zPtr = z;
}

/* parse one surface specification; return 1 if valid, 0 otherwise */

int decodeSurface (char *string,int *nxzPtr, float **xPtr, float **zPtr)
/**************************************************************************
decodeSurface - parse one surface specification
***************************************************************************
Input:
string          string representing surface

Output:
nxzPtr          pointer to number of x,z pairs
xPtr            array of x values for one surface
zPtr            array of z values for one surface
**************************************************************************/
{
        int nxz,ixz;
        float *x,*z;
        char *s,*t;

        s = string;
        s = strtok(s,",;\0");
        /* count x and z values, while splitting string into tokens */
        /*把字符串分解为字符*/
        for (t=s,nxz=0; t!=NULL; ++nxz)
                t = strtok(NULL,",;\0");

        /* total number of values must be even */
        if (nxz%2) return 0;

        /* number of (x,z) pairs */
        nxz /= 2;

        /* 2 or more (x,z) pairs are required */
        if (nxz<2) return 0;

        /* allocate space */
        x = alloc1(nxz,sizeof(float));
        z = alloc1(nxz,sizeof(float));

        /* convert (x,z) values */
        for (ixz=0; ixz<nxz; ++ixz) {
                x[ixz] = atof(s);//字符串转为浮点数
                s += strlen(s)+1;
                z[ixz] = atof(s);
                s += strlen(s)+1;
        }

        /* set output parameters before returning */
        *nxzPtr = nxz;
        *xPtr = x;
        *zPtr = z;
        return 1;
}

void makesurf(float dsmax,int nr,int *nu,float **xu,float **zu,
        Surface **r)
/*****************************************************************************
Make piecewise cubic surfaces
******************************************************************************
Input:
dsmax           maximum length of surface segment
nr              number of surfaces = 2
nu              array[nr] of numbers of (x,z) pairs; u = 0, 1, ..., nu[ir]
xu              array[nr][nu[ir]] of surface x coordinates x(u)
zu              array[nr][nu[ir]] of surface z coordinates z(u)

Output:
r               array[nr] of surfaces
******************************************************************************
Notes:
Space for the nu, xu, and zu arrays is freed by this function, since
they are only used to construct the surfaces.
*****************************************************************************/
{
   void yxtoxy (int nx, float dx, float fx, float y[], 
	int ny, float dy, float fy, float xylo, float xyhi, float x[]);
   void intcub (int ideriv, int nin, float xin[], float ydin[][4], 
	int nout, float xout[], float yout[]);
   void csplin (int n, float x[], float y[], float yd[][4]);

        int ir,iu,nuu,iuu,ns,is;
        float x,z,xlast,zlast,dx,dz,duu,uu,ds,fs,rsx,rsz,rsxd,rszd,
                *u,*s,(*xud)[4],(*zud)[4],*us;
        SurfaceSegment *ss;
        Surface *rr;

        /* allocate space for surfaces */
        *r = rr = alloc1(nr,sizeof(Surface));

        /* loop over surfaces */
        for (ir=0; ir<nr; ++ir) {

                /* compute cubic spline coefficients for uniformly sampled u */
                u = alloc1float(nu[ir]);
                for (iu=0; iu<nu[ir]; ++iu)
                        u[iu] = iu;
                xud = (float(*)[4])alloc1float(4*nu[ir]);
                csplin(nu[ir],u,xu[ir],xud);
                zud = (float(*)[4])alloc1float(4*nu[ir]);
                csplin(nu[ir],u,zu[ir],zud);

                /* finely sample x(u) and z(u) and compute length s(u) */
                nuu = 20*nu[ir];
                duu = (u[nu[ir]-1]-u[0])/(nuu-1);
                s = alloc1float(nuu);
                s[0] = 0.0;
                xlast = xu[ir][0];
                zlast = zu[ir][0];
                for (iuu=1,uu=duu; iuu<nuu; ++iuu,uu+=duu) {
                        intcub(0,nu[ir],u,xud,1,&uu,&x);
                        intcub(0,nu[ir],u,zud,1,&uu,&z);
                        dx = x-xlast;
                        dz = z-zlast;
                        s[iuu] = s[iuu-1]+sqrt(dx*dx+dz*dz);
                        xlast = x;
                        zlast = z;
                }

                /* compute u(s) from s(u) */
                ns = 1+s[nuu-1]/dsmax;
                ds = s[nuu-1]/ns;
                fs = 0.5*ds;
                us = alloc1float(ns);
                yxtoxy(nuu,duu,0.0,s,ns,ds,fs,0.0,(float)(nu[ir]-1),us);

                /* compute surface segments uniformly sampled in s */
                ss = alloc1(ns,sizeof(SurfaceSegment));
                for (is=0; is<ns; ++is) {
                        intcub(0,nu[ir],u,xud,1,&us[is],&rsx);
                        intcub(0,nu[ir],u,zud,1,&us[is],&rsz);
                        intcub(1,nu[ir],u,xud,1,&us[is],&rsxd);
                        intcub(1,nu[ir],u,zud,1,&us[is],&rszd);
                        ss[is].x = rsx;
                        ss[is].z = rsz;
                        ss[is].c = rsxd/sqrt(rsxd*rsxd+rszd*rszd);
                        ss[is].s = -rszd/sqrt(rsxd*rsxd+rszd*rszd);
                }

                /* fill in surface structure */
                rr[ir].ns = ns;
                rr[ir].ds = ds;
                rr[ir].ss = ss;

                /* free workspace */
                free1float(us);
                free1float(s);
                free1float(u);
                free1float((float*)xud);
                free1float((float*)zud);

                /* free space replaced by surface segments */
                free1(xu[ir]);
                free1(zu[ir]);
        }

        /* free space replaced by surface segments */
        free1(nu);
        free1(xu);
        free1(zu);
}
void yxtoxy (int nx, float dx, float fx, float y[], 
	int ny, float dy, float fy, float xylo, float xyhi, float x[])
/*****************************************************************************
Compute a regularly-sampled, monotonically increasing function x(y) from a 
regularly-sampled, monotonically increasing function y(x) by inverse linear 
interpolation.
******************************************************************************
Input:
nx		number of samples of y(x)
dx		x sampling interval; dx>0.0 is required
fx		first x
y		array[nx] of y(x) values; y[0] < y[1] < ... < y[nx-1] required
ny		number of samples of x(y)
dy		y sampling interval; dy>0.0 is required
fy		first y
xylo		x value assigned to x(y) when y is less than smallest y(x)
xyhi		x value assigned to x(y) when y is greater than largest y(x)

Output:
x		array[ny] of x(y) values
******************************************************************************
Notes:
User must ensure that:
(1) dx>0.0 && dy>0.0
(2) y[0] < y[1] < ... < y[nx-1]
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
{
	int nxi,nyo,jxi1,jxi2,jyo;
	float dxi,fxi,dyo,fyo,fyi,yo,xi1,yi1,yi2; 

	nxi = nx; dxi = dx; fxi = fx;
	nyo = ny; dyo = dy; fyo = fy;
	fyi = y[0];

	/* loop over output y less than smallest input y */
	for (jyo=0,yo=fyo; jyo<nyo; jyo++,yo+=dyo) {
		if (yo>=fyi) break;
		x[jyo] = xylo;
	}

	/* loop over output y between smallest and largest input y */
	if (jyo==nyo-1 && yo==fyi) {
		x[jyo++] = fxi;
		yo += dyo;
	}
	jxi1 = 0;
	jxi2 = 1;
	xi1 = fxi;
	while (jxi2<nxi && jyo<nyo) {
		yi1 = y[jxi1];
		yi2 = y[jxi2];
		if (yi1<=yo && yo<=yi2) {
			x[jyo++] = xi1+dxi*(yo-yi1)/(yi2-yi1);
			yo += dyo;
		} else {
			jxi1++;
			jxi2++;
			xi1 += dxi;
		}
	}

	/* loop over output y greater than largest input y */
	while (jyo<nyo)
		x[jyo++] = xyhi;
}
/* Estimation of the Z coordinate of the topography and its normal vector */
/*估计Z坐标的地形和法向量*/
  void zcoorTopog(float fxs,float dxs,int nxs,Surface *srf,float *sz,
		  float *nangl)
/*****************************************************************************
From the cubic spline calculation, the Z coor. of the surface are estimated
三次样条曲线的计算,Z坐标。表面的估计
******************************************************************************
Input:
fxs             x-coordinate of the first shot (travel-time tables)
dxs             x-coordinate increment of shots (travel-time tables)
nxs             number of shots (travel-time tables)
srf             surface structure

Output:
sz              z[] coordinates of the current surface (shot positions)
nangl           array of angles that the normal form with the vertical
*****************************************************************************
Author: Trino Salinas, CSM, 1996.
*****************************************************************************/
{
        int   ik,jx,ns,ixi,is,ix;
        float x,ax,ax0,az,temp;
        float *ssx,*ssz,*sss;
        SurfaceSegment *ss;

        ns = srf[0].ns;
        ss = srf[0].ss;
        ssx = alloc1float(ns);
        ssz = alloc1float(ns);
        sss = alloc1float(ns);

        /* number of segments, segment length */
        for (jx=0; jx<ns ; jx++) {
            ssx[jx] = ss[jx].x;
            ssz[jx] = ss[jx].z;
            sss[jx] = ss[jx].s;
        }
        ixi = 0;
        for (ix=0; ix<nxs; ++ix) {
            x = fxs + ix*dxs;
            for (ik=ixi; ik<ns-1; ++ik) {
                if (ssx[ik]<=x && ssx[ik+1]>=x) {
                   az = ssz[ik] - ssz[ik+1];
                   ax0 = ssx[ik+1] - x;
                   ax = ssx[ik+1] - ssx[ik];
                   sz[ix] = ax0*az/ax+ssz[ik+1];
                   az = sss[ik] - sss[ik+1];
                   temp = ax0*az/ax + sss[ik+1];
                   nangl[ix] = asin(temp);
                   ixi = ik;
                }
            }
        }
        sz[0] = 2*sz[1] - sz[2];
        nangl[0] = 2*nangl[1] - nangl[2];
        sz[nxs-1] = 2*sz[nxs-2] - sz[nxs-3];
        nangl[nxs-1] = 2*nangl[nxs-2] - nangl[nxs-3];

         for(is=0;is<nxs;is++)
           warn("sz[%d]=%g,nangl[%d]=%g",is,sz[is],is,nangl[is]);

        free1float(ssx);
        free1float(ssz);
        free1float(sss);
}

/* Prototypes of functions used interally */
 void dfrungk(float v,float vx,float vz,float vxx,float vxz,float vzz,
     	float x,float z,float px,float pz,float p2,float q2,
      	float *dx,float *dz,float *dpx,float *dpz,float *dp2,float *dq2);
 void sum2(float h,float a1,float a2,float a3,float a4,float a5,float a6,
	float b1,float b2,float b3,float b4,float b5,float b6,
	float *c1,float *c2,float *c3,float *c4,float *c5,float *c6);

void makeRay(Geo2d geo2dv,float *vel,float *vxx,float *vxz,float *vzz,
    Raytr raytr,float a0,int *nrs,Ray *ray,int npv,float *pv) 
/*****************************************************************************
 Trace rays in uniformly sampled v(nz,nx)
*****************************************************************************
Input:
geo2dv		grid parameters of velocity 
v2		sampled velocity array
raytr		geometry parameters of ray tracing 
vxx,vxz,vzz	sampled second derivatives of velocity  
a0  		shooting angles at source point
npv		flag of velocity analysis
pv		velocity variation if npv>0

Output:
ray		pointer to ray parameters sampled at discrete ray steps
nrs		number of points on ray 
*****************************************************************************
Note:
The ray ends when it runs out of time or with the first step that is out
of the velocity boudary or out of emergence angle boundary
*****************************************************************************
Author: Zhenyue Liu, CSM, 1995.
****************************************************************************/ 
{
	int nx=geo2dv.nx,nz=geo2dv.nz,nt=raytr.nt;
	float fx=geo2dv.fx,fz=geo2dv.fz;
	float dt=raytr.dt,amin=raytr.amin,amax=raytr.amax;

	int it;
	float tzmin,tzmax,ov,lx,lz,h,h2,h6;

	float v,dvdx,dvdz,uxx,uxz,uzz,tzt,dtv=0.0;
	float x,z,px,pz,p2,q2,
	      dx,dz,dpx,dpz,dp2,dq2;
	float xt,zt,pxt,pzt,p2t,q2t,
	    dxt,dzt,dpxt,dpzt,dp2t,dq2t;


	/* velocity and derivatives at source	*/
	x = raytr.xs;
	z = raytr.zs;
	vel2Interp(geo2dv,vel,vxx,vxz,vzz,x,z,&v,&dvdx,&dvdz,&uxx,&uxz,&uzz);
	ov = 1.0/v;

	p2 = 1.0;
	q2 = 0.0;
	*nrs = nt;

/*	compute slowness at source	*/
	px = ov*sin(a0);
	pz = ov*cos(a0);

/*	first ray step output	*/
	tzt = v*pz;
	ray[0].x = x;
	ray[0].z = z;
	ray[0].px = px;
	ray[0].pz = pz;
	ray[0].q2 = q2;
	ray[0].p2 = p2;
	ray[0].v = v;
	ray[0].dvdx = dvdx;
	ray[0].dvdz = dvdz;
	if(npv) ray[0].dtv = 0;

/*	compute minimum and maximum z-component of unit ray vector	*/
	tzmin = cos(amax)-0.01;
	tzmax = cos(amin)+0.01;

/*	determine fraction steps for Rung-Kuta	*/
	h = dt;		h2 = h*0.5;	h6 = h/6.0;

	lx = fx+(nx-1)*geo2dv.dx;
	lz = fz+(nz-1)*geo2dv.dz;

/*	loop over time steps	*/
	for(it=1; it<nt; ++it){
		if(x<fx || x>lx || z<fz || z>lz ||
      		   tzt<tzmin || tzt>tzmax) {
			it -= 1; 
			break;
		}

/*	step 1 of 4th-order Rung-Kuta	*/

/*	  compute k1 = f(y0)	*/
	    dfrungk(v,dvdx,dvdz,uxx,uxz,uzz,
      		x,z,px,pz,p2,q2,
      		&dx,&dz,&dpx,&dpz,&dp2,&dq2);
          
/*	  compute y1 = y0+0.5*h*k1	*/
	    sum2(h2,x,z,px,pz,p2,q2,
      		dx,dz,dpx,dpz,dp2,dq2,
      		&xt,&zt,&pxt,&pzt,&p2t,&q2t);

/*	velocity interpolation	*/
	vel2Interp(geo2dv,vel,vxx,vxz,vzz,
      	  xt,zt,&v,&dvdx,&dvdz,&uxx,&uxz,&uzz);

/*	step 2 of 4th-order Rung-Kuta	*/

/*	  compute k2 = f(y1)	*/
	    dfrungk(v,dvdx,dvdz,uxx,uxz,uzz,
      		xt,zt,pxt,pzt,p2t,q2t,
      		&dxt,&dzt,&dpxt,&dpzt,&dp2t,&dq2t);

/*	  compute y2 = y0+0.5*h*k2	*/
	    sum2(h2,x,z,px,pz,p2,q2,
      		dxt,dzt,dpxt,dpzt,dp2t,dq2t,
      		&xt,&zt,&pxt,&pzt,&p2t,&q2t);

/*	  compute k = k1+2.0*k2	*/
	    sum2(2.0,dx,dz,dpx,dpz,dp2,dq2,
      		dxt,dzt,dpxt,dpzt,dp2t,dq2t,
      		&dx,&dz,&dpx,&dpz,&dp2,&dq2);
 
/*	velocity interpolation	*/
	vel2Interp(geo2dv,vel,vxx,vxz,vzz,
      	  xt,zt,&v,&dvdx,&dvdz,&uxx,&uxz,&uzz);

/*	step 3 of 4th-order Rung-Kuta	*/

/*	  compute k3 = f(y2)	*/
	    dfrungk(v,dvdx,dvdz,uxx,uxz,uzz,
      		xt,zt,pxt,pzt,p2t,q2t,
      		&dxt,&dzt,&dpxt,&dpzt,&dp2t,&dq2t);

/*	  compute y3 = y0+h*k2	*/
 	    sum2(h,x,z,px,pz,p2,q2,
      		dxt,dzt,dpxt,dpzt,dp2t,dq2t,
      		&xt,&zt,&pxt,&pzt,&p2t,&q2t);

/*	  compute k = k1+2.0*k2+2.0*k3	*/
	    sum2(2.0,dx,dz,dpx,dpz,dp2,dq2,
      		dxt,dzt,dpxt,dpzt,dp2t,dq2t,
      		&dx,&dz,&dpx,&dpz,&dp2,&dq2);
 
/*	velocity interpolation	*/
	vel2Interp(geo2dv,vel,vxx,vxz,vzz,
      	  xt,zt,&v,&dvdx,&dvdz,&uxx,&uxz,&uzz);

/*	step 4 of 4th-order Rung-Kuta	*/

/*	  compute k4 = f(y3)	*/
	    dfrungk(v,dvdx,dvdz,uxx,uxz,uzz,
      		xt,zt,pxt,pzt,p2t,q2t,
      		&dxt,&dzt,&dpxt,&dpzt,&dp2t,&dq2t);

/*	  compute k = k1+2.0*k2+2.0*k3+k4	*/
	    sum2(1.0,dx,dz,dpx,dpz,dp2,dq2,
      		dxt,dzt,dpxt,dpzt,dp2t,dq2t,
      		&dx,&dz,&dpx,&dpz,&dp2,&dq2);
 
/*	  compute y4 = y0+h*k/6	*/
 	    sum2(h6,x,z,px,pz,p2,q2,
      		dx,dz,dpx,dpz,dp2,dq2,
      		&x,&z,&px,&pz,&p2,&q2);
         
/*	velocity interpolation	*/
	    vel2Interp(geo2dv,vel,vxx,vxz,vzz,
      	  	x,z,&v,&dvdx,&dvdz,&uxx,&uxz,&uzz);

/*	velocity itself interpolation	*/
	if(npv)	vintp(geo2dv,pv,x,z,&dtv);

		tzt = v*pz;

/*	save ray parameters */
		ray[it].x = x;
		ray[it].z = z;
		ray[it].px = px;
		ray[it].pz = pz;
		ray[it].p2 = p2;
		ray[it].q2 = q2;
		ray[it].v = v;
		ray[it].dvdx = dvdx;
		ray[it].dvdz = dvdz;
		if(npv)  ray[it].dtv = dtv/v;
	}
	*nrs = it; 

}

/***********************************************************************/
  void dfrungk(float v,float vx,float vz,float vxx,float vxz,float vzz,
     	float x,float z,float px,float pz,float p2,float q2,
      	float *dx,float *dz,float *dpx,float *dpz,float *dp2,float *dq2)
{

	float ov,vv,c,s,v11;

		ov = 0.0*z+ 0.0*x + 1.0/v;
		vv = v*v;
		*dx = vv*px;
		*dz = vv*pz;
		*dpx = -ov*vx;
		*dpz = -ov*vz;

		c = v*pz;
		s = v*px;
		v11 = vxx*c*c-2.0*vxz*c*s+vzz*s*s;
		*dp2 = -ov*v11*q2;
		*dq2 = vv*p2;

}


/************************ c = a+h*b ***********************************/
 void sum2(float h,float a1,float a2,float a3,float a4,float a5,float a6,
	float b1,float b2,float b3,float b4,float b5,float b6,
	float *c1,float *c2,float *c3,float *c4,float *c5,float *c6)
{

	*c1 = a1+h*b1;
	*c2 = a2+h*b2;
	*c3 = a3+h*b3;
	*c4 = a4+h*b4;
	*c5 = a5+h*b5;
	*c6 = a6+h*b6;

}

/* Prototype of function used internally */
void ddt(float p,float q,float c,float s, float *dv,float v,float *d2t,
	float *cuv);

void raytime2d(Raytr raytr,Geo2d geo2dv,float *v,Geo2d geo2dt,float *t,
	int npv,float *vo,float *pv,float *tv,float *cs,FILE *fpraypath) 
/****************************************************************************
 raytime2d - calculate traveltimes by paraxial ray tracing
*****************************************************************************
Input:
geo2dt		grid parameters of traveltime 
geo2dv		grid parameters of velocity 
v2		sampled velocity array
raytr		geometry parameters of ray tarcing 

Output:
t		traveltime 
*****************************************************************************
Note: when traveltime has multiple values, the one has the shortest 
 ray path is chosen.
*****************************************************************************
Author: Zhenyue Liu, CSM 1995.
****************************************************************************/ 
{
    

	int na=raytr.na,nt=raytr.nt,nxo=geo2dt.nx,nzo=geo2dt.nz,
		nx=geo2dv.nx,nz=geo2dv.nz;
	float fxo=geo2dt.fx,fzo=geo2dt.fz,dxo=geo2dt.dx,dzo=geo2dt.dz,
		odxo=geo2dt.odx,odzo=geo2dt.odz;
	float fac=raytr.fac,dt=raytr.dt,fa=raytr.fa,da=raytr.da,xs=raytr.xs;
	float *vxx,*vxz,*vzz,*s,zs=raytr.zs;
	int nrs;
	Ray *ray;

 
 	int i,ia,ixo,izo;
	float a,xo,zo,exo,ezo;

        /*variables used for extrapolation
            用于外推的变量     */
	int it,jt,nxf,nxe,nzf,nze;
	float tc,xoc,zoc,xc,zc,r1,v0=0.0,norm2,terr,vd1,cuv,
      		sd,sind,cosd,vd,pxd,pzd,r2,t1,t2,r1x,r2x,t1x,t2x,t2xz,
      		p2d,q2d,gradv[2],d2t[3],dtvd,odt=1.0/dt,xcosd=0.0,*tvd;

	vxx = alloc1float(nx*nz);
	vxz = alloc1float(nx*nz);
	vzz = alloc1float(nx*nz);
	tvd = alloc1float(nt);
	s = alloc1float(nxo*nzo);
	ray = (Ray*)alloc1(nt,sizeof(Ray)); 

	/* compute second derivatives of velocity
             计算速度的二阶导数    */
 	dv2(geo2dv,v,vxx,vxz,vzz); 


       /*maximum range of output points
           输出点的最大范围       */
	exo = fxo+(nxo-1)*dxo;
	ezo = fzo+(nzo-1)*dzo;

      /*initialize traveltime and raypath length
           初始化走时和射线路径长度        */
	for(ixo=0; ixo<nxo; ++ixo) 
		for(izo=0; izo<nzo; ++izo){
			i = izo+ixo*nzo;
			t[i] = INITIAL_T;
			s[i] = INITIAL_T;
		     if(npv) tv[i] = cs[i] = 0.0;
	}

/* 	loop over shooting angles at source point 	*/
	for(ia=0,a=fa; ia<na; ++ia,a+=da){  
	
/*		trace rays	*/
		makeRay(geo2dv,v,vxx,vxz,vzz,raytr,a,&nrs,ray,npv,pv);

/*		extropolate to grids near central rays	*/
		    v0 = ray[0].v;
		    r2 = 0.0;
		    xc = raytr.xs;
		    zc = raytr.zs;
		    sd = 0;
		    vd1 = v0;
		    if(npv) {
		      dtvd = 0.0;
		      for(it=1; it<nrs; ++it) {
			dtvd += 0.5*dt*(ray[it].dtv+ray[it-1].dtv);
			tvd[it] = ray[it].v*dtvd;
		      }
		    }
		    for (it=1; it<nrs; ++it) {
			xo = ray[it].x;
			zo = ray[it].z;
			vd = ray[it].v;
			sd = sd+0.5*dt*(vd+vd1);
			vd1 = vd;

                     /*if the point is within the previous circle
                         如果点是在前面的圆圈           */
			if(r2 > (xc-xo)*(xc-xo)+(zc-zo)*(zc-zo)
      			  && it != nrs-1) continue;

                     /*if the point is out of output range
                         如果点超出输出范围       */
			if(xo<fxo || xo>exo || zo<fzo || zo>ezo) continue;
								
			xc = xo;
			zc = zo;
			q2d = ray[it].q2; 

                     /*if caustics*/
			if(q2d == 0.0) {
				r2 = 0.0;
				continue;
			}

			vd = ray[it].v;
			pxd = ray[it].px;
			pzd = ray[it].pz;

			p2d = ray[it].p2;
			sind = vd*pxd;
			cosd = vd*pzd;
			gradv[0] = ray[it].dvdx;
			gradv[1] = ray[it].dvdz;

                     /*calculate second derivatives of traveltime
                           计算走时的二阶导数      */
			ddt(p2d,q2d,cosd,sind,gradv,vd,d2t,&cuv);

                     /*determine radius for extrapolation
                         确定外推的半径      */
			tc = it*dt;
			terr = tc*fac;
			norm2 = sqrt(d2t[0]*d2t[0]+d2t[2]*d2t[2]+
      				2.0*d2t[1]*d2t[1]);

			r2 = terr/norm2;
			r1 = sqrt(r2);
                    /*keep ray-centered coordinate system regular 
                          保持射线中心坐标系统正常       */
			if(r1*cuv > 0.1) r1 = 0.1/cuv;

                    /*radius cannot be too large
                          半径不能太大      */
			if(r1 > 0.1*sd) r1 = 0.1*sd;

			r2 = r1*r1;

			nxf = (xc-r1-fxo)*odxo+0.9;
			if(nxf < 0) nxf = 0;
			nxe = (xc+r1-fxo)*odxo+0.1;
			if(nxe >= nxo) nxe = nxo-1;

                   /* output raypath*/
                     fprintf(fpraypath,"%g %g %g\n",
                                       xc,zc,tc); 

			for(ixo=nxf; ixo<=nxe; ++ixo){
				xoc = fxo-xc+ixo*dxo;
				r2x = r2-xoc*xoc;
				if(r2x < 0) continue;

				r1x = sqrt(r2x);
				t1x = tc+xoc*pxd;
				t2x = tc*xoc*xoc*d2t[0];
				t2xz = 2.0*xoc*d2t[1];
				if(npv) xcosd = pzd+xoc*d2t[1];

				nzf = (zc-r1x-fzo)*odzo+0.9;
				if(nzf < 0) nzf = 0;
				nze = (zc+r1x-fzo)*odzo+0.1;
				if(nze>=nzo) nze = nzo-1;

				i = ixo*nzo;
				for(izo=nzf; izo<=nze; ++izo){ 
                                    /*if ray path is shorter
                                        如果射线路径更短   */
				    if(sd < s[i+izo] ) {
					zoc = fzo-zc+izo*dzo;
					t1 = t1x+zoc*pzd;
					t2 = t2x+zoc*(zoc*d2t[2]+t2xz);
					s[i+izo] = sd;
					t[i+izo] = t1*t1+tc*t2;
				      if(npv) {
					jt = (t1+0.5*t2)*odt+0.5;
					if(jt<0) jt = 0;
					if(jt>nrs-1) jt = nrs-1;
					tv[i+izo] = tvd[jt];
					cs[i+izo] = xcosd+d2t[2]*zoc;
				      }
				    }
 				}
 			}
 		}

	}

      /*square root of traveltimes 
         走时开方  */
	for (ixo=0; ixo<nxo; ++ixo){ 
		i = ixo*nzo;
		for (izo=0; izo<nzo; ++izo){  
			t[i+izo] = MAX(0.0,t[i+izo]);
			if(t[i+izo] < 999) t[i+izo] = sqrt(t[i+izo]);
  			if(npv) cs[i+izo] *= vo[i+izo];
		}
	}

       /*compute traveltime near source
          计算震源附近走时       */
	nxf = (xs-fxo)/dxo-1.5;
	if(nxf<0) nxf = 0;
	nxe = (xs-fxo)/dxo+2.5;
	if(nxe>nxo-1) nxe = nxo-1;
	nzf = (zs-fzo)/dzo-0.5;
	if(nzf<0) nzf = 0;
	nze = (zs-fzo)/dzo+1.5;
	if(nze>=nzo) nze = nzo-1;
	for (ixo=nxf; ixo<=nxe; ++ixo){
		xo = fxo+ixo*dxo;
		i = ixo*nzo;
		for (izo=nzf; izo<=nze; ++izo){
			zo = fzo+izo*dzo;
			if(t[i+izo] > 999) 
			  t[i+izo] = sqrt((xo-xs)*(xo-xs)+(zo-zs)*(zo-zs))/v0;
 		}
	}

	free1float(vxx);
	free1float(vxz);
	free1float(vzz);
	free1float(s);
	free1((void*)ray);

}

void ddt(float p,float q,float c,float s, float *dv,float v,float *d2t,
	float *cuv)
/************************************************************************
ddt - calculate second derivatives of traveltime with respect to x,y
************************************************************************/
{
 	float ov2,m11,m12,m22;

	ov2 = 1.0/(v*v);
	m11 = p/q;

/*	calculate 2nd column of m	*/
	m12 = -(dv[0]*c-dv[1]*s)*ov2;
	m22 = -(dv[0]*s+dv[1]*c)*ov2;

/*	compute h*m*h^T	*/
	d2t[0] = m11*c*c+2.0*m12*c*s+m22*s*s;
	d2t[1] = (m22-m11)*c*s+m12*(c*c-s*s);
	d2t[2] = m11*s*s-2.0*m12*c*s+m22*c*c;

/*	compute the curvature of raypath	*/
	*cuv = fabs(m12)*v;

}


void trans(int nx,int nz,int nxt,int nx0,float *v,float *vt) 
{
	int ix,iz,i,i0;

	for(ix=0*nx; ix<nxt; ++ix){
		i = ix*nz;
		i0 = (ix+nx0)*nz;
		for(iz=0; iz<nz; ++iz)
			vt[i+iz] = v[i0+iz];
	}
}

void voint(Geo2d geo2dv,float *v,Geo2d geo2dt,float *ov2,int npv,float *vo)
/*interpolate velocity and compute slowness squares 
           插入速度和计算慢度场*/
{

	int nx=geo2dv.nx,nz=geo2dv.nz,nxo=geo2dt.nx,nzo=geo2dt.nz;
	int i,io,jx,jz,ixo,izo;
	float fx=geo2dv.fx,fz=geo2dv.fz,odx=geo2dv.odx,odz=geo2dv.odz;
	float fxo=geo2dt.fx,fzo=geo2dt.fz,dxo=geo2dt.dx,dzo=geo2dt.dz;
	float x,z,ax,az,sx,sz,temp;


	for(ixo=0,x=fxo; ixo<nxo; ++ixo,x+=dxo){
		ax = (x-fx)*odx;
		jx = ax;
		sx = ax-jx;
		if(jx < 0) jx = 0;
		if(jx > nx-2) jx = nx-2;
		for(izo=0,z=fzo; izo<nzo; ++izo,z+=dzo){
			az = (z-fz)*odz;
			jz = az;
			sz = az-jz;
			if(jz < 0) jz = 0;
			if(jz > nz-2) jz = nz-2;
			io = ixo*nzo+izo;
			i = jx*nz+jz;
		    	temp = (1.0-sx)*((1.0-sz)*v[i]+sz*v[i+1])
      			    	+sx*((1.0-sz)*v[i+nz]+sz*v[i+nz+1]);
			if(npv) vo[io] = temp;
			ov2[io] = 1.0/(temp*temp);
		}
	}
}

 
void dv2(Geo2d geo2dv,float *v,float *vxx,float *vxz,float *vzz) 
/*calculate second DeriVatives from a 2D VELocity grid by finite difference*/
/*  	input: 	velocity v 
  	output: vxx,vxz and vzz  */	 
{ 
  	int  nx=geo2dv.nx, nz=geo2dv.nz, ix, iz, i; 	 
 	float   odx=geo2dv.odx,odz=geo2dv.odz,odxx,odzz,odxz; 
   	
	odxx = odx*odx; 
	odzz = odz*odz; 
	odxz = 0.25*odx*odz; 

/*	initialize	*/
	for(ix=0; ix<nx; ++ix)    
 	    for(iz=0; iz<nz; ++iz){
 		i = ix*nz+iz;
		vxx[i] = 0.; 
		vxz[i] = 0.; 
		vzz[i] = 0.; 
	}
	
/*	finite difference 	*/
	for(ix=1; ix<nx-1; ++ix)    
 	    for(iz=1; iz<nz-1; ++iz){
 		i = ix*nz+iz;
 		vxx[i] = odxx*(v[i+nz]-2.0*v[i]+v[i-nz]); 
 		vxz[i] = odxz*(v[i+nz+1]-v[i+nz-1]-v[i-nz+1]+v[i-nz-1]); 
      		vzz[i] = odzz*(v[i+1]-2.0*v[i]+v[i-1]);
	}

}

void eiknl(Geo2d geo2dt,float *t,float *ov2,float tmax) 
/* 	compute traveltime in shadow zones by eikonal equation
            程函方程计算走时在阴影区
  input
   t 		traveltimes form ray tracing
   ov2 		slowness squares
   tmax	 	upper limit of ordinary traveltime value
  output:
   t	traveltime (unchanged if t<=tmax)
*/
{
	int nx=geo2dt.nx,nz=geo2dt.nz,ix,iz;
	float dx=geo2dt.dx,dz=geo2dt.dz;
	float tx2,tz2,t0,tl,tr,temp,odx2,*tt1,*tt2;

	tt1 = alloc1float(nx);
	tt2 = alloc1float(nx);

	odx2 = 1.0/(dx*dx);

	for(ix=0; ix<nx; ++ix)
		tt1[ix] = t[ix*nz];

	for(iz=1; iz<nz; ++iz){ 
	    for(ix=0; ix<nx; ++ix) 
		tt2[ix] = t[ix*nz+iz];

	    for(ix=1; ix<nx; ++ix){
        /*if traveltime is abnormal and the upper is normal*/
		t0 = tt1[ix];
		if(tt2[ix] > tmax && t0 <= tmax) {
			tl = tr = 0.0; 
			if(ix > 0) tl = tt1[ix-1];
			if(ix < nx-1) tr = tt1[ix+1];
                        tx2 = (tl-t0)*(tl-t0);
                        temp = (tr-t0)*(tr-t0);
                        if(tx2>temp) tx2 = temp;
                        temp = 0.25*(tl-tr)*(tl-tr);
                        if(tx2>temp) tx2 = temp;
                        tx2 *= odx2;

			tz2 = 0.5*(ov2[ix*nz+iz-1]+ov2[ix*nz+iz])-tx2;
			if(tz2 >= 0) tt2[ix] = t0+dz*sqrt(tz2);
		}
	    }
	    for(ix=0; ix<nx; ++ix){
		t[ix*nz+iz] = tt2[ix];
		tt1[ix] = tt2[ix];
	    }
	}

	free1float(tt1);
	free1float(tt2);

}

void vel2Interp(Geo2d geo2dv,float *v,float *vxx,float *vxz,float *vzz,
     float x,float z,
     float *u,float *ux,float *uz,float *uxx,float *uxz,float *uzz) 
/*************************************************************************
vel2Interp - Function to support interpolation of velocity and its
		derivatives	 
*************************************************************************
Input:
x,z		2D coordinate at which to interpolate 
v		velocity array(nz,nx) on unoform grids 
vxx,vxz,vzz	second derivaitve arrays(nz,nx) on uniform grids 

Output:  
u 	v(x,z)
ux  	dv/dx
uz  	dv/dz
uxx	ddv/dxdx
uxz 	ddv/dzdx
uzz	ddv/dzdz
*************************************************************************
Author: Zhenyue Liu, CSM 1995
*************************************************************************/
{
	int k0,k1,k2,k3,jx,jz,nx=geo2dv.nx,nz=geo2dv.nz;
	float fx=geo2dv.fx,dx=geo2dv.dx,odx=geo2dv.odx;
	float fz=geo2dv.fz,dz=geo2dv.dz,odz=geo2dv.odz;
	float ax,az,sx,sz,sxx,szz,a0,a1,a2,a3;
	float g0,g1,g2,g3,gx0,gx1,gx2,gx3,gz0,gz1,gz2,gz3;

/*	determine interpolate coefficients	*/
	    ax = (x-fx)*odx;
	    jx = ax;
	    if(jx<0) jx = 0;
	    if(jx>nx-2) jx = nx-2;
	    sx = ax-jx;
	    az = (z-fz)*odz;
	    jz = az;
	    if(jz<0) jz = 0;
	    if(jz>nz-2) jz = nz-2;
	    sz = az-jz;

	    sxx = 0.5*sx*(1.0-sx)*dx*dx;
	    szz = 0.5*sz*(1.0-sz)*dz*dz;

	    a0 = (1.0-sx)*(1.0-sz);
	    a1 = (1.0-sx)*sz;
	    a2 = sx*(1.0-sz);
	    a3 = sx*sz;

/*	    set the table of indices for interpolation	*/
	    k0 = nz*jx+jz;
	    k1 = k0+1;
 	    k2 = k0+nz;
	    k3 = k2+1;

	    g0 = v[k0];
	    g1 = v[k1];
	    g2 = v[k2];
	    g3 = v[k3];
	    gx0 = vxx[k0];
	    gx1 = vxx[k1];
	    gx2 = vxx[k2];
	    gx3 = vxx[k3];
	    gz0 = vzz[k0];
	    gz1 = vzz[k1];
	    gz2 = vzz[k2];
	    gz3 = vzz[k3];

/*	interpolation	*/
	    *uxx = a0*gx0+a1*gx1+a2*gx2+a3*gx3;
	    *uxz = a0*vxz[k0]+a1*vxz[k1]+a2*vxz[k2]+a3*vxz[k3];
	    *uzz = a0*gz0+a1*gz1+a2*gz2+a3*gz3;

	    *u = a0*g0+a1*g1+a2*g2+a3*g3-(sxx*(*uxx)+szz*(*uzz));
	    *ux = ((1.0-sz)*(g2-g0-sxx*(gx2-gx0)-szz*(gz2-gz0))
      		    +sz*(g3-g1-sxx*(gx3-gx1)-szz*(gz3-gz1)))*odx
      		    +(sx-0.5)*dx*(*uxx);
	    *uz = ((1.0-sx)*(g1-g0-sxx*(gx1-gx0)-szz*(gz1-gz0))
      		    +sx*(g3-g2-sxx*(gx3-gx2)-szz*(gz3-gz2)))*odz
      		    +(sz-0.5)*dz*(*uzz);
}

void vintp(Geo2d geo2dv,float *v,float x,float z,float *u) 
/*************************************************************************
 Function to support interpolation of a single fuction 
*************************************************************************
Input:
x,z		2D coordinate at which to interpolate
v		array(nz,nx) on unoform grids 
Output: 
u		v(x,z)
*************************************************************************
Author: Zhenyue Liu, CSM 1995
*************************************************************************/
{

	int k0,k1,k2,k3,jx,jz,nx=geo2dv.nx,nz=geo2dv.nz;
	float fx=geo2dv.fx,odx=geo2dv.odx;
	float fz=geo2dv.fz,odz=geo2dv.odz;
	float ax,az,sx,sz;

/*	determine interpolate coefficients	*/
	    ax = (x-fx)*odx;
	    jx = ax;
	    if(jx<0) jx = 0;
	    if(jx>nx-2) jx = nx-2;
	    sx = ax-jx;
	    az = (z-fz)*odz;
	    jz = az;
	    if(jz<0) jz = 0;
	    if(jz>nz-2) jz = nz-2;
	    sz = az-jz;

/*	    set the table of indices for interpolation	*/
	    k0 = nz*jx+jz;
	    k1 = k0+1;
 	    k2 = k0+nz;
	    k3 = k2+1;

/*	interpolation	*/
	    *u = (1.0-sx)*((1.0-sz)*v[k0]+sz*v[k1])
      			+sx*((1.0-sz)*v[k2]+sz*v[k3]);

}
/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/*****************************************************************************
INTCUB - evaluate y(x), y'(x), y''(x), ... via piecewise cubic interpolation

intcub		evaluate y(x), y'(x), y''(x), ... via piecewise spline
			interpolation

******************************************************************************
Function Prototype:
void intcub (int ideriv, int nin, float xin[], float ydin[][4],

******************************************************************************
Input:
ideriv		=0 if y(x) desired; =1 if y'(x) desired, ...
nin		length of xin and ydin arrays
xin		array[nin] of monotonically increasing or decreasing x values
ydin		array[nin][4] of y(x), y'(x), y''(x), and y'''(x)
nout		length of xout and yout arrays
xout		array[nout] of x values at which to evaluate y(x), y'(x), ...

Output:
yout		array[nout] of y(x), y'(x), ... values

******************************************************************************
Notes:
xin values must be monotonically increasing or decreasing.

Extrapolation of the function y(x) for xout values outside the range
spanned by the xin values is performed using the derivatives in 
ydin[0][0:3] or ydin[nin-1][0:3], depending on whether xout is closest
to xin[0] or xin[nin-1], respectively.

******************************************************************************
Reference:
See: Greville, T. N., Theory and Application of Spline Functions, Academic Press
for a general discussion of spline functions.

******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
/**************** end self doc ********************************/


#define OO2 0.5000000
#define OO6 0.1666667

void intcub (int ideriv, int nin, float xin[], float ydin[][4], 
	int nout, float xout[], float yout[])
/*****************************************************************************
evaluate y(x), y'(x), y''(x), ... via piecewise cubic interpolation
******************************************************************************
Input:
ideriv		=0 if y(x) desired; =1 if y'(x) desired, ...
nin		length of xin and ydin arrays
xin		array[nin] of monotonically increasing or decreasing x values
ydin		array[nin][4] of y(x), y'(x), y''(x), and y'''(x)
nout		length of xout and yout arrays
xout		array[nout] of x values at which to evaluate y(x), y'(x), ...

Output:
yout		array[nout] of y(x), y'(x), ... values
******************************************************************************
Notes:
xin values must be monotonically increasing or decreasing.

Extrapolation of the function y(x) for xout values outside the range
spanned by the xin values is performed using the derivatives in 
ydin[0][0:3] or ydin[nin-1][0:3], depending on whether xout is closest
to xin[0] or xin[nin-1], respectively.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
{
    void xindex (int nx, float ax[], float x, int *index);

	static int idx;
	int iout;
	float delx;

	/* y(x) is desired, then */
	if (ideriv==0) {
		for (iout=0; iout<nout; iout++) {
			xindex(nin,xin,xout[iout],&idx);
			delx = xout[iout]-xin[idx];
			yout[iout] = (ydin[idx][0]+delx*
				(ydin[idx][1]+delx*
				(ydin[idx][2]*OO2+delx*
				(ydin[idx][3]*OO6))));
		}
		
	/* else, if y'(x) is desired, then */
	} else if (ideriv==1) {
		for (iout=0; iout<nout; iout++) {
			xindex(nin,xin,xout[iout],&idx);
			delx = xout[iout]-xin[idx];
			yout[iout] = (ydin[idx][1]+delx*
				(ydin[idx][2]+delx*
				(ydin[idx][3]*OO2)));
		}
		
	/* else, if y''(x) is desired, then */
	} else if (ideriv==2) {
		for (iout=0; iout<nout; iout++) {
			xindex(nin,xin,xout[iout],&idx);
			delx = xout[iout]-xin[idx];
			yout[iout] = (ydin[idx][2]+delx*
				(ydin[idx][3]));
		}
		
	/* else, if y'''(x) is desired, then */
	} else if (ideriv==3) {
		for (iout=0; iout<nout; iout++) {
			xindex(nin,xin,xout[iout],&idx);
			delx = xout[iout]-xin[idx];
			yout[iout] = (ydin[idx][3]);
		}
		
	/* else, if any other derivative is desired, then */
	} else {
		for (iout=0; iout<nout; iout++)
			yout[iout] = 0.0;
	}
}
/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/*****************************************************************************
XINDEX - determine index of x with respect to an array of x values

xindex		determine index of x with respect to an array of x values

******************************************************************************
Input:
nx		number of x values in array ax
ax		array[nx] of monotonically increasing or decreasing x values
x		the value for which index is to be determined
index		index determined previously (used to begin search)

Output:
index		for monotonically increasing ax values, the largest index
		for which ax[index]<=x, except index=0 if ax[0]>x;
		for monotonically decreasing ax values, the largest index
		for which ax[index]>=x, except index=0 if ax[0]<x

******************************************************************************
Notes:
This function is designed to be particularly efficient when called
repeatedly for slightly changing x values; in such cases, the index 
returned from one call should be used in the next.

******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 12/25/89
*****************************************************************************/
/**************** end self doc ********************************/

void xindex (int nx, float ax[], float x, int *index)
/*****************************************************************************
determine index of x with respect to an array of x values
******************************************************************************
Input:
nx		number of x values in array ax
ax		array[nx] of monotonically increasing or decreasing x values
x		the value for which index is to be determined
index		index determined previously (used to begin search)

Output:
index		for monotonically increasing ax values, the largest index
		for which ax[index]<=x, except index=0 if ax[0]>x;
		for monotonically decreasing ax values, the largest index
		for which ax[index]>=x, except index=0 if ax[0]<x
******************************************************************************
Notes:
This function is designed to be particularly efficient when called
repeatedly for slightly changing x values; in such cases, the index 
returned from one call should be used in the next.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 12/25/89
*****************************************************************************/
{
	int lower,upper,middle,step;

	/* initialize lower and upper indices and step */
	lower = *index;
	if (lower<0) lower = 0;
	if (lower>=nx) lower = nx-1;
	upper = lower+1;
	step = 1;

	/* if x values increasing */
	if (ax[nx-1]>ax[0]) {

		/* find indices such that ax[lower] <= x < ax[upper] */
		while (lower>0 && ax[lower]>x) {
			upper = lower;
			lower -= step;
			step += step;
		}
		if (lower<0) lower = 0;
		while (upper<nx && ax[upper]<=x) {
			lower = upper;
			upper += step;
			step += step;
		}
		if (upper>nx) upper = nx;

		/* find index via bisection */
		while ((middle=(lower+upper)>>1)!=lower) {
			if (x>=ax[middle])
				lower = middle;
			else
				upper = middle;
		}

	/* else, if not increasing */
	} else {

		/* find indices such that ax[lower] >= x > ax[upper] */
		while (lower>0 && ax[lower]<x) {
			upper = lower;
			lower -= step;
			step += step;
		}
		if (lower<0) lower = 0;
		while (upper<nx && ax[upper]>=x) {
			lower = upper;
			upper += step;
			step += step;
		}
		if (upper>nx) upper = nx;

		/* find index via bisection */
		while ((middle=(lower+upper)>>1)!=lower) {
			if (x<=ax[middle])
				lower = middle;
			else
				upper = middle;
		}
	}

	/* return lower index */
	*index = lower;
}

void csplin (int n, float x[], float y[], float yd[][4])
/*****************************************************************************
compute cubic spline interpolation coefficients for interpolation with
continuous second derivatives
******************************************************************************
Input:
n		number of samples
x  		array[n] of monotonically increasing or decreasing abscissae
y		array[n] of ordinates

Output:
yd		array[n][4] of cubic interpolation coefficients (see notes)
******************************************************************************
Notes:
The computed cubic spline coefficients are as follows:
yd[i][0] = y(x[i])    (the value of y at x = x[i])
yd[i][1] = y'(x[i])   (the 1st derivative of y at x = x[i])
yd[i][2] = y''(x[i])  (the 2nd derivative of y at x = x[i])
yd[i][3] = y'''(x[i]) (the 3rd derivative of y at x = x[i])

To evaluate y(x) for x between x[i] and x[i+1] and h = x-x[i],
use the computed coefficients as follows:
y(x) = yd[i][0]+h*(yd[i][1]+h*(yd[i][2]/2.0+h*yd[i][3]/6.0))
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 10/03/89
Modified:  Dave Hale, Colorado School of Mines, 02/28/91
	changed to work for n=1.
Modified:  Dave Hale, Colorado School of Mines, 08/04/91
	fixed bug in computation of left end derivative
*****************************************************************************/
{
	int i;
	float h1,h2,del1,del2,dmax,hsum,w1,w2,divdf3,sleft,sright,alpha,t;

	/* if n=1, then use constant interpolation */
	if (n==1) {
		yd[0][0] = y[0];
		yd[0][1] = 0.0;
		yd[0][2] = 0.0;
		yd[0][3] = 0.0;
		return;

	/* else, if n=2, then use linear interpolation */
	} else if (n==2) {
		yd[0][0] = y[0];  yd[1][0] = y[1];
		yd[0][1] = yd[1][1] = (y[1]-y[0])/(x[1]-x[0]);
		yd[0][2] = yd[1][2] = 0.0;
		yd[0][3] = yd[1][3] = 0.0;
		return;
	}
	
	/* set left end derivative via shape-preserving 3-point formula */
	h1 = x[1]-x[0];
	h2 = x[2]-x[1];
	hsum = h1+h2;
	del1 = (y[1]-y[0])/h1;
	del2 = (y[2]-y[1])/h2;
	w1 = (h1+hsum)/hsum;
	w2 = -h1/hsum;
	sleft = w1*del1+w2*del2;
	if (sleft*del1<=0.0)
		sleft = 0.0;
	else if (del1*del2<0.0) {
		dmax = 3.0*del1;
		if (ABS(sleft)>ABS(dmax)) sleft = dmax;
	}

	/* set right end derivative via shape-preserving 3-point formula */
	h1 = x[n-2]-x[n-3];
	h2 = x[n-1]-x[n-2];
	hsum = h1+h2;
	del1 = (y[n-2]-y[n-3])/h1;
	del2 = (y[n-1]-y[n-2])/h2;
	w1 = -h2/hsum;
	w2 = (h2+hsum)/hsum;
	sright = w1*del1+w2*del2;
	if (sright*del2<=0.0)
		sright = 0.0;
	else if (del1*del2<0.0) {
		dmax = 3.0*del2;
		if (ABS(sright)>ABS(dmax)) sright = dmax;
	}
	
	/* compute tridiagonal system coefficients and right-hand-side */	
	yd[0][0] = 1.0;
	yd[0][2] = 2.0*sleft;
	for (i=1; i<n-1; i++) {
		h1 = x[i]-x[i-1];
		h2 = x[i+1]-x[i];
		del1 = (y[i]-y[i-1])/h1;
		del2 = (y[i+1]-y[i])/h2;
		alpha = h2/(h1+h2);
		yd[i][0] = alpha;
		yd[i][2] = 3.0*(alpha*del1+(1.0-alpha)*del2);
	}
	yd[n-1][0] = 0.0;
	yd[n-1][2] = 2.0*sright;
	
	/* solve tridiagonal system for slopes */
	t = 2.0;
	yd[0][1] = yd[0][2]/t;
	for (i=1; i<n; i++) {
		yd[i][3] = (1.0-yd[i-1][0])/t;
		t = 2.0-yd[i][0]*yd[i][3];
		yd[i][1] = (yd[i][2]-yd[i][0]*yd[i-1][1])/t;
	}
	for (i=n-2; i>=0; i--)
		yd[i][1] -= yd[i+1][3]*yd[i+1][1];

	/* copy ordinates into output array */
	for (i=0; i<n; i++)
		yd[i][0] = y[i];

	/* compute 2nd and 3rd derivatives of cubic polynomials */
	for (i=0; i<n-1; i++) {
		h2 = x[i+1]-x[i];
		del2 = (y[i+1]-y[i])/h2;
		divdf3 = yd[i][1]+yd[i+1][1]-2.0*del2;
		yd[i][2] = 2.0*(del2-yd[i][1]-divdf3)/h2;
		yd[i][3] = (divdf3/h2)*(6.0/h2);
	}
	yd[n-1][2] = yd[n-2][2]+(x[n-1]-x[n-2])*yd[n-2][3];
	yd[n-1][3] = yd[n-2][3];
}
