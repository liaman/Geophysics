/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* RAYT2DAN -- RAYTracing for P- and SV-waves in a 2D ANisotropic media*/
#include "hc/cjbsegy.h"
#include "hc/fft.c"
#include "hc/alloc.c"
#include "hc/complex.c"

/*//#include "ealloc.c"
//#include "par.h"
//#include "getpars.c"
//#include "intcub.c"
//#include "cubicspline.c"
//#include "xindex.c"
//#include "cwp.h"*/

//#include "par.h"

/*********************** self documentation **********************/
char *sdoc[] = {
"									",
" RAYT2DAN -- P- and SV-wave raytracing in 2D anisotropic media		",
"									",
" rayt2dan > ttime parameterfiles= nt= nx= nz= [optional parameters]	",
"									",
" Required Parameters:							",
" VP0file=		 name of file containing VP0(x,z)		",
" nt=                    number of time samples	for each ray		",
" nx=                    number of samples (x) for the parameter fields	",
" nz=                    number of samples (z) for the parameter fields	",
" 									",
" Optional Parameters:							",
" SV=0			 for P-waves, SV=1 for Shear waves		",
"									",
" Parameters defining the velocity field				",
" dt=0.008                 time sampling interval			",	
" fx=0                     first lateral sample (x) in parameter field  ",
" fz=0                     first lateral sample (z) in parameter field  ",
" dx=100.0                 sample spacing (x) for the parameter fields	",
" dz=100.0                 sample spacing (z) for the parameter fields	",
" 									",	
" Parameters defining the takeoff angle of a ray at a source position	",	
" fa=-60                 first take-off angle of rays (degrees)         ",
" na=61                  number of rays					",      
" da=2                   increment of take-off angle			",
" amin=0                 minimum emergence angle; must be > -90 degrees	",
" amax=90                maximum emergence angle; must be < 90 degrees	",
" 									",
" Parameters defining the output traveltime table			",
" fxo=fx                 first lateral sample in traveltime table	",
" nxo=nx                 number of later samples in traveltime table	",
" dxo=dx                 lateral interval in traveltime table		",
" fzo=fz                 first depth sample in traveltime table		",
" nzo=nz                 number of depth samples in traveltime table	",
" dzo=dz                 depth interval in traveltime table		",
" fac=0.01               factor to determine the radius of extrap. 	",
" 									",
" Parameters defining the source positions			        ",
" fsx=fx                 x-coordinate of first source                   ",
" nsx=1                  number of sources                              ",
" dsx=2*dxo              x-coordinate increment of sources              ",
" aperx=0.5*nx*dx        ray tracing aperature in x-direction           ",
"									",
" Files for general anisotropic parameters confined to a vertical plane:",
" a1111file		name of file containing a1111(x,z)		",
" a1133file          	name of file containing a1133(x,z)		",
" VS0file	       	name of file containing VS0(x,z)		",
" a1113file          	name of file containing a1113(x,z)		",
" a3313file          	name of file containing a3313(x,z)		",
"									",
" For transversely isotropic media Thomsen's parameters could be used:	",
" deltafile		name of file containing delta(x,z)		",
" epsilonfile		name of file containing epsilon(x,z)		",
"									",
" if anisotropy parameters are not given the program considers		", 
" isotropic media.							", 
"									",
NULL};

/* Credits:
 *      Debashish Sarkar                  
 *		 
 *		
 *   Technical Reference:
 *
 *	Cerveny, V., 1972, Seismic rays and ray intensities 
 *	in inhomogeneous anisotropic media: 
 *	Geophys. J. R. Astr. Soc., 29, 1--13.
 *
 * 	Hangya, A., 1986, Gaussian beams in anisotropic elastic media:
 *      Geophys. J. R. Astr. Soc., 85, 473--563.
 *	
 *	Gajewski, D. and Psencik, I., 1987, Computation of high frequency 
 *	seismic wavefields in 3-D  laterally inhomogeneous anisotropic 
 * 	media: Geophys. J. R. Astr. Soc., 91, 383-411.
 *
 */

/**************** end self doc ********************************/

/* initial value for traveltime as an arbitary large number */
# define INITIAL 99999999

/* functions defined and used internally */
/* one step along ray */
typedef struct RayStepStruct {
	float t;		/* time */
	float x,z;		/* x,z coordinates */
	float q1,p1,q2,p2;	/* Cerveny's dynamic ray tracing solution */
	int kmah;		/* KMAH index */
	float c,s;		/* cos(angle) and sin(angle) */
	float v,dvdx,dvdz;	/* velocity and its derivatives */
} RayStep;

/* one ray */
typedef struct RayStruct {
	int nrs;		/* number of ray steps */
	RayStep *rs;		/* array[nrs] of ray steps */
	int nc;			/* number of circles */
	int ic;			/* index of circle containing nearest step */
	void *c;		/* array[nc] of circles */
} Ray;

/* functions (raytracer)*/
Ray *makeRay (int SV, float x0, float z0, float a0, int nt, float dt, float ft,
	float **a1111xz, float **a3333xz, float **a1133xz, float **a1313xz, 
	float **a1113xz, float **a3313xz,
	int nx, float dx, float fx, int nz, float dz, float fz, float amax, 
	float amin,FILE *fp);
void freeRay (Ray *ray);

/* functions (velocity interpolation)*/
void* vel2Alloc (int nx, float dx, float fx,
	int nz, float dz, float fz, float **v);
void vel2Free (void *vel2);
void vel2Interp (void *vel2, float x, float z,
	float *v, float *vx, float *vz, float *vxx, float *vxz, float *vzz);

int main(int argc, char **argv)
{
	int nx,nz,nt,ix,iz,ia1111=1,ia1133=1,iVS0=1,ia1113=1,ia3313=1,it,
		nxf,nxe,nzf,nze,ixo,izo,i,NRS,SV,
		idelta=1,iepsilon=1,cntt,isx,ia,nxo,nzo,nsx,na;
	float sx,dx,sz,dz,dt,a,fxo,dxo,fzo,dzo,fx,fz,ft,v0,vd,vd1,px,pz,exo,r2,
		r1,r1x,t1x,t2x,t2xz,xoc,zoc,xc,zc,terr,tc,norm2,curv,x1,c1,s1,
		z1,Det,fac,r2x,v1,t1,t2,a1,aperx,*dvdx1,*dvdz1,
		ezo,xo,zo,sd,fsx,dsx,amin,amax,da,fa,**VP0,**a1111,**a3333,
		**VS0,**a1133,**a1313,**a1113,**a3313,**delta,**epsilon,
		*tt,*xx,*zz,*t,*s,*vv,*cc,*ss,**X,**Y,**N;
	char *VS0file="",*a1111file="",*a1133file="",
		*a1113file="",*a3313file="",*deltafile="",*epsilonfile="";
	FILE *a1111fp,*a1133fp,*VP0fp,*VS0fp,
		*a1113fp,*a3313fp,*deltafp,*epsilonfp,*tfp;
	Ray *ray;

	/* hook up getpar */
	//initargs(argc,argv);
	//requestdoc(0);

	/* get required parameters */
	char VP0file[250]={"vel100100.dat"};
        char tfile[250]={"time.dat"};
        char rayfile[250]={"Raypath.txt"};

	nt = 10000;
	nx = 100;
	nz = 100;

	/* get optional parameters */
	/* get parameters for the parameter field*/
	SV=0;
	fx = 0.0;
	fz = 0.0;
	dx = 5.0;
	dz = 5.0;
	dt = 0.0008;
	fac = 0.01;

	/* get parameters for the take-off angle of the ray at the source*/
	amin = 0;
	amax = 180;
	fa=-180;
	na=361;
	da=1;

	/* get parameters defining the output traveltime table*/
	fxo = fx;
	nxo = nx;
	dxo = dx;
	fzo = fz;
	nzo = nz;
	dzo = dz;

        /* get parameters defining source positions */
        fsx = nx*dx/2.0;
        nsx = 1;
        dsx = 2*dxo;
        aperx = 0.5*nx*dx;

	/* output traveltime file */
	tfp = fopen(tfile,"w");

        FILE *fpraypath;
        fpraypath=fopen(rayfile,"w");
	
	/* ensure takeoff point is within model */
        if (fsx<fx || fsx>(fx+dx*(nx-1)) || (fsx+dsx*(nsx-1))<fx 
			|| (fsx+dsx*(nsx-1))>(fx+dx*(nx-1))) 
		err("Shot point locations exceed model dimensions");	
		
	/* allocate workspace */
	VP0 = alloc2float(nz,nx);
	VS0 = alloc2float(nz,nx);
	a1111 = alloc2float(nz,nx);
	a3333 = alloc2float(nz,nx);
	a1133 = alloc2float(nz,nx);
	a1313 = alloc2float(nz,nx);
	a1113 = alloc2float(nz,nx);
	a3313 = alloc2float(nz,nx);
	delta = alloc2float(nz,nx);
	epsilon = alloc2float(nz,nx);
	tt = alloc1float(nt);
	xx = alloc1float(nt);
	zz = alloc1float(nt);
	t  = alloc1float(nxo*nzo);
	s  = alloc1float(nxo*nzo);
	cc = alloc1float(nt);
	ss = alloc1float(nt);
	vv = alloc1float(nt);
	dvdx1 = alloc1float(nt);
	dvdz1 = alloc1float(nt);
	N  = alloc2float(2,2);
	X  = alloc2float(2,2);
	Y  = alloc2float(2,2);

	/* anisotropic parameters */
	ia1111=0;
	ia1133=0;
	iVS0=0;
	ia1113=0;
	ia3313=0;
	idelta=0;      //whether read file
	iepsilon=0;    //whether read file


        //checkpars();

	/* read required velocities*/
	if ((VP0fp=fopen(VP0file,"r"))==NULL)
		err("error opening VP0file=%s!\n",VP0file);
	if (fread(VP0[0],sizeof(float),nz*nx,VP0fp)!=nz*nx)
		err("error reading VP0file=%s!\n",VP0file);

	/* compute a3333 */
	for (ix=0; ix<nx; ++ix)
		for (iz=0; iz<nz; ++iz)
			a3333[ix][iz] = VP0[ix][iz]*VP0[ix][iz];
	

	/* read velocities*/
	if(ia1111){
		if ((a1111fp=fopen(a1111file,"r"))==NULL)
			err("error opening a1111file=%s!\n",a1111file);
		if (fread(a1111[0],sizeof(float),nz*nx,a1111fp)!=nz*nx)
			err("error reading a1111file=%s!\n",a1111file);
	}else{
		for (ix=0; ix<nx; ++ix)
			for (iz=0; iz<nz; ++iz)
				a1111[ix][iz] = a3333[ix][iz];
	}		
	
	if(ia1133){
		if ((a1133fp=fopen(a1133file,"r"))==NULL)
			err("error opening a1133file=%s!\n",a1133file);
		if (fread(a1133[0],sizeof(float),nz*nx,a1133fp)!=nz*nx)
			err("error reading a1133file=%s!\n",a1133file);
	}else{
		for (ix=0; ix<nx; ++ix)
			for (iz=0; iz<nz; ++iz)
				a1133[ix][iz] = .4*a3333[ix][iz];
	}

	if(iVS0){
		if ((VS0fp=fopen(VS0file,"r"))==NULL)
			err("error opening VS0file=%s!\n",VS0file);
		if (fread(VS0[0],sizeof(float),nz*nx,VS0fp)!=nz*nx)
			err("error reading VS0file=%s!\n",VS0file);
	}else{
		for (ix=0; ix<nx; ++ix)
			for (iz=0; iz<nz; ++iz)
				VS0[ix][iz] = .5477*VP0[ix][iz];
	}

	/* compute a1313 */
	for (ix=0; ix<nx; ++ix)
		for (iz=0; iz<nz; ++iz)
			a1313[ix][iz] = VS0[ix][iz]*VS0[ix][iz];

	if(ia1113){
		if ((a1113fp=fopen(a1113file,"r"))==NULL)
			err("error opening a1113file=%s!\n",a1113file);
		if (fread(a1113[0],sizeof(float),nz*nx,a1113fp)!=nz*nx)
			err("error reading a1113file=%s!\n",a1113file);
	}else{
		for (ix=0; ix<nx; ++ix)
			for (iz=0; iz<nz; ++iz)
				a1113[ix][iz] = 0;
	}

	if(ia3313){
		if ((a3313fp=fopen(a3313file,"r"))==NULL)
			err("error opening a3313file=%s!\n",a3313file);
		if (fread(a3313[0],sizeof(float),nz*nx,a3313fp)!=nz*nx)
			err("error reading a3313file=%s!\n",a3313file);
	}else{
		for (ix=0; ix<nx; ++ix)
			for (iz=0; iz<nz; ++iz)
				a3313[ix][iz] = 0;
	}

	/* if specified read Thomsen's parameters*/
	//if(idelta==1 || iepsilon==1) {
		if(idelta){
			if ((deltafp=fopen("delta","r"))==NULL)
				err("error opening deltafile=%s!\n",deltafile);
			if (fread(delta[0],sizeof(float),nz*nx,deltafp)!=nz*nx)
				err("error reading deltafile=%s!\n",deltafile);
		}else{
			for (ix=0; ix<nx; ++ix)
				for (iz=0; iz<nz; ++iz)
					delta[ix][iz] = 0.4;
		}
		if(iepsilon) {
			if ((epsilonfp=fopen("epsilon","r"))==NULL)
			err("error opening epsilonfile=%s!\n",epsilonfile);
			if (fread(epsilon[0],sizeof(float),nz*nx,epsilonfp)!=nz*nx)
			err("error reading epsilonfile=%s!\n",epsilonfile);
		}else{
			for (ix=0; ix<nx; ++ix)
				for (iz=0; iz<nz; ++iz)
					epsilon[ix][iz] = 0.4;
		}

		for (ix=0; ix<nx; ++ix)
			for (iz=0; iz<nz; ++iz){
			a1111[ix][iz] = a3333[ix][iz]+2*epsilon[ix][iz]
					*a3333[ix][iz];
			a1133[ix][iz] = sqrt(2*delta[ix][iz]*a3333[ix][iz]*
				(a3333[ix][iz]-a1313[ix][iz])+(a3333[ix][iz]-
				a1313[ix][iz])*(a3333[ix][iz]-a1313[ix][iz]))
					-a1313[ix][iz];
			}
	//}

	/* iterate over sources	*/
        for(isx=0;isx<nsx;isx++) {
	   sx=fsx+isx*dsx;

	    /* initialize time and raypath length*/
	    for (ixo=0;ixo<nxo;ixo++)
		for(izo=0;izo<nzo;izo++) {
		     i = izo+ixo*nzo;
		     t[i] = INITIAL;
		     s[i] = INITIAL;
		}

		/* iterate over takeoff angle */
		for(ia=0;ia<na;ia++) {
				a=fa+ia*da;
                                if(ia%10==0)	
				    printf("sx = %f angle = %f\n",sx,a);
				/* Assume 0s as initial time */
        			ft=0; sz=nz*dz/2.0;
				ray = makeRay(SV,sx,sz,a,nt,dt,ft,a1111,a3333,
					a1133,a1313,a1113,a3313,nx,
					dx,fx,nz,dz,fz,amax,amin,fpraypath);
					/* store the ray */
					for(cntt=0;cntt<nt;cntt++){
						tt[cntt]=ray->rs[cntt].t;
						xx[cntt]=ray->rs[cntt].x;
						zz[cntt]=ray->rs[cntt].z;
						cc[cntt]=ray->rs[cntt].c;
						ss[cntt]=ray->rs[cntt].s;
						vv[cntt]=ray->rs[cntt].v;
						dvdx1[cntt]=ray->rs[cntt].dvdx;
						dvdz1[cntt]=ray->rs[cntt].dvdz;
					      }
		NRS=ray->nrs;

		/* free ray */
		freeRay(ray);

		/* Extrapolate traveltimes around the central ray 
		(Similar to Zhenyue Liu's code, rayt2d)*/
		v0 = vv[0];
		r2 = 0.0;
		xc = sx;
		zc = sz;
		sd = 0;
		vd1= v0;
		
		/* trace the auxillary ray */		
		a1=a+0.1;
		ray=makeRay(SV,sx,sz,a1,NRS,dt,ft,a1111,a3333,a1133,a1313,
			a1113,a3313,nx,dx,fx,nz,dz,fz,amax,amin,fpraypath);

		for (it=1;it<NRS;++it) {
			xo = xx[it];
			zo = zz[it];
			vd = vv[it];
			sd = sd+0.5*dt*(vd+vd1);
			vd1= vd;
			curv = fabs(-(dvdx1[it]*cc[it]-dvdz1[it]*ss[it]))/vv[it];
		
		/* Check if the point is in the previous circle */
		if (r2 > (xc-xo)*(xc-xo)+(zc-zo)*(zc-zo) && it != NRS-1) continue;

		exo = fxo+(nxo-1)*dxo;
		ezo = fzo+(nzo-1)*dzo;

		/* Check if the point is out of output range */
		if (xo<fxo || xo>exo || zo<fzo || zo>ezo) continue;
	
		xc = xo; zc = zo;
		px = ss[it]/vd1;
		pz = cc[it]/vd1;
		
		/* set parameters of the first auxilliary ray */
		x1 = ray->rs[it].x;
		z1 = ray->rs[it].z;
		c1 = ray->rs[it].c;
		s1 = ray->rs[it].s;
		v1 = ray->rs[it].v;

		/* Compute the second derivatives of traveltimes */	
		/* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

		/* compute Xij and Yij */
		X[0][0] = (x1-xc)/0.0017;
		X[1][0] = (z1-zc)/0.0017;
		Y[0][0] = ((s1/v1)-px)/0.0017;
		Y[1][0] = ((c1/v1)-pz)/0.0017;
		X[0][1] = (xx[it+1]-xx[it])/dt;
		X[1][1] = (zz[it+1]-zz[it])/dt;
		Y[0][1] = ((ss[it+1]/vv[it+1])-(ss[it]/vv[it]))/dt;
		Y[1][1] = ((cc[it+1]/vv[it+1])-(cc[it]/vv[it]))/dt;
	
		/* compute Nij */
		Det = X[0][0]*X[1][1]-X[0][1]*X[1][0];
		N[0][0] = (Y[0][0]*X[1][1]-Y[0][1]*X[1][0])/Det;
		N[0][1] = (-Y[0][0]*X[0][1]+Y[0][1]*X[0][0])/Det;
		N[1][0] = (Y[1][0]*X[1][1]-Y[1][1]*X[1][0])/Det;
		N[1][1] = (-Y[1][0]*X[0][1]+Y[1][1]*X[0][0])/Det;	
	 	/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/	
		
		/* determine the radius of extrapolation */
		tc = it*dt;
		terr = tc*fac;
		norm2 = sqrt(N[0][0]*N[0][0]+N[1][1]*N[1][1]+2*N[0][1]*N[0][1]);
		r2 = terr/norm2;
		r1 = sqrt(r2);
		
		/* Make sure that the coordinates are regular */	
		if (r1*curv > 0.1) r1 = 0.1/curv;
		/* Ensure that the radius is not large */
		if (r1 > 0.1*sd) r1 = 0.1*sd;
		
		r2 = r1*r1;
		
		nxf = ((xc-r1-fxo)/dxo)+0.9;
		if (nxf < 0) nxf=0;
		nxe = ((xc+r1-fxo)/dxo)+0.1;
		if (nxe >= nxo) nxe = nxo-1;

/***************** output raypath*************************/
                     

				for(ixo=nxf; ixo<=nxe; ++ixo) {
					xoc = fxo-xc+ixo*dxo;
					r2x = r2-xoc*xoc;
					if(r2x<0) continue;
					
					r1x = sqrt(r2x);
					t1x = tc+xoc*px;
					t2x = xoc*xoc*N[0][0];
					t2xz= 2*xoc*N[0][1];
					nzf = (zc-r1x-fzo)/dzo+0.9;
					if (nzf < 0) nzf = 0;
					nze = (zc+r1x-fzo)/dzo+0.1;
					if (nze>=nzo) nze = nzo-1;
					
					for(izo=nzf;izo<=nze;++izo) {
						if (sd<s[nzo*ixo+izo]) {
						zoc = fzo-zc+izo*dzo;
						t1 = t1x+zoc*pz;
						t2 = t2x+zoc*(zoc*N[1][1]+t2xz);
						s[nzo*ixo+izo] = sd;
						t[nzo*ixo+izo] = t1+0.5*t2;
							}
						}
				}
		}
		
		/* free ray */
		freeRay(ray);
				
		}	
                for(i=0;i<nxo*nzo;i++)
                      if((t[i]==999999)||(t[i]>99999))
                             t[i]=0.0;
		fwrite(t,sizeof(float),nxo*nzo,tfp);
	}
		 		

	/* free workspace */
	free2float(delta);
	free2float(epsilon);

	/* free workspace */
	free2float(VP0);
	free2float(VS0);
	free2float(a1111);
	free2float(a3333);
	free2float(a1133);
	free2float(a1313);
	free2float(a1113);
	free2float(a3313);

        fclose(VP0fp);
        fclose(tfp);
        fclose(fpraypath);
	
	return 0;
}

/*****************************************************************
* Functions for Ray tracing					 *	
*****************************************************************/

Ray *makeRay (int SV, float x0, float z0, float a0, int nt, float dt, float ft,
	float **a1111xz, float **a3333xz, float **a1133xz, float **a1313xz, 
	float **a1113xz, float **a3313xz, 
	int nx, float dx, float fx, int nz, float dz, float fz, 
        float amax, float amin,FILE *fpraypath)
/*****************************************************************************
Trace a ray for uniformly sampled v(x,z).
******************************************************************************
Input:
x0		x coordinate of takeoff point
z0		z coordinate of takeoff point
a0		takeoff angle (radians)
nt		number of time samples
dt		time sampling interval
ft		first time sample
nx		number of x samples
dx		x sampling interval
fx		first x sample
nz		number of z samples
dz		z sampling interval
fz		first z sample
amax            maximum emergence angle
amin            minimum emergence angle
a1111		array[nx][nz] of uniformly sampled density normalized elastic coef.
a3333		array[nx][nz] of uniformly sampled density normalized elastic coef.
a1133           array[nx][nz] of uniformly sampled density normalized elastic coef.
a1313           array[nx][nz] of uniformly sampled density normalized elastic coef.
a1113           array[nx][nz] of uniformly sampled density normalized elastic coef.
a3313           array[nx][nz] of uniformly sampled density normalized elastic coef.
******************************************************************************
Returned:	pointer to ray parameters sampled at discrete ray steps
******************************************************************************
Notes:
The ray ends when it runs out of time (after nt steps) or with the first 
step that is out of the (x,z) bounds of the velocity function v(x,z).
*****************************************************************************
Technical Reference:

Cerveny, V., 1972, Seismic rays and ray intensities 
	in inhomogeneous anisotropic media: 
	Geophys. J. R. Astr. Soc., 29, 1--13.

*****************************************************************************
 Credits: CWP: Tariq Alkhalifah
*****************************************************************************/
{
	int it,kmah;
	float t,x,z,c,s,p1,q1,p2,q2,px,pz,px2,pz2,pxz,
		lx,lz,cc,ss,
		a1111,da1111dx,da1111dz,dda1111dxdx,dda1111dxdz,dda1111dzdz,
		a3333,da3333dx,da3333dz,dda3333dxdx,dda3333dxdz,dda3333dzdz,
		a1133,da1133dx,da1133dz,dda1133dxdx,dda1133dxdz,dda1133dzdz,
		a1313,da1313dx,da1313dz,dda1313dxdx,dda1313dxdz,dda1313dzdz,
		a1113,da1113dx,da1113dz,dda1113dxdx,dda1113dxdz,dda1113dzdz,
		a3313,da3313dx,da3313dz,dda3313dxdx,dda3313dxdz,dda3313dzdz,
		da1111dn,dda1111dndn,da3333dn,dda3333dndn,da1133dn,dda1133dndn,
		da1313dn,dda1313dndn,da1113dn,dda1113dndn,da3313dn,dda3313dndn,
		gamm11,gamm13,gamm33,vp2,vp,ovp,sqr;
	Ray *ray;
	RayStep *rs;
	void *a11112;
	void *a33332;
	void *a11332;
	void *a13132;
	void *a11132;
	void *a33132;

	/*Convert degrees to radians   度转换为弧度 */
	a0=a0*3.1415926/180; amax=amax*3.1415926/180; amin=amin*3.1415926/180;
	/* allocate and initialize velocities interpolator 
               分配和初始化速度插值       */
	a11112 = vel2Alloc(nx,dx,fx,nz,dz,fz,a1111xz);
	a33332 = vel2Alloc(nx,dx,fx,nz,dz,fz,a3333xz);
	a11332 = vel2Alloc(nx,dx,fx,nz,dz,fz,a1133xz);
	a13132 = vel2Alloc(nx,dx,fx,nz,dz,fz,a1313xz);
	a11132 = vel2Alloc(nx,dx,fx,nz,dz,fz,a1113xz);
	a33132 = vel2Alloc(nx,dx,fx,nz,dz,fz,a3313xz);

	/* last x and z in velocity model */
	lx = fx+(nx-1)*dx;
	lz = fz+(nz-1)*dz;

	/* ensure takeoff point is within model */
	if (x0<fx || x0>lx || z0<fz || z0>lz ) return NULL;

	/* allocate space for ray and raysteps */
	ray = (Ray*)alloc1(1,sizeof(Ray));
	rs = (RayStep*)alloc1(nt,sizeof(RayStep));

	/* cosine and sine of takeoff angle */
	c = cos(a0);
	s = sin(a0);
	cc = c*c;
	ss = s*s;

	/* velocities and derivatives at takeoff point */
	vel2Interp(a11112,x0,z0,&a1111,&da1111dx,&da1111dz,&dda1111dxdx,
		&dda1111dxdz,&dda1111dzdz);
	da1111dn    = da1111dx*c-da1111dz*s;
	dda1111dndn = dda1111dxdx*cc-2.0*dda1111dxdz*s*c+dda1111dzdz*ss;

	vel2Interp(a33332,x0,z0,&a3333,&da3333dx,&da3333dz,&dda3333dxdx,
		&dda3333dxdz,&dda3333dzdz);
	da3333dn    = da3333dx*c-da3333dz*s;
	dda3333dndn = dda3333dxdx*cc-2.0*dda3333dxdz*s*c+dda3333dzdz*ss;
	
	vel2Interp(a11332,x0,z0,&a1133,&da1133dx,&da1133dz,&dda1133dxdx,
		&dda1133dxdz,&dda1133dzdz);
	da1133dn    = da1133dx*c-da1133dz*s;
	dda1133dndn = dda1133dxdx*cc-2.0*dda1133dxdz*s*c+dda1133dzdz*ss;

	vel2Interp(a13132,x0,z0,&a1313,&da1313dx,&da1313dz,&dda1313dxdx,
		&dda1313dxdz,&dda1313dzdz);
	da1313dn    = da1313dx*c-da1313dz*s;
	dda1313dndn = dda1313dxdx*cc-2.0*dda1313dxdz*s*c+dda1313dzdz*ss;

	vel2Interp(a11132,x0,z0,&a1113,&da1113dx,&da1113dz,&dda1113dxdx,
		&dda1113dxdz,&dda1113dzdz);
	da1113dn    = da1113dx*c-da1113dz*s;
	dda1113dndn = dda1113dxdx*cc-2.0*dda1113dxdz*s*c+dda1113dzdz*ss;

	vel2Interp(a33132,x0,z0,&a3313,&da3313dx,&da3313dz,&dda3313dxdx,
		&dda3313dxdz,&dda3313dzdz);
	da3313dn    = da3313dx*c-da3313dz*s;
	dda3313dndn = dda3313dxdx*cc-2.0*dda3313dxdz*s*c+dda3313dzdz*ss;

	/*computing the phase velocity for a0 angle */
	gamm11 = a1111*ss+ a1313*cc +2*a1113*s*c;
	gamm33 = a3333*cc + a1313*ss+2*a3313*s*c;
	gamm13 = (a1133+a1313)*s*c+ a1113*ss+ a3313*cc;
	sqr    = sqrt((gamm11+gamm33)*(gamm11+gamm33)-
			4*(gamm11*gamm33-gamm13*gamm13));
	if(SV)
	vp2    = gamm11+gamm33-sqr;
	else
	vp2    = gamm11+gamm33+sqr;
	vp     = sqrt(vp2*.5);
	ovp    = 1/vp;
	px     = s*ovp;
	pz     = c*ovp;

	/* first ray step */
	rs[0].t = t = ft;
	rs[0].x = x = x0;
	rs[0].z = z = z0;
	rs[0].q1 = q1 = 1.0;
	rs[0].p1 = p1 = 0.0;
	rs[0].q2 = q2 = 0.0;
	rs[0].p2 = p2 = 1.0;
	rs[0].kmah = kmah = 0;
	rs[0].c = c;
	rs[0].s = s;
	rs[0].v = vp;
	rs[0].dvdx = .5*da3333dx*vp/a3333;
	rs[0].dvdz = .5*da3333dz*vp/a3333;

	/* loop over time steps */
	for (it=1; it<nt; ++it) {

		/* variables used for Runge-Kutta integration */
		float h=dt,hhalf=dt/2.0,hsixth=dt/6.0,
			q2old,xt,zt,p1t,q1t,p2t,q2t,
			dx,dz,dp1,dq1,dp2,dq2,
			dxt,dzt,dp1t,dq1t,dp2t,dq2t,
			dxm,dzm,dp1m,dq1m,dp2m,dq2m,
			gamma11,gamma33,gamma13,g11,g13,g33,den,
			sxfactor,szfactor,snfact,dpx,dpz,pxt,pzt,dpxt,
			dpzt,dpxm,dpzm,dxdn,dzdn,snfactor,
			dxx,dzz,dcdp1,dcdp3,dcdp13,ddcdnn,ddcdqq,
			ddcdpn,dgamma11dpx,dgamma11dpz,dgamma33dpx,
			dgamma33dpz,dgamma13dpx,dgamma13dpz,dg11dpx,
			dg11dpz,dg33dpx,dg33dpz,dg13dpx,dg13dpz,ddxdpx,
			ddzdpz,ddxdpz,dgamma11dn,dgamma33dn,dgamma13dn,
			dg11dn,dg33dn,dg13dn,dsnfactordn,ddxdn,ddzdn;

		/* if ray is out of bounds, break */
		if (x<fx || x>lx || z<fz || z>lz || c>(cos(amin)+0.01) || c<(cos(amax))-0.01) break;

		/* remember old q2 */
		q2old = q2;
		
	/* step 1 of 4th-order Runge-Kutta */
		px2   = px*px;
		pz2   = pz*pz;
		pxz   = px*pz;

		/*anisotropy parameters*/
		gamma11 = a1111*px2+ a1313*pz2 +2*a1113*pxz;
		gamma33 = a3333*pz2 + a1313*px2+2*a3313*pxz;
		gamma13 = (a1133+a1313)*pxz+ a1113*px2+ a3313*pz2;
		den     = 1/(gamma11+gamma33-2);
		g11     = (gamma33-1)*den;
		g33     = (gamma11-1)*den;
		g13     = -gamma13*den;
		sxfactor = da1111dx*px2*g11+da3333dx*pz2*g33+
			2*(da1133dx+da1313dx)*pxz*g13+da1313dx*(px2*g33+pz2*g11)+
			2*da3313dx*(pz2*g13+pxz*g33)+2*da1113dx*(pxz*g11+px2*g13);
		szfactor = da1111dz*px2*g11+da3333dz*pz2*g33+
			2*(da1133dz+da1313dz)*pxz*g13+da1313dz*(px2*g33+pz2*g11)+
			2*da3313dz*(pz2*g13+pxz*g33)+2*da1113dz*(pxz*g11+px2*g13);
		snfact = sxfactor*c-szfactor*s;
		
		/*computing ray velocities and derivatives*/
		dx =  (a1111*px*g11+(a1133+a1313)*pz*g13+a3313*pz*g33
			+a1113*(pz*g11+2*px*g13)+a1313*g33*px);
		dz =  (a3333*pz*g33+(a1133+a1313)*px*g13+a1113*px*g11
			+a3313*(px*g33+2*pz*g13)+a1313*g11*pz);

		dgamma11dpx = 2*a1111*px+2*a1113*pz;
		dgamma11dpz = 2*a1313*pz+2*a1113*px;
		dgamma33dpx = 2*a1313*px+2*a3313*pz;
		dgamma33dpz = 2*a3333*pz+2*a3313*px;
		dgamma13dpx= (a1133+a1313)*pz+2*a1113*px;
		dgamma13dpz= (a1133+a1313)*px+2*a3313*pz;
		dgamma11dn = da1111dn*px2+ da1313dn*pz2 +2*da1113dn*pxz;
		dgamma33dn = da3333dn*pz2 + da1313dn*px2+2*da3313dn*pxz;
		dgamma13dn = (da1133dn+da1313dn)*pxz+ da1113dn*px2+ da3313dn*pz2;
		dg11dpx = -(gamma33-1)*(dgamma11dpx+dgamma33dpx-4*dx)*den*den+
			(dgamma33dpx-2*dx)*den;
		dg11dpz = -(gamma33-1)*(dgamma11dpz+dgamma33dpz-4*dz)*den*den+
			(dgamma33dpz-2*dz)*den;
		dg33dpx = -(gamma11-1)*(dgamma11dpx+dgamma33dpx-4*dx)*den*den+
			(dgamma11dpx-2*dx)*den;
		dg33dpz = -(gamma11-1)*(dgamma11dpz+dgamma33dpz-4*dz)*den*den+
			(dgamma11dpz-2*dz)*den;
		dg13dpx = gamma13*(dgamma11dpx+dgamma33dpx-4*dx)*den*den-
			dgamma13dpx*den;
		dg13dpz = gamma13*(dgamma11dpz+dgamma33dpz-4*dz)*den*den-
			dgamma13dpz*den;
		dg11dn = -(gamma33-1)*(dgamma11dn+dgamma33dn-2*snfact)*den*den+
			(dgamma33dn-snfact)*den;
		dg33dn = -(gamma11-1)*(dgamma11dn+dgamma33dn-2*snfact)*den*den+
			(dgamma11dn-snfact)*den;
		dg13dn = gamma13*(dgamma11dn+dgamma33dn-2*snfact)*den*den-
			dgamma13dn*den;
		ddxdpx=   a1111*px*dg11dpx+(a1133+a1313)*pz*dg13dpx+
			a3313*pz*dg33dpx+a1113*(pz*dg11dpx+2*px*dg13dpx)
			+a1313*dg33dpx*px;
		ddzdpz= a3333*pz*dg33dpz+(a1133+a1313)*px*dg13dpz+
			a1113*px*dg11dpz+a3313*(px*dg33dpz+2*pz*dg13dpz)+
			a1313*dg11dpz*pz;
		ddxdpz= a1111*px*dg11dpz+(a1133+a1313)*pz*dg13dpz+
			a3313*pz*dg33dpz+a1113*(pz*dg11dpz+2*px*dg13dpz)+
			a1313*dg33dpz*px;
		dsnfactordn = da1111dn*px2*dg11dn+da3333dn*pz2*dg33dn+
			2*(da1133dn+da1313dn)*pxz*dg13dn+da1313dn*(px2*dg33dn+pz2*dg11dn)+
			2*da3313dn*(pz2*dg13dn+pxz*dg33dn)+2*da1113dn*(pxz*dg11dn+px2*dg13dn);
		ddxdn =  (a1111*px*dg11dn+(a1133+a1313)*pz*dg13dn+a3313*pz*dg33dn
			+a1113*(pz*dg11dn+2*px*dg13dn)+a1313*dg33dn*px);
		ddzdn =  (a3333*pz*dg33dn+(a1133+a1313)*px*dg13dn+a1113*px*dg11dn
			+a3313*(px*dg33dn+2*pz*dg13dn)+a1313*dg11dn*pz);
		

		/*evaluating change in slowness and amplitude along ray*/
		dpx = -0.5*sxfactor;
		dpz = -0.5*szfactor;

		dcdp1  = a1111*g11+a1313*g33+2*a1113*g13+ddxdpx-dx*dx;
		dcdp3  = a3333*g33+a1313*g11+2*a3313*g13+ddzdpz-dz*dz;
		dcdp13 = a1133*g13+a1313*g13+a1113*g11+a3313*g33+ddxdpz-dx*dz;
		ddcdqq = dcdp1*cc-2.0*dcdp13*s*c+dcdp3*ss;
		dxdn   =  (da1111dn*px*g11+(da1133dn+da1313dn)*pz*g13+da3313dn*pz*g33
			+da1113dn*(pz*g11+2*px*g13)+da1313dn*g33*px);
		dzdn   =  (da3333dn*pz*g33+(da1133dn+da1313dn)*px*g13+da1113dn*px*g11
			+da3313dn*(px*g33+2*pz*g13)+da1313dn*g11*pz);
		ddcdpn = dxdn*c-dzdn*s-.5*dx*sxfactor*cc+
			.5*(dx*szfactor+dz*sxfactor)*s*c-.5*dz*szfactor*ss
			+ddxdn*c-ddzdn*s;
		snfactor = dda1111dndn*px2*g11+dda3333dndn*pz2*g33+
			2*(dda1133dndn+dda1313dndn)*pxz*g13+
			dda1313dndn*(px2*g33+pz2*g11)+
			2*dda3313dndn*(pz2*g13+pxz*g33)+
			2*dda1113dndn*(pxz*g11+px2*g13);
		ddcdnn = 0.5*snfactor-.25*sxfactor*sxfactor*cc+
			.5*sxfactor*szfactor*s*c-.25*szfactor*szfactor*ss
			+.5*dsnfactordn;

		dp1 = -ddcdnn*q1-ddcdpn*p1;
		dq1 = ddcdqq*p1+ddcdpn*q1;
		dp2 = -ddcdnn*q2-ddcdpn*p2;
		dq2 = ddcdqq*p2+ddcdpn*q2;
		xt = x+hhalf*dx;
		zt = z+hhalf*dz;
		pxt = px+hhalf*dpx;
		pzt = pz+hhalf*dpz;
		p1t = p1+hhalf*dp1;
		q1t = q1+hhalf*dq1;
		p2t = p2+hhalf*dp2;
		q2t = q2+hhalf*dq2;
		vp  = 1/sqrt(pxt*pxt+pzt*pzt);
		s   = pxt*vp;
		c   = pzt*vp;
		ss  = s*s;
		cc  = c*c;
		
		vel2Interp(a11112,xt,zt,&a1111,&da1111dx,&da1111dz,&dda1111dxdx,
			&dda1111dxdz,&dda1111dzdz);
		da1111dn    = da1111dx*c-da1111dz*s;
		dda1111dndn = dda1111dxdx*cc-2.0*dda1111dxdz*s*c+dda1111dzdz*ss;

		vel2Interp(a33332,xt,zt,&a3333,&da3333dx,&da3333dz,&dda3333dxdx,
		&dda3333dxdz,&dda3333dzdz);
		da3333dn    = da3333dx*c-da3333dz*s;
		dda3333dndn = dda3333dxdx*cc-2.0*dda3333dxdz*s*c+dda3333dzdz*ss;
	
		vel2Interp(a11332,xt,zt,&a1133,&da1133dx,&da1133dz,&dda1133dxdx,
			&dda1133dxdz,&dda1133dzdz);
		da1133dn    = da1133dx*c-da1133dz*s;
		dda1133dndn = dda1133dxdx*cc-2.0*dda1133dxdz*s*c+dda1133dzdz*ss;

		vel2Interp(a13132,xt,zt,&a1313,&da1313dx,&da1313dz,&dda1313dxdx,
			&dda1313dxdz,&dda1313dzdz);
		da1313dn    = da1313dx*c-da1313dz*s;
		dda1313dndn = dda1313dxdx*cc-2.0*dda1313dxdz*s*c+dda1313dzdz*ss;

		vel2Interp(a11132,xt,zt,&a1113,&da1113dx,&da1113dz,&dda1113dxdx,
			&dda1113dxdz,&dda1113dzdz);
		da1113dn    = da1113dx*c-da1113dz*s;
		dda1113dndn = dda1113dxdx*cc-2.0*dda1113dxdz*s*c+dda1113dzdz*ss;

		vel2Interp(a33132,xt,zt,&a3313,&da3313dx,&da3313dz,&dda3313dxdx,
			&dda3313dxdz,&dda3313dzdz);
		da3313dn    = da3313dx*c-da3313dz*s;
		dda3313dndn = dda3313dxdx*cc-2.0*dda3313dxdz*s*c+dda3313dzdz*ss;
		
	/* step 2 of 4th-order Runge-Kutta */
		px2   = pxt*pxt;
		pz2   = pzt*pzt;
		pxz   = pxt*pzt;

		/*anisotropy parameters*/
		gamma11 = a1111*px2+ a1313*pz2 +2*a1113*pxz;
		gamma33 = a3333*pz2 + a1313*px2+2*a3313*pxz;
		gamma13 = (a1133+a1313)*pxz+ a1113*px2+ a3313*pz2;
		den     = 1/(gamma11+gamma33-2);
		g11     = (gamma33-1)*den;
		g33     = (gamma11-1)*den;
		g13     = -gamma13*den;
		sxfactor = da1111dx*px2*g11+da3333dx*pz2*g33+
			2*(da1133dx+da1313dx)*pxz*g13+da1313dx*(px2*g33+pz2*g11)+
			2*da3313dx*(pz2*g13+pxz*g33)+2*da1113dx*(pxz*g11+px2*g13);
		szfactor = da1111dz*px2*g11+da3333dz*pz2*g33+
			2*(da1133dz+da1313dz)*pxz*g13+da1313dz*(px2*g33+pz2*g11)+
			2*da3313dz*(pz2*g13+pxz*g33)+2*da1113dz*(pxz*g11+px2*g13);
		snfact = sxfactor*c-szfactor*s;
		
		/*computing ray velocities and derivatives*/
		dxt =  (a1111*pxt*g11+(a1133+a1313)*pzt*g13+a3313*pzt*g33
			+a1113*(pzt*g11+2*pxt*g13)+a1313*g33*pxt);
		dzt =  (a3333*pzt*g33+(a1133+a1313)*pxt*g13+a1113*pxt*g11
			+a3313*(pxt*g33+2*pzt*g13)+a1313*g11*pzt);
		dpxt = -0.5*sxfactor;
		dpzt = -0.5*szfactor;

		dgamma11dpx = 2*a1111*pxt+2*a1113*pzt;
		dgamma11dpz = 2*a1313*pzt+2*a1113*pxt;
		dgamma33dpx = 2*a1313*pxt+2*a3313*pzt;
		dgamma33dpz = 2*a3333*pzt+2*a3313*pxt;
		dgamma13dpx= (a1133+a1313)*pzt+2*a1113*pxt;
		dgamma13dpz= (a1133+a1313)*pxt+2*a3313*pzt;
		dgamma11dn = da1111dn*px2+ da1313dn*pz2 +2*da1113dn*pxz;
		dgamma33dn = da3333dn*pz2 + da1313dn*px2+2*da3313dn*pxz;
		dgamma13dn = (da1133dn+da1313dn)*pxz+ da1113dn*px2+ da3313dn*pz2;
		dg11dpx = -(gamma33-1)*(dgamma11dpx+dgamma33dpx-4*dxt)*den*den+
			(dgamma33dpx-2*dxt)*den;
		dg11dpz = -(gamma33-1)*(dgamma11dpz+dgamma33dpz-4*dzt)*den*den+
			(dgamma33dpz-2*dzt)*den;
		dg33dpx = -(gamma11-1)*(dgamma11dpx+dgamma33dpx-4*dxt)*den*den+
			(dgamma11dpx-2*dxt)*den;
		dg33dpz = -(gamma11-1)*(dgamma11dpz+dgamma33dpz-4*dzt)*den*den+
			(dgamma11dpz-2*dzt)*den;
		dg13dpx = gamma13*(dgamma11dpx+dgamma33dpx-4*dxt)*den*den-
			dgamma13dpx*den;
		dg13dpz = gamma13*(dgamma11dpz+dgamma33dpz-4*dzt)*den*den-
			dgamma13dpz*den;
		dg11dn = -(gamma33-1)*(dgamma11dn+dgamma33dn-2*snfact)*den*den+
			(dgamma33dn-snfact)*den;
		dg33dn = -(gamma11-1)*(dgamma11dn+dgamma33dn-2*snfact)*den*den+
			(dgamma11dn-snfact)*den;
		dg13dn = gamma13*(dgamma11dn+dgamma33dn-2*snfact)*den*den-
			dgamma13dn*den;
		ddxdpx=   a1111*pxt*dg11dpx+(a1133+a1313)*pzt*dg13dpx+
			a3313*pzt*dg33dpx+a1113*(pzt*dg11dpx+2*pxt*dg13dpx)
			+a1313*dg33dpx*pxt;
		ddzdpz= a3333*pzt*dg33dpz+(a1133+a1313)*pxt*dg13dpz+
			a1113*pxt*dg11dpz+a3313*(pxt*dg33dpz+2*pzt*dg13dpz)+
			a1313*dg11dpz*pzt;
		ddxdpz= a1111*pxt*dg11dpz+(a1133+a1313)*pzt*dg13dpz+
			a3313*pzt*dg33dpz+a1113*(pzt*dg11dpz+2*pxt*dg13dpz)+
			a1313*dg33dpz*pxt;
		dsnfactordn = da1111dn*px2*dg11dn+da3333dn*pz2*dg33dn+
			2*(da1133dn+da1313dn)*pxz*dg13dn+da1313dn*(px2*dg33dn+pz2*dg11dn)+
			2*da3313dn*(pz2*dg13dn+pxz*dg33dn)+2*da1113dn*(pxz*dg11dn+px2*dg13dn);
		ddxdn =  (a1111*pxt*dg11dn+(a1133+a1313)*pzt*dg13dn+a3313*pzt*dg33dn
			+a1113*(pzt*dg11dn+2*pxt*dg13dn)+a1313*dg33dn*pxt);
		ddzdn =  (a3333*pzt*dg33dn+(a1133+a1313)*pxt*dg13dn+a1113*pxt*dg11dn
			+a3313*(pxt*dg33dn+2*pzt*dg13dn)+a1313*dg11dn*pzt);
		
		dcdp1  = a1111*g11+a1313*g33+2*a1113*g13+ddxdpx-dxt*dxt;
		dcdp3  = a3333*g33+a1313*g11+2*a3313*g13+ddzdpz-dzt*dzt;
		dcdp13 = a1133*g13+a1313*g13+a1113*g11+a3313*g33+ddxdpz-dxt*dzt;
		ddcdqq = dcdp1*cc-2.0*dcdp13*s*c+dcdp3*ss;
		dxdn   =  (da1111dn*pxt*g11+(da1133dn+da1313dn)*pzt*g13+
			da3313dn*pzt*g33+da1113dn*(pzt*g11+2*pxt*g13)+
			da1313dn*g33*pxt);
		dzdn   =  (da3333dn*pzt*g33+(da1133dn+da1313dn)*pxt*g13+
			da1113dn*pxt*g11+da3313dn*(pxt*g33+2*pzt*g13)+
			da1313dn*g11*pzt);
		ddcdpn = dxdn*c-dzdn*s-.5*dxt*sxfactor*cc+
			.5*(dxt*szfactor+dzt*sxfactor)*s*c-.5*dzt*szfactor*ss
			+ddxdn*c-ddzdn*s;
		snfactor = dda1111dndn*px2*g11+dda3333dndn*pz2*g33+
			2*(dda1133dndn+dda1313dndn)*pxz*g13+
			dda1313dndn*(px2*g33+pz2*g11)+
			2*dda3313dndn*(pz2*g13+pxz*g33)+
			2*dda1113dndn*(pxz*g11+px2*g13);
		ddcdnn = 0.5*snfactor-.25*sxfactor*sxfactor*cc+
			.5*sxfactor*szfactor*s*c-.25*szfactor*szfactor*ss
			+.5*dsnfactordn;

		dp1t = -ddcdnn*q1t-ddcdpn*p1t;
		dq1t = ddcdqq*p1t+ddcdpn*q1t;
		dp2t = -ddcdnn*q2t-ddcdpn*p2t;
		dq2t = ddcdqq*p2t+ddcdpn*q2t;
		xt = x+hhalf*dxt;
		zt = z+hhalf*dzt;
		pxt = px+hhalf*dpxt;
		pzt = pz+hhalf*dpzt;
		p1t = p1+hhalf*dp1t;
		q1t = q1+hhalf*dq1t;
		p2t = p2+hhalf*dp2t;
		q2t = q2+hhalf*dq2t;
		vp  = 1/sqrt(pxt*pxt+pzt*pzt);
		s   = pxt*vp;
		c   = pzt*vp;
		ss  = s*s;
		cc  = c*c;
		
		vel2Interp(a11112,xt,zt,&a1111,&da1111dx,&da1111dz,&dda1111dxdx,
			&dda1111dxdz,&dda1111dzdz);
		da1111dn    = da1111dx*c-da1111dz*s;
		dda1111dndn = dda1111dxdx*cc-2.0*dda1111dxdz*s*c+dda1111dzdz*ss;

		vel2Interp(a33332,xt,zt,&a3333,&da3333dx,&da3333dz,&dda3333dxdx,
		&dda3333dxdz,&dda3333dzdz);
		da3333dn    = da3333dx*c-da3333dz*s;
		dda3333dndn = dda3333dxdx*cc-2.0*dda3333dxdz*s*c+dda3333dzdz*ss;
	
		vel2Interp(a11332,xt,zt,&a1133,&da1133dx,&da1133dz,&dda1133dxdx,
			&dda1133dxdz,&dda1133dzdz);
		da1133dn    = da1133dx*c-da1133dz*s;
		dda1133dndn = dda1133dxdx*cc-2.0*dda1133dxdz*s*c+dda1133dzdz*ss;

		vel2Interp(a13132,xt,zt,&a1313,&da1313dx,&da1313dz,&dda1313dxdx,
			&dda1313dxdz,&dda1313dzdz);
		da1313dn    = da1313dx*c-da1313dz*s;
		dda1313dndn = dda1313dxdx*cc-2.0*dda1313dxdz*s*c+dda1313dzdz*ss;

		vel2Interp(a11132,xt,zt,&a1113,&da1113dx,&da1113dz,&dda1113dxdx,
			&dda1113dxdz,&dda1113dzdz);
		da1113dn    = da1113dx*c-da1113dz*s;
		dda1113dndn = dda1113dxdx*cc-2.0*dda1113dxdz*s*c+dda1113dzdz*ss;

		vel2Interp(a33132,xt,zt,&a3313,&da3313dx,&da3313dz,&dda3313dxdx,
			&dda3313dxdz,&dda3313dzdz);
		da3313dn    = da3313dx*c-da3313dz*s;
		dda3313dndn = dda3313dxdx*cc-2.0*dda3313dxdz*s*c+dda3313dzdz*ss;
		
	/* step 3 of 4th-order Runge-Kutta */
		px2   = pxt*pxt;
		pz2   = pzt*pzt;
		pxz   = pxt*pzt;

		/*anisotropy parameters*/
		gamma11 = a1111*px2+ a1313*pz2 +2*a1113*pxz;
		gamma33 = a3333*pz2 + a1313*px2+2*a3313*pxz;
		gamma13 = (a1133+a1313)*pxz+ a1113*px2+ a3313*pz2;
		den     = 1/(gamma11+gamma33-2);
		g11     = (gamma33-1)*den;
		g33     = (gamma11-1)*den;
		g13     = -gamma13*den;
		sxfactor = da1111dx*px2*g11+da3333dx*pz2*g33+
			2*(da1133dx+da1313dx)*pxz*g13+da1313dx*(px2*g33+pz2*g11)+
			2*da3313dx*(pz2*g13+pxz*g33)+2*da1113dx*(pxz*g11+px2*g13);
		szfactor = da1111dz*px2*g11+da3333dz*pz2*g33+
			2*(da1133dz+da1313dz)*pxz*g13+da1313dz*(px2*g33+pz2*g11)+
			2*da3313dz*(pz2*g13+pxz*g33)+2*da1113dz*(pxz*g11+px2*g13);
		snfact = sxfactor*c-szfactor*s;
		
		/*computing ray velocities and derivatives*/
		dxm =  (a1111*pxt*g11+(a1133+a1313)*pzt*g13+a3313*pzt*g33
			+a1113*(pzt*g11+2*pxt*g13)+a1313*g33*pxt);
		dzm =  (a3333*pzt*g33+(a1133+a1313)*pxt*g13+a1113*pxt*g11
			+a3313*(pxt*g33+2*pzt*g13)+a1313*g11*pzt);
		dpxm = -0.5*sxfactor;
		dpzm = -0.5*szfactor;

		dgamma11dpx = 2*a1111*pxt+2*a1113*pzt;
		dgamma11dpz = 2*a1313*pzt+2*a1113*pxt;
		dgamma33dpx = 2*a1313*pxt+2*a3313*pzt;
		dgamma33dpz = 2*a3333*pzt+2*a3313*pxt;
		dgamma13dpx= (a1133+a1313)*pzt+2*a1113*pxt;
		dgamma13dpz= (a1133+a1313)*pxt+2*a3313*pzt;
		dgamma11dn = da1111dn*px2+ da1313dn*pz2 +2*da1113dn*pxz;
		dgamma33dn = da3333dn*pz2 + da1313dn*px2+2*da3313dn*pxz;
		dgamma13dn = (da1133dn+da1313dn)*pxz+ da1113dn*px2+ da3313dn*pz2;
		dg11dpx = -(gamma33-1)*(dgamma11dpx+dgamma33dpx-4*dxm)*den*den+
			(dgamma33dpx-2*dxm)*den;
		dg11dpz = -(gamma33-1)*(dgamma11dpz+dgamma33dpz-4*dzm)*den*den+
			(dgamma33dpz-2*dzm)*den;
		dg33dpx = -(gamma11-1)*(dgamma11dpx+dgamma33dpx-4*dxm)*den*den+
			(dgamma11dpx-2*dxm)*den;
		dg33dpz = -(gamma11-1)*(dgamma11dpz+dgamma33dpz-4*dzm)*den*den+
			(dgamma11dpz-2*dzm)*den;
		dg13dpx = gamma13*(dgamma11dpx+dgamma33dpx-4*dxm)*den*den-
			dgamma13dpx*den;
		dg13dpz = gamma13*(dgamma11dpz+dgamma33dpz-4*dzm)*den*den-
			dgamma13dpz*den;
		dg11dn = -(gamma33-1)*(dgamma11dn+dgamma33dn-2*snfact)*den*den+
			(dgamma33dn-snfact)*den;
		dg33dn = -(gamma11-1)*(dgamma11dn+dgamma33dn-2*snfact)*den*den+
			(dgamma11dn-snfact)*den;
		dg13dn = gamma13*(dgamma11dn+dgamma33dn-2*snfact)*den*den-
			dgamma13dn*den;
		ddxdpx=   a1111*pxt*dg11dpx+(a1133+a1313)*pzt*dg13dpx+
			a3313*pzt*dg33dpx+a1113*(pzt*dg11dpx+2*pxt*dg13dpx)
			+a1313*dg33dpx*pxt;
		ddzdpz= a3333*pzt*dg33dpz+(a1133+a1313)*pxt*dg13dpz+
			a1113*pxt*dg11dpz+a3313*(pxt*dg33dpz+2*pzt*dg13dpz)+
			a1313*dg11dpz*pzt;
		ddxdpz= a1111*pxt*dg11dpz+(a1133+a1313)*pzt*dg13dpz+
			a3313*pzt*dg33dpz+a1113*(pzt*dg11dpz+2*pxt*dg13dpz)+
			a1313*dg33dpz*pxt;
		dsnfactordn = da1111dn*px2*dg11dn+da3333dn*pz2*dg33dn+
			2*(da1133dn+da1313dn)*pxz*dg13dn+da1313dn*(px2*dg33dn+pz2*dg11dn)+
			2*da3313dn*(pz2*dg13dn+pxz*dg33dn)+2*da1113dn*(pxz*dg11dn+px2*dg13dn);
		ddxdn =  (a1111*pxt*dg11dn+(a1133+a1313)*pzt*dg13dn+a3313*pzt*dg33dn
			+a1113*(pzt*dg11dn+2*pxt*dg13dn)+a1313*dg33dn*pxt);
		ddzdn =  (a3333*pzt*dg33dn+(a1133+a1313)*pxt*dg13dn+a1113*pxt*dg11dn
			+a3313*(pxt*dg33dn+2*pzt*dg13dn)+a1313*dg11dn*pzt);

		
		dcdp1  = a1111*g11+a1313*g33+2*a1113*g13+ddxdpx-dxm*dxm;
		dcdp3  = a3333*g33+a1313*g11+2*a3313*g13+ddzdpz-dzm*dzm;
		dcdp13 = a1133*g13+a1313*g13+a1113*g11+a3313*g33+ddxdpz-dxm*dzm;
		ddcdqq = dcdp1*cc-2.0*dcdp13*s*c+dcdp3*ss;
		dxdn   =  (da1111dn*pxt*g11+(da1133dn+da1313dn)*pzt*g13+
			da3313dn*pzt*g33+da1113dn*(pzt*g11+2*pxt*g13)+
			da1313dn*g33*pxt);
		dzdn   =  (da3333dn*pzt*g33+(da1133dn+da1313dn)*pxt*g13+
			da1113dn*pxt*g11+da3313dn*(pxt*g33+2*pzt*g13)+
			da1313dn*g11*pzt);
		ddcdpn = dxdn*c-dzdn*s-.5*dxm*sxfactor*cc+
			.5*(dxm*szfactor+dzm*sxfactor)*s*c-.5*dzm*szfactor*ss
			+ddxdn*c-ddzdn*s;
		snfactor = dda1111dndn*px2*g11+dda3333dndn*pz2*g33+
			2*(dda1133dndn+dda1313dndn)*pxz*g13+
			dda1313dndn*(px2*g33+pz2*g11)+
			2*dda3313dndn*(pz2*g13+pxz*g33)+
			2*dda1113dndn*(pxz*g11+px2*g13);
		ddcdnn = 0.5*snfactor-.25*sxfactor*sxfactor*cc+
			.5*sxfactor*szfactor*s*c-.25*szfactor*szfactor*ss
			+.5*dsnfactordn;

		dp1m = -ddcdnn*q1t-ddcdpn*p1t;
		dq1m = ddcdqq*p1t+ddcdpn*q1t;
		dp2m = -ddcdnn*q2t-ddcdpn*p2t;
		dq2m = ddcdqq*p2t+ddcdpn*q2t;
		xt = x+hhalf*dx;
		zt = z+hhalf*dz;
		pxt = px+h*dpxm;
		pzt = pz+h*dpzm;
		p1t = p1+h*dp1m;
		q1t = q1+h*dq1m;
		p2t = p2+h*dp2m;
		q2t = q2+h*dq2m;
		dxm += dxt;
		dzm += dzt;
		dpxm += dpxt;
		dpzm += dpzt;
		dp1m += dp1t;
		dq1m += dq1t;
		dp2m += dp2t;
		dq2m += dq2t;
		vp  = 1/sqrt(pxt*pxt+pzt*pzt);
		s   = pxt*vp;
		c   = pzt*vp;
		ss  = s*s;
		cc  = c*c;
		
		vel2Interp(a11112,xt,zt,&a1111,&da1111dx,&da1111dz,&dda1111dxdx,
			&dda1111dxdz,&dda1111dzdz);
		da1111dn    = da1111dx*c-da1111dz*s;
		dda1111dndn = dda1111dxdx*cc-2.0*dda1111dxdz*s*c+dda1111dzdz*ss;

		vel2Interp(a33332,xt,zt,&a3333,&da3333dx,&da3333dz,&dda3333dxdx,
		&dda3333dxdz,&dda3333dzdz);
		da3333dn    = da3333dx*c-da3333dz*s;
		dda3333dndn = dda3333dxdx*cc-2.0*dda3333dxdz*s*c+dda3333dzdz*ss;
	
		vel2Interp(a11332,xt,zt,&a1133,&da1133dx,&da1133dz,&dda1133dxdx,
			&dda1133dxdz,&dda1133dzdz);
		da1133dn    = da1133dx*c-da1133dz*s;
		dda1133dndn = dda1133dxdx*cc-2.0*dda1133dxdz*s*c+dda1133dzdz*ss;

		vel2Interp(a13132,xt,zt,&a1313,&da1313dx,&da1313dz,&dda1313dxdx,
			&dda1313dxdz,&dda1313dzdz);
		da1313dn    = da1313dx*c-da1313dz*s;
		dda1313dndn = dda1313dxdx*cc-2.0*dda1313dxdz*s*c+dda1313dzdz*ss;

		vel2Interp(a11132,xt,zt,&a1113,&da1113dx,&da1113dz,&dda1113dxdx,
			&dda1113dxdz,&dda1113dzdz);
		da1113dn    = da1113dx*c-da1113dz*s;
		dda1113dndn = dda1113dxdx*cc-2.0*dda1113dxdz*s*c+dda1113dzdz*ss;

		vel2Interp(a33132,xt,zt,&a3313,&da3313dx,&da3313dz,&dda3313dxdx,
			&dda3313dxdz,&dda3313dzdz);
		da3313dn    = da3313dx*c-da3313dz*s;
		dda3313dndn = dda3313dxdx*cc-2.0*dda3313dxdz*s*c+dda3313dzdz*ss;
		
		/* step 4 of 4th-order Runge-Kutta */
		px2   = pxt*pxt;
		pz2   = pzt*pzt;
		pxz   = pxt*pzt;

		/*anisotropy parameters*/
		gamma11 = a1111*px2+ a1313*pz2 +2*a1113*pxz;
		gamma33 = a3333*pz2 + a1313*px2+2*a3313*pxz;
		gamma13 = (a1133+a1313)*pxz+ a1113*px2+ a3313*pz2;
		den     = 1/(gamma11+gamma33-2);
		g11     = (gamma33-1)*den;
		g33     = (gamma11-1)*den;
		g13     = -gamma13*den;
		sxfactor = da1111dx*px2*g11+da3333dx*pz2*g33+
			2*(da1133dx+da1313dx)*pxz*g13+da1313dx*(px2*g33+pz2*g11)+
			2*da3313dx*(pz2*g13+pxz*g33)+2*da1113dx*(pxz*g11+px2*g13);
		szfactor = da1111dz*px2*g11+da3333dz*pz2*g33+
			2*(da1133dz+da1313dz)*pxz*g13+da1313dz*(px2*g33+pz2*g11)+
			2*da3313dz*(pz2*g13+pxz*g33)+2*da1113dz*(pxz*g11+px2*g13);
		snfact = sxfactor*c-szfactor*s;
		
		/*computing ray velocities and derivatives*/
		dxt =  (a1111*pxt*g11+(a1133+a1313)*pzt*g13+a3313*pzt*g33
			+a1113*(pzt*g11+2*pxt*g13)+a1313*g33*pxt);
		dzt =  (a3333*pzt*g33+(a1133+a1313)*pxt*g13+a1113*pxt*g11
			+a3313*(pxt*g33+2*pzt*g13)+a1313*g11*pzt);
		dpxt = -0.5*sxfactor;
		dpzt = -0.5*szfactor;

		dgamma11dpx = 2*a1111*pxt+2*a1113*pzt;
		dgamma11dpz = 2*a1313*pzt+2*a1113*pxt;
		dgamma33dpx = 2*a1313*pxt+2*a3313*pzt;
		dgamma33dpz = 2*a3333*pzt+2*a3313*pxt;
		dgamma13dpx= (a1133+a1313)*pzt+2*a1113*pxt;
		dgamma13dpz= (a1133+a1313)*pxt+2*a3313*pzt;
		dgamma11dn = da1111dn*px2+ da1313dn*pz2 +2*da1113dn*pxz;
		dgamma33dn = da3333dn*pz2 + da1313dn*px2+2*da3313dn*pxz;
		dgamma13dn = (da1133dn+da1313dn)*pxz+ da1113dn*px2+ da3313dn*pz2;
		dg11dpx = -(gamma33-1)*(dgamma11dpx+dgamma33dpx-4*dxt)*den*den+
			(dgamma33dpx-2*dxt)*den;
		dg11dpz = -(gamma33-1)*(dgamma11dpz+dgamma33dpz-4*dzt)*den*den+
			(dgamma33dpz-2*dzt)*den;
		dg33dpx = -(gamma11-1)*(dgamma11dpx+dgamma33dpx-4*dxt)*den*den+
			(dgamma11dpx-2*dxt)*den;
		dg33dpz = -(gamma11-1)*(dgamma11dpz+dgamma33dpz-4*dzt)*den*den+
			(dgamma11dpz-2*dzt)*den;
		dg13dpx = gamma13*(dgamma11dpx+dgamma33dpx-4*dxt)*den*den-
			dgamma13dpx*den;
		dg13dpz = gamma13*(dgamma11dpz+dgamma33dpz-4*dzt)*den*den-
			dgamma13dpz*den;
		dg11dn = -(gamma33-1)*(dgamma11dn+dgamma33dn-2*snfact)*den*den+
			(dgamma33dn-snfact)*den;
		dg33dn = -(gamma11-1)*(dgamma11dn+dgamma33dn-2*snfact)*den*den+
			(dgamma11dn-snfact)*den;
		dg13dn = gamma13*(dgamma11dn+dgamma33dn-2*snfact)*den*den-
			dgamma13dn*den;
		ddxdpx=   a1111*pxt*dg11dpx+(a1133+a1313)*pzt*dg13dpx+
			a3313*pzt*dg33dpx+a1113*(pzt*dg11dpx+2*pxt*dg13dpx)
			+a1313*dg33dpx*pxt;
		ddzdpz= a3333*pzt*dg33dpz+(a1133+a1313)*pxt*dg13dpz+
			a1113*pxt*dg11dpz+a3313*(pxt*dg33dpz+2*pzt*dg13dpz)+
			a1313*dg11dpz*pzt;
		ddxdpz= a1111*pxt*dg11dpz+(a1133+a1313)*pzt*dg13dpz+
			a3313*pzt*dg33dpz+a1113*(pzt*dg11dpz+2*pxt*dg13dpz)+
			a1313*dg33dpz*pxt;
		dsnfactordn = da1111dn*px2*dg11dn+da3333dn*pz2*dg33dn+
			2*(da1133dn+da1313dn)*pxz*dg13dn+da1313dn*(px2*dg33dn+pz2*dg11dn)+
			2*da3313dn*(pz2*dg13dn+pxz*dg33dn)+2*da1113dn*(pxz*dg11dn+px2*dg13dn);
		ddxdn =  (a1111*pxt*dg11dn+(a1133+a1313)*pzt*dg13dn+a3313*pzt*dg33dn
			+a1113*(pzt*dg11dn+2*pxt*dg13dn)+a1313*dg33dn*pxt);
		ddzdn =  (a3333*pzt*dg33dn+(a1133+a1313)*pxt*dg13dn+a1113*pxt*dg11dn
			+a3313*(pxt*dg33dn+2*pzt*dg13dn)+a1313*dg11dn*pzt);

		
		dcdp1  = a1111*g11+a1313*g33+2*a1113*g13+ddxdpx-dxt*dxt;
		dcdp3  = a3333*g33+a1313*g11+2*a3313*g13+ddzdpz-dzt*dzt;
		dcdp13 = a1133*g13+a1313*g13+a1113*g11+a3313*g33+ddxdpz-dxt*dzt;
		ddcdqq = dcdp1*cc-2.0*dcdp13*s*c+dcdp3*ss;
		dxdn   =  (da1111dn*pxt*g11+(da1133dn+da1313dn)*pzt*g13+
			da3313dn*pzt*g33+da1113dn*(pzt*g11+2*pxt*g13)+
			da1313dn*g33*pxt);
		dzdn   =  (da3333dn*pzt*g33+(da1133dn+da1313dn)*pxt*g13+
			da1113dn*pxt*g11+da3313dn*(pxt*g33+2*pzt*g13)+
			da1313dn*g11*pzt);
		ddcdpn = dxdn*c-dzdn*s-.5*dxt*sxfactor*cc+
			.5*(dxt*szfactor+dzt*sxfactor)*s*c-.5*dzt*szfactor*ss
			+ddxdn*c-ddzdn*s;
		snfactor = dda1111dndn*px2*g11+dda3333dndn*pz2*g33+
			2*(dda1133dndn+dda1313dndn)*pxz*g13+
			dda1313dndn*(px2*g33+pz2*g11)+
			2*dda3313dndn*(pz2*g13+pxz*g33)+
			2*dda1113dndn*(pxz*g11+px2*g13);
		ddcdnn = 0.5*snfactor-.25*sxfactor*sxfactor*cc+
			.5*sxfactor*szfactor*s*c-.25*szfactor*szfactor*ss
			+.5*dsnfactordn;

		dp1t = -ddcdnn*q1t-ddcdpn*p1t;
		dq1t = ddcdqq*p1t+ddcdpn*q1t;
		dp2t = -ddcdnn*q2t-ddcdpn*p2t;
		dq2t = ddcdqq*p2t+ddcdpn*q2t;
		dxx  = hsixth*(dx+dxt+2.0*dxm);
		dzz  = hsixth*(dz+dzt+2.0*dzm);
		x += dxx;
		z += dzz;
		px += hsixth*(dpx+dpxt+2.0*dpxm);
		pz += hsixth*(dpz+dpzt+2.0*dpzm);
		p1 += hsixth*(dp1+dp1t+2.0*dp1m);
		q1 += hsixth*(dq1+dq1t+2.0*dq1m);
		p2 += hsixth*(dp2+dp2t+2.0*dp2m);
		q2 += hsixth*(dq2+dq2t+2.0*dq2m);
		vp  = 1/sqrt(px*px+pz*pz);
		s   = px*vp;
		c   = pz*vp;
		ss  = s*s;
		cc  = c*c;
		
		vel2Interp(a11112,x,z,&a1111,&da1111dx,&da1111dz,&dda1111dxdx,
			&dda1111dxdz,&dda1111dzdz);
		da1111dn    = da1111dx*c-da1111dz*s;
		dda1111dndn = dda1111dxdx*cc-2.0*dda1111dxdz*s*c+dda1111dzdz*ss;

		vel2Interp(a33332,x,z,&a3333,&da3333dx,&da3333dz,&dda3333dxdx,
		&dda3333dxdz,&dda3333dzdz);
		da3333dn    = da3333dx*c-da3333dz*s;
		dda3333dndn = dda3333dxdx*cc-2.0*dda3333dxdz*s*c+dda3333dzdz*ss;
	
		vel2Interp(a11332,x,z,&a1133,&da1133dx,&da1133dz,&dda1133dxdx,
			&dda1133dxdz,&dda1133dzdz);
		da1133dn    = da1133dx*c-da1133dz*s;
		dda1133dndn = dda1133dxdx*cc-2.0*dda1133dxdz*s*c+dda1133dzdz*ss;

		vel2Interp(a13132,x,z,&a1313,&da1313dx,&da1313dz,&dda1313dxdx,
			&dda1313dxdz,&dda1313dzdz);
		da1313dn    = da1313dx*c-da1313dz*s;
		dda1313dndn = dda1313dxdx*cc-2.0*dda1313dxdz*s*c+dda1313dzdz*ss;

		vel2Interp(a11132,x,z,&a1113,&da1113dx,&da1113dz,&dda1113dxdx,
			&dda1113dxdz,&dda1113dzdz);
		da1113dn    = da1113dx*c-da1113dz*s;
		dda1113dndn = dda1113dxdx*cc-2.0*dda1113dxdz*s*c+dda1113dzdz*ss;

		vel2Interp(a33132,x,z,&a3313,&da3313dx,&da3313dz,&dda3313dxdx,
			&dda3313dxdz,&dda3313dzdz);
		da3313dn    = da3313dx*c-da3313dz*s;
		dda3313dndn = dda3313dxdx*cc-2.0*dda3313dxdz*s*c+dda3313dzdz*ss;


		/* update time */
		t  += dt;

		/* update kmah index */
                if ((q2<=0.0 && q2old>0.0) || (q2>=0.0 && q2old<0.0)) kmah++;

		/* save ray parameters */
		rs[it].t = t;
		rs[it].x = x;
		rs[it].z = z;
		rs[it].c = c;
		rs[it].s = s;
		rs[it].q1 = q1;
		rs[it].p1 = p1;
		rs[it].q2 = q2;
		rs[it].p2 = p2;
		rs[it].kmah = kmah;
		rs[it].v = vp;
		rs[it].dvdx = .5*da3333dx*vp/a3333;
		rs[it].dvdz = .5*da3333dz*vp/a3333;
	//printf("%d %f %f %f %f %f %f %f \n",it,t,x,z,fx,fz,lx,lz);
         fprintf(fpraypath,"%g %g %g\n",x,z,t); 
	}


	/* free velocity interpolator */
	vel2Free(a11112);
	vel2Free(a33332);
	vel2Free(a11332);
	vel2Free(a13132);
	vel2Free(a11132);
	vel2Free(a33132);

	/* return ray */
	ray->nrs = it;
	ray->rs = rs;
	ray->nc = 0;
	ray->c = NULL;
	return ray;
}

void freeRay (Ray *ray)
/*****************************************************************************
Free a ray.
******************************************************************************
Input:
ray		ray to be freed
*****************************************************************************/
{
	if (ray->c!=NULL) free1((void*)ray->c);
	free1((void*)ray->rs);
	free1((void*)ray);
}

/*****************************************************************************
Functions to support interpolation of velocity and its derivatives.
******************************************************************************
Functions:
vel2Alloc	allocate and initialize an interpolator for v(x,z)
vel2Interp	interpolate v(x,z) and its derivatives
******************************************************************************
Notes:
Interpolation is performed by piecewise cubic Hermite polynomials,
so that velocity and first derivatives are continuous.  Therefore,
velocity v, first derivatives dv/dx and dv/dz, and the mixed
derivative ddv/dxdz are continuous.  However, second derivatives
ddv/dxdx and ddv/dzdz are discontinuous.
*****************************************************************************
Technical Reference:
	Hale, D., 1992, Migration by the Kirchhoff, 
	slant stack, and Gaussian beam methods:
	Colorado School of Mines.
*****************************************************************************
 Credits: 	CWP: Dave Hale
*****************************************************************************/

/* number of pre-computed, tabulated interpolators */
#define NTABLE 101

/* length of each interpolator in table (4 for piecewise cubic) */
#define LTABLE 4

/* table of pre-computed interpolators, for 0th, 1st, and 2nd derivatives */
static float tbl[3][NTABLE][LTABLE];

/* constants */
static int ix=1-LTABLE/2-LTABLE,iz=1-LTABLE/2-LTABLE;
static float ltable=LTABLE,ntblm1=NTABLE-1;

/* indices for 0th, 1st, and 2nd derivatives */
static int kx[6]={0,1,0,2,1,0};
static int kz[6]={0,0,1,0,1,2};

/* function to build interpolator tables; sets tabled=1 when built */
static void buildTables (void);
static int tabled=0;

/* interpolator for velocity function v(x,z) of two variables */
typedef struct Vel2Struct {
	int nx;		/* number of x samples */
	int nz;		/* number of z samples */
	int nxm;	/* number of x samples minus LTABLE */
	int nzm;	/* number of x samples minus LTABLE */
	float xs,xb,zs,zb,sx[3],sz[3],**vu;
} Vel2;

void* vel2Alloc (int nx, float dx, float fx,
	int nz, float dz, float fz, float **v)
/*****************************************************************************
Allocate and initialize an interpolator for v(x,z) and its derivatives.
******************************************************************************
Input:
nx		number of x samples
dx		x sampling interval
fx		first x sample
nz		number of z samples
dz		z sampling interval
fz		first z sample
v		array[nx][nz] of uniformly sampled v(x,z)

*****************************************************************************
Returned:	pointer to interpolator
*****************************************************************************/
{
	Vel2 *vel2;

	/* allocate space */
	vel2 = (Vel2*)alloc1(1,sizeof(Vel2));

	/* set state variables used for interpolation 
                设置状态变量用于插值         */
	vel2->nx = nx;
	vel2->nxm = nx-LTABLE;
	vel2->xs = 1.0/dx;
	vel2->xb = ltable-fx*vel2->xs;
	vel2->sx[0] = 1.0;
	vel2->sx[1] = vel2->xs;
	vel2->sx[2] = vel2->xs*vel2->xs;
	vel2->nz = nz;
	vel2->nzm = nz-LTABLE;
	vel2->zs = 1.0/dz;
	vel2->zb = ltable-fz*vel2->zs;
	vel2->sz[0] = 1.0;
	vel2->sz[1] = vel2->zs;
	vel2->sz[2] = vel2->zs*vel2->zs;
	vel2->vu = v;
	
	/* if necessary, build interpolator coefficient tables 
                如果有必要,构建插值系数表          */
	if (!tabled) buildTables();

	return vel2;
}

void vel2Free (void *vel2)
/*****************************************************************************
Free an interpolator for v(x,z) and its derivatives.
******************************************************************************
Input:
vel2		pointer to interpolator as returned by vel2Alloc()
*****************************************************************************/
{
	free1(vel2);
}

void vel2Interp (void *vel2, float x, float z,
	float *v, float *vx, float *vz, float *vxx, float *vxz, float *vzz)
/*****************************************************************************
Interpolation of a velocity function v(x,z) and its derivatives.
           速度插值
******************************************************************************
Input:
vel2		pointer to interpolator as returned by vel2Alloc()
x		x coordinate at which to interpolate v(x,z) and derivatives
z		z coordinate at which to interpolate v(x,z) and derivatives
*****************************************************************************
Output:
v		v(x,z)
vx		dv/dx
vz		dv/dz
vxx		ddv/dxdx
vxz		ddv/dxdz
vzz		ddv/dzdz
*****************************************************************************/
{
	Vel2 *v2=vel2;
	int nx=v2->nx,nz=v2->nz,nxm=v2->nxm,nzm=v2->nzm;
	float xs=v2->xs,xb=v2->xb,zs=v2->zs,zb=v2->zb,
		*sx=v2->sx,*sz=v2->sz,**vu=v2->vu;
	int i,jx,lx,mx,jz,lz,mz,jmx,jmz,mmx,mmz;
	float ax,bx,*px,az,bz,*pz,sum,vd[6];

	/* determine offsets into vu and interpolation coefficients */
	ax = xb+x*xs;
	jx = (int)ax;
	bx = ax-jx;
	lx = (bx>=0.0)?bx*ntblm1+0.5:(bx+1.0)*ntblm1-0.5;
	lx *= LTABLE;
	mx = ix+jx;
	az = zb+z*zs;
	jz = (int)az;
	bz = az-jz;
	lz = (bz>=0.0)?bz*ntblm1+0.5:(bz+1.0)*ntblm1-0.5;
	lz *= LTABLE;
	mz = iz+jz;
	
	/* if totally within input array, use fast method */
	if (mx>=0 && mx<=nxm && mz>=0 && mz<=nzm) {
		for (i=0; i<6; ++i) {
			px = &(tbl[kx[i]][0][0])+lx;
			pz = &(tbl[kz[i]][0][0])+lz;
			vd[i] = sx[kx[i]]*sz[kz[i]]*(
				vu[mx][mz]*px[0]*pz[0]+
				vu[mx][mz+1]*px[0]*pz[1]+
				vu[mx][mz+2]*px[0]*pz[2]+
				vu[mx][mz+3]*px[0]*pz[3]+
				vu[mx+1][mz]*px[1]*pz[0]+
				vu[mx+1][mz+1]*px[1]*pz[1]+
				vu[mx+1][mz+2]*px[1]*pz[2]+
				vu[mx+1][mz+3]*px[1]*pz[3]+
				vu[mx+2][mz]*px[2]*pz[0]+
				vu[mx+2][mz+1]*px[2]*pz[1]+
				vu[mx+2][mz+2]*px[2]*pz[2]+
				vu[mx+2][mz+3]*px[2]*pz[3]+
				vu[mx+3][mz]*px[3]*pz[0]+
				vu[mx+3][mz+1]*px[3]*pz[1]+
				vu[mx+3][mz+2]*px[3]*pz[2]+
				vu[mx+3][mz+3]*px[3]*pz[3]);
		}
		
	/* else handle end effects with constant extrapolation */
	} else {
		for (i=0; i<6; ++i) {
			px = &(tbl[kx[i]][0][0])+lx;
			pz = &(tbl[kz[i]][0][0])+lz;
			for (jx=0,jmx=mx,sum=0.0; jx<4; ++jx,++jmx) {
				mmx = jmx;
				if (mmx<0) mmx = 0;
				else if (mmx>=nx) mmx = nx-1;
				for (jz=0,jmz=mz; jz<4; ++jz,++jmz) {
					mmz = jmz;
					if (mmz<0) mmz = 0;
					else if (mmz>=nz) mmz = nz-1;
					sum += vu[mmx][mmz]*px[jx]*pz[jz];
				}
			}
			vd[i] = sx[kx[i]]*sz[kz[i]]*sum;
		}
	}

	/* set output variables */
	*v = vd[0];
	*vx = vd[1];
	*vz = vd[2];
	*vxx = vd[3];
	*vxz = vd[4];
	*vzz = vd[5];
}

/* hermite polynomials */
static float h00 (float x) {return 2.0*x*x*x-3.0*x*x+1.0;}
static float h01 (float x) {return 6.0*x*x-6.0*x;}
static float h02 (float x) {return 12.0*x-6.0;}
static float h10 (float x) {return -2.0*x*x*x+3.0*x*x;}
static float h11 (float x) {return -6.0*x*x+6.0*x;}
static float h12 (float x) {return -12.0*x+6.0;}
static float k00 (float x) {return x*x*x-2.0*x*x+x;}
static float k01 (float x) {return 3.0*x*x-4.0*x+1.0;}
static float k02 (float x) {return 6.0*x-4.0;}
static float k10 (float x) {return x*x*x-x*x;}
static float k11 (float x) {return 3.0*x*x-2.0*x;}
static float k12 (float x) {return 6.0*x-2.0;}

/* function to build interpolation tables */
static void buildTables(void)
{
	int i;
	float x;
	
	/* tabulate interpolator for 0th derivative */
	for (i=0; i<NTABLE; ++i) {
		x = (float)i/(NTABLE-1.0);
		tbl[0][i][0] = -0.5*k00(x);
		tbl[0][i][1] = h00(x)-0.5*k10(x);
		tbl[0][i][2] = h10(x)+0.5*k00(x);
		tbl[0][i][3] = 0.5*k10(x);
		tbl[1][i][0] = -0.5*k01(x);
		tbl[1][i][1] = h01(x)-0.5*k11(x);
		tbl[1][i][2] = h11(x)+0.5*k01(x);
		tbl[1][i][3] = 0.5*k11(x);
		tbl[2][i][0] = -0.5*k02(x);
		tbl[2][i][1] = h02(x)-0.5*k12(x);
		tbl[2][i][2] = h12(x)+0.5*k02(x);
		tbl[2][i][3] = 0.5*k12(x);
	}
	
	/* remember that tables have been built */
	tabled = 1;
}

#ifdef TEST2
#include "cwp.h"
#define NZ 2
#define NX 2
#define NXOUT 21
#define NZOUT 21
main()
{
	int nx=NX,nz=NZ,nxout=NXOUT,nzout=NZOUT,i,j;
	float dx=2.0,fx=-1.0,dxout=0.2,fxout=-2.0;
	float dz=4.0,fz=-2.0,dzout=0.4,fzout=-4.0;
	float x,z,v,vx,vz,vxx,vxz,vzz,**vu,**vo;
	void *vel2;

	vu = alloc2float(nz,nx);
	vo = alloc2float(nzout,nxout);

	vu[0][0] = 1.0;
	vu[1][0] = 2.0;
	vu[0][1] = 1.0;
	vu[1][1] = 2.0;

	vel2 = vel2Alloc(nx,dx,fx,nz,dz,fz,vu);

	for (i=0; i<nxout; ++i) {
		x = fxout+i*dxout;
		for (j=0; j<nzout; ++j) {
			z = fzout+j*dzout;
			vel2Interp(vel2,x,z,&v,&vx,&vz,&vxx,&vxz,&vzz);
			vo[i][j] = vz;
		}
	}
	vel2Free(vel2);
	fwrite(vo[0],sizeof(float),nxout*nzout,stdout);
}
#endif
