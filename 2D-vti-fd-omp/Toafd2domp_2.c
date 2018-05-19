//a#########################################################
//a## 2D Acoustic VTI Medium FD
//a## Ps : P + sv wave and get rid of sv
//a## OpenMP 
//a##
//a##/*a****************************************************/
//a##Function for VTI medium modeling,2017.2.13
//a##
//a## Ps: the function of modeling following:
//a##
//a## du/dt=1/rho*dp/dx ,
//a## dw/dt=1/rho*dq/dz ,
//a## dp/dt=rho*vpx^2*du/dx+rho*vp*vpn*dw/dz ,
//a## dq/dt=rho*vp*vpn*du/dx+rho*vp^2*dw/dz ,
//a## vpx^2=vp^2*(1+2*epsilon);
//a## vpn^2=vp^2*(1+2*delta);
//a##*********a*********************************************/
//a##
//a## programming by Rong Tao
//a#########################################################
#include<stdio.h>
#include<malloc.h>
#include<math.h>
#include<stdlib.h>
#include "/home/Toa/hc/cjbsegy.h"
#include "/home/Toa/hc/fft.c"
#include "/home/Toa/hc/alloc.c"
#include "/home/Toa/hc/complex.c"

#ifdef _OPENMP
#include <omp.h>
#endif

#define pi 3.141592653

float c[4] = {1.196289, -0.0797526, 0.009570313, -0.0006975447};

void pad_vv(int nx, int nz, int npd, float **ee);

void main(int argc,char *argv[])
{

        //omp_set_num_threads(20);

        #ifdef _OPENMP
        /* Testing for OpenMP */
        double start_time, end_time;
        #endif

        int i, j, k, nx, nz, nt, npd, it, sx, sz;
        float dx, dz, dt, pfac, favg, t, **p_cal;



        char FN1[250] = {"v_601_301.dat"};
        char FN2[250] = {"e_601_301.dat"};
        char FN3[250] = {"d_601_301.dat"};
        char FN4[250] = {"shot.dat"};
        char FN5[250] = {"snap.dat"};


        nx = 601;               
        nz = 301;              
        dx = 5.0;   
        dz = 5.0;

        npd = 0;

        favg = 30;
        pfac = 1000.0;

        nt = 4001;    
        dt = 0.0005;

        sx = 100;
        sz = 100;


        #ifdef _OPENMP
        start_time = omp_get_wtime();
        #endif

        p_cal = alloc2float(nt,nx);      
        zero2float(p_cal,nt,nx);


        float **vp0,  **eps,  **del,  **s;
        float **u, **w, **Px, **Pz, **Qx, **Qz, **P, **Q;


        vp0 = alloc2float(nz+2*npd,nx+2*npd);   zero2float(vp0,nz+2*npd,nx+2*npd); 
        eps = alloc2float(nz+2*npd,nx+2*npd);   zero2float(eps,nz+2*npd,nx+2*npd);
        del = alloc2float(nz+2*npd,nx+2*npd);   zero2float(del,nz+2*npd,nx+2*npd);


        FILE *fp1, *fp2, *fp3;
        if ( (fp1=fopen(FN1,"rb")) == NULL ) { printf("readfile error\n"); exit(0); }
        if ( (fp2=fopen(FN2,"rb")) == NULL ) { printf("readfile error\n"); exit(0); }
        if ( (fp3=fopen(FN3,"rb")) == NULL ) { printf("readfile error\n"); exit(0); }
        for ( i = npd; i < nx + npd; i ++ ) {
                for ( j = npd; j < nz + npd; j ++ ) {
                        fread( &vp0[i][j], 4, 1, fp1);//vp0[i][j] = 2000.0;
                        fread( &eps[i][j], 4, 1, fp2);//eps[i][j] = 0.5;
                        fread( &del[i][j], 4, 1, fp3);//del[i][j] = 0.1;
                }
        }
        fclose(fp1);
        fclose(fp2);
        fclose(fp3);


        pad_vv(nx, nz, npd, eps);
        pad_vv(nx, nz, npd, del);
        pad_vv(nx, nz, npd, vp0);

        u = alloc2float(nz+2*npd,nx+2*npd);   zero2float(u,nz+2*npd,nx+2*npd);
        w = alloc2float(nz+2*npd,nx+2*npd);   zero2float(w,nz+2*npd,nx+2*npd);

        P = alloc2float(nz+2*npd,nx+2*npd);  zero2float(P,nz+2*npd,nx+2*npd);
        Q = alloc2float(nz+2*npd,nx+2*npd);  zero2float(Q,nz+2*npd,nx+2*npd);

        Px = alloc2float(nz+2*npd,nx+2*npd); zero2float(Px,nz+2*npd,nx+2*npd);
        Pz = alloc2float(nz+2*npd,nx+2*npd); zero2float(Pz,nz+2*npd,nx+2*npd);
        Qx = alloc2float(nz+2*npd,nx+2*npd); zero2float(Qx,nz+2*npd,nx+2*npd);
        Qz = alloc2float(nz+2*npd,nx+2*npd); zero2float(Qz,nz+2*npd,nx+2*npd);

        s = alloc2float(nz+2*npd,nx+2*npd);  

	 float source;
	 float tdelay, ts;
        int ixs,izs,lx,lz;


        int ii, im;
        float dtx, dtz, xx = 0.0, zz = 0.0;

        dtx = dt / dx;
        dtz = dt / dz;
        tdelay = 1.0 / favg;


        FILE *fp5 = fopen(FN5,"wb");

        for ( it = 0, t = 0.0; it < nt; it++, t += dt ) { 

                if(it%100==0)printf("it = %d\n",it);

                zero2float(s,nz+2*npd,nx+2*npd); 
                source = 0.0;
                if( t <= 2 * tdelay ) {
                        ts = t - tdelay;
                        source = (1-2*pi*pi*(favg*ts*favg*ts))*exp(-(pi*pi*favg*ts*favg*ts));
                        ixs = (int)(sx + 0.5) + npd - 1;
                        izs = (int)(sz + 0.5) + npd - 1;
                        for(j = izs-3;j <= izs+3; j++) { 
                                for(i = ixs-3;i <= ixs+3; i++) {  
                                        lx = i - ixs;
                                        lz = j - izs;
                                        s[i][j] = pfac*source*exp(-lz*lz-lx*lx);
                                }
                        }
                }
                #ifdef _OPENMP
                #pragma omp parallel for default(none) num_threads(20) \
                    private(i, j) \
                    shared(nx, nz, npd, dtx, dtz, P, Q, u, w, c)
                #endif
                for( j =  4 ; j <= (2*npd+nz- 4 -1); j ++ ) { 
                        for( i =  4 ; i <= (2*npd+nx- 4 -1); i++ ) {

                                u[i][j] = u[i][j] 
                                        - dtx * (  c[0]*(P[i+0+1][j] - P[i-0][j]) +
                                                   c[1]*(P[i+1+1][j] - P[i-1][j]) +
                                                   c[2]*(P[i+2+1][j] - P[i-2][j]) +
                                                   c[3]*(P[i+3+1][j] - P[i-3][j]) );

                                w[i][j] = w[i][j] 
                                        - dtz * (  c[0]*(Q[i][j+0+1] - Q[i][j-0]) +
                                                   c[1]*(Q[i][j+1+1] - Q[i][j-1]) +
                                                   c[2]*(Q[i][j+2+1] - Q[i][j-2]) +
                                                   c[3]*(Q[i][j+3+1] - Q[i][j-3]) );
                        }
                }
                #ifdef _OPENMP
                #pragma omp parallel for default(none) num_threads(20) \
                    private(i, j, xx) \
                    shared(nx, nz, npd, dtx, dtz, P, Q, u, w, c, vp0, delta, epsilon, Qx, Qz, Px, Pz)
                #endif
                for( i =  4 ; i <= (2*npd+nx- 4 -1); i ++ ) {
                        for( j =  4 ; j <= (2*npd+nz- 4 -1); j ++ )  {

                                xx = 0.0;
                                zz = 0.0;

                                for( im = 0; im <  4 ; im ++ ) {
                                        xx += c[im]*(u[i+im][j] - u[i-im-1][j]);

                                        zz += c[im]*(w[i][j+im] - w[i][j-im-1]);
                                }
                                Px[i][j] = Px[i][j] - 
                                           (vp0[i][j] * vp0[i][j]) * (1+2*eps[i][j])*dtx*xx;  

                                Pz[i][j] = Pz[i][j] -
                                           (vp0[i][j] * vp0[i][j]) * (pow((1+2*del[i][j]),0.5))*dtz*zz;

                                Qx[i][j] = Qx[i][j] -
                                           (vp0[i][j] * vp0[i][j]) * (pow((1+2*del[i][j]),0.5))*dtx*xx;

                                Qz[i][j] = Qz[i][j] -
                                           (vp0[i][j] * vp0[i][j]) * dtz * zz;

                                P[i][j] = Px[i][j] + Pz[i][j] + s[i][j];

                                Q[i][j] = Qx[i][j] + Qz[i][j] + s[i][j];
                        }
                }

                for( i = npd; i < npd+nx; i++ )   {   
                        p_cal[i-npd][it] = Q[i][npd] + P[i][npd];
                }
                if(it%50==0) {
                        fseek(fp5,(int)(it/50)*(nx)*(nz)*4L,0);
                        for(i=npd;i<nx+npd;i++)
                                for(j=npd;j<nz+npd;j++)
                                        fwrite(&P[i][j],4L,1,fp5);
                }
        }//it

        FILE *fp4;
        fp4 = fopen(FN4,"wb");
        for(i=0;i<nx;i++)
                for(j=0;j<nt;j++)
                        fwrite(&p_cal[i][j],4L,1,fp4);
        fclose(fp4);
        fclose(fp5);

        #ifdef _OPENMP
        end_time = omp_get_wtime();
        printf("Totally %f(s).\n",end_time - start_time);
        #endif

}

void pad_vv(int nx, int nz, int npd, float **ee)
{
        int i, j;

        for ( j = npd; j <= (nz+npd-1); j++ ) 
        for ( i = 0;   i <= npd-1;      i++ )  
                ee[i][j] = ee[npd][j];

        for ( j = npd;    j <= (nz+npd-1);   j++ )
        for ( i = nx+npd; i <= (nx+2*npd-1); i++ )
                ee[i][j] = ee[nx+npd-1][j];

        for ( j = 0; j <= (npd-1);      j++ )
        for ( i = 0; i <= (nx+2*npd-1); i++ )
                ee[i][j] = ee[i][npd];

        for ( j = nz+npd; j <= (nz+2*npd-1); j++ )
        for ( i = 0;      i <= (nx+2*npd-1); i++ )
                ee[i][j] = ee[i][nz+npd-1];	
}
