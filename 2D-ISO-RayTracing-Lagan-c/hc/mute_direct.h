/************* mute direct-wave on one trace  **********/

void mute_direct(float *trace, int nt, float dt, float vel_direct,
     float offset, float t0, int nttaper)
{
    /*** t0--------- direct-wave starting time of zero-offset ****/
    /*** nttaper---- taper point number beyond direct-wave starting time ****/

    int     it, jt, kt, jt1;
    float   t, rc;

    t=fabs(offset)/vel_direct+t0;
    jt=t/dt;
    kt=jt+nttaper;
    jt1=jt-1;
    rc=1.0/nttaper;

    for(it=0;it<jt1;it++) trace[it]=0.0;
    for(it=jt-1;it<kt;it++) trace[it]=trace[it]*rc;
}    
