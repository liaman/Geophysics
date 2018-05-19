#include<gtk/gtk.h>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
//#include "model_VTI.c"
void model(int nx, int nz,float dx,float dz,int npd,int mm,
           const char FNv[],const char FNe[],const char FNd[],
           float favg,int ns,int fs,int ds,int zs,
           const char FNshot[],const char FNsnap[])
{

  printf("##### model start #####\n");
  printf("#  nx=%2d, dx=%.2f, npd=%d\n",nx,dx,npd);
  printf("#  nz=%2d, dz=%.2f, mm=%d\n",nz,dz,mm);
  printf("#     vel=<%s>\n",FNv);
  printf("#  epsilu=<%s>\n",FNe);
  printf("#    deta=<%s>\n",FNd);
  printf("#  favg=%.2f\n",favg);
  printf("#  ns=%3d\n",ns);
  printf("#  fs=%3d\n",fs);
  printf("#  ds=%3d\n",ds);
  printf("#  zs=%3d\n",zs);
  printf("#    shot=<%s>\n",FNshot);
  printf("#    snap=<%s>\n",FNsnap);
}
GtkWidget *_nx, *_nz, *_dx, *_dz, *_npd, *_mm;
GtkWidget *_vel, *_epsilu, *_deta, *_shot, *_snap;
GtkWidget *_favg, *_fs, *_ds, *_ns, *_zs;

gint nx, nz, npd, mm, fs, ds, ns, zs;
gfloat dx, dz, favg;

void debug(GtkWidget *window,gpointer data)
{  
    gtk_main_quit();
}
 
void write_txt(GtkWidget *widget,gpointer data)
{

    const char*str_nx=gtk_entry_get_text(GTK_ENTRY((GtkWidget *)_nx));
    const char*str_nz=gtk_entry_get_text(GTK_ENTRY((GtkWidget *)_nz));
    const char*str_dx=gtk_entry_get_text(GTK_ENTRY((GtkWidget *)_dx));
    const char*str_dz=gtk_entry_get_text(GTK_ENTRY((GtkWidget *)_dz));
    const char*str_npd=gtk_entry_get_text(GTK_ENTRY((GtkWidget *)_npd));
    const char*str_mm=gtk_entry_get_text(GTK_ENTRY((GtkWidget *)_mm));
    const char*str_vel=gtk_entry_get_text(GTK_ENTRY((GtkWidget *)_vel));
    const char*str_epsilu=gtk_entry_get_text(GTK_ENTRY((GtkWidget *)_epsilu));
    const char*str_deta=gtk_entry_get_text(GTK_ENTRY((GtkWidget *)_deta));
    const char*str_favg=gtk_entry_get_text(GTK_ENTRY((GtkWidget *)_favg));
    const char*str_fs=gtk_entry_get_text(GTK_ENTRY((GtkWidget *)_fs));
    const char*str_ns=gtk_entry_get_text(GTK_ENTRY((GtkWidget *)_ns));
    const char*str_ds=gtk_entry_get_text(GTK_ENTRY((GtkWidget *)_ds));
    const char*str_zs=gtk_entry_get_text(GTK_ENTRY((GtkWidget *)_zs));
    const char*str_shot=gtk_entry_get_text(GTK_ENTRY((GtkWidget *)_shot));
    const char*str_snap=gtk_entry_get_text(GTK_ENTRY((GtkWidget *)_snap));


    nx=strtod(str_nx,NULL);
    nz=strtod(str_nz,NULL);
    dx=strtod(str_dx,NULL);
    dz=strtod(str_dz,NULL);
    npd=strtod(str_npd,NULL);
    mm=strtod(str_mm,NULL);
    favg=strtod(str_favg,NULL);
    fs=strtod(str_fs,NULL);
    ns=strtod(str_ns,NULL);
    ds=strtod(str_ds,NULL);
    zs=strtod(str_zs,NULL);

    model(nx,nz,dx,dz,npd,mm,str_vel,str_epsilu,str_deta,
          favg,ns,fs,ds,zs,str_shot,str_snap);
}
int main(int argc,char* argv[]){
 
    GtkWidget *window;
    GtkWidget *button;

    GtkWidget *label;
    GtkWidget *table;
   
    gtk_init(&argc,&argv);
       
    window=gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(window),"RT's VTI model");
    gtk_window_set_position(GTK_WINDOW(window),GTK_WIN_POS_CENTER);
    //gtk_window_set_default_size(GTK_WINDOW(window),400,400);
    g_signal_connect(GTK_OBJECT(window),"delete_event",GTK_SIGNAL_FUNC(debug),NULL);
     
    table=gtk_table_new(13,6,FALSE);
    gtk_container_add(GTK_CONTAINER(window),table);
    label=gtk_label_new("Model parameter:");
    gtk_table_attach_defaults(GTK_TABLE(table),label,1,2,0,1);

/***   nx  ***/
    label=gtk_label_new("nx=");
    gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,1,2);
    _nx=gtk_entry_new_with_max_length(10);
    gtk_entry_set_text (GTK_ENTRY (_nx), "100");
    gtk_table_attach_defaults(GTK_TABLE(table),_nx,1,2,1,2);
/***   dx  ***/
    label=gtk_label_new("dx=");
    gtk_table_attach_defaults(GTK_TABLE(table),label,2,3,1,2);
    _dx=gtk_entry_new_with_max_length(10);
    gtk_entry_set_text (GTK_ENTRY (_dx), "5");
    gtk_table_attach_defaults(GTK_TABLE(table),_dx,3,4,1,2);
/***   npd  ***/
    label=gtk_label_new("npd=");
    gtk_table_attach_defaults(GTK_TABLE(table),label,4,5,1,2);
    _npd=gtk_entry_new_with_max_length(10);
    gtk_entry_set_text (GTK_ENTRY (_npd), "20");
    gtk_table_attach_defaults(GTK_TABLE(table),_npd,5,6,1,2);
/***   nz  ***/
    label=gtk_label_new("nz=");
    gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,2,3);
    _nz=gtk_entry_new_with_max_length(10);
    gtk_entry_set_text (GTK_ENTRY (_nz), "100");
    gtk_table_attach_defaults(GTK_TABLE(table),_nz,1,2,2,3);
/***   dz  ***/
    label=gtk_label_new("dz=");
    gtk_table_attach_defaults(GTK_TABLE(table),label,2,3,2,3);
    _dz=gtk_entry_new_with_max_length(10);
    gtk_entry_set_text (GTK_ENTRY (_dz), "5");
    gtk_table_attach_defaults(GTK_TABLE(table),_dz,3,4,2,3);
/***   mm  ***/
    label=gtk_label_new("mm=");
    gtk_table_attach_defaults(GTK_TABLE(table),label,4,5,2,3);
    _mm=gtk_entry_new_with_max_length(10);
    gtk_entry_set_text (GTK_ENTRY (_mm), "4");
    gtk_table_attach_defaults(GTK_TABLE(table),_mm,5,6,2,3);

    label=gtk_label_new("Input Filename:");
    gtk_table_attach_defaults(GTK_TABLE(table),label,1,2,3,4);
/***   vel  ***/
    label=gtk_label_new("vel:");
    gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,4,5);
    _vel=gtk_entry_new_with_max_length(10);
    gtk_entry_set_text (GTK_ENTRY (_vel), "vel.dat");
    gtk_table_attach_defaults(GTK_TABLE(table),_vel,1,4,4,5);
/***   epsilu  ***/
    label=gtk_label_new("epsilu:");
    gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,5,6);
    _epsilu=gtk_entry_new_with_max_length(10);
    gtk_entry_set_text (GTK_ENTRY (_epsilu), "epsilu.dat");
    gtk_table_attach_defaults(GTK_TABLE(table),_epsilu,1,4,5,6);
/***   deta  ***/
    label=gtk_label_new("deta:");
    gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,6,7);
    _deta=gtk_entry_new_with_max_length(10);
    gtk_entry_set_text (GTK_ENTRY (_deta), "deta.dat");
    gtk_table_attach_defaults(GTK_TABLE(table),_deta,1,4,6,7);

    label=gtk_label_new("ShotParas:");
    gtk_table_attach_defaults(GTK_TABLE(table),label,1,2,7,8);
/***   ns  ***/
    label=gtk_label_new("ns=");
    gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,8,9);
    _ns=gtk_entry_new_with_max_length(10);
    gtk_entry_set_text (GTK_ENTRY (_ns), "1");
    gtk_table_attach_defaults(GTK_TABLE(table),_ns,1,2,8,9);
/***   favg  ***/
    label=gtk_label_new("favg=");
    gtk_table_attach_defaults(GTK_TABLE(table),label,2,3,8,9);
    _favg=gtk_entry_new_with_max_length(10);
    gtk_entry_set_text (GTK_ENTRY (_favg), "20");
    gtk_table_attach_defaults(GTK_TABLE(table),_favg,3,4,8,9);
/***   fs  ***/
    label=gtk_label_new("fs=");
    gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,9,10);
    _fs=gtk_entry_new_with_max_length(10);
    gtk_entry_set_text (GTK_ENTRY (_fs), "1");
    gtk_table_attach_defaults(GTK_TABLE(table),_fs,1,2,9,10);
/***   ds  ***/
    label=gtk_label_new("ds=");
    gtk_table_attach_defaults(GTK_TABLE(table),label,2,3,9,10);
    _ds=gtk_entry_new_with_max_length(10);
    gtk_entry_set_text (GTK_ENTRY (_ds), "0");
    gtk_table_attach_defaults(GTK_TABLE(table),_ds,3,4,9,10);
/***   zs  ***/
    label=gtk_label_new("zs=");
    gtk_table_attach_defaults(GTK_TABLE(table),label,4,5,9,10);
    _zs=gtk_entry_new_with_max_length(10);
    gtk_entry_set_text (GTK_ENTRY (_zs), "1");
    gtk_table_attach_defaults(GTK_TABLE(table),_zs,5,6,9,10);


    label=gtk_label_new("Output Filename:");
    gtk_table_attach_defaults(GTK_TABLE(table),label,1,2,10,11);
/***   shot  ***/
    label=gtk_label_new("shot=");
    gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,11,12);
    _shot=gtk_entry_new_with_max_length(10);
    gtk_entry_set_text (GTK_ENTRY (_shot), "shot.dat");
    gtk_table_attach_defaults(GTK_TABLE(table),_shot,1,4,11,12);
/***   snap  ***/
    label=gtk_label_new("snap=");
    gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,12,13);
    _snap=gtk_entry_new_with_max_length(10);
    gtk_entry_set_text (GTK_ENTRY (_snap), "snap.dat");
    gtk_table_attach_defaults(GTK_TABLE(table),_snap,1,4,12,13);

    label=gtk_label_new("Please ensure all parameters input accuratlly !");
    gtk_table_attach_defaults(GTK_TABLE(table),label,0,10,13,14);
/***   run button  ***/
    button=gtk_button_new_with_label("RUN");
    g_signal_connect(G_OBJECT(button),"clicked",G_CALLBACK(write_txt),NULL);
    gtk_table_attach_defaults(GTK_TABLE(table),button,5,6,12,13);
 
    gtk_widget_show_all(window);
    gtk_main();
    return 0;
}
