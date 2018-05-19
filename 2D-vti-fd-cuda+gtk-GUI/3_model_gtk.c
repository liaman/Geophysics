#include<gtk/gtk.h>
#include<stdlib.h>
//#include"model_VTI.c"
//#include"kernels.c"
//###################################model#######################################
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
//###################################parameter of model#########################################
GtkWidget *file,*_vel,*_deta,*_epsilu;
GtkWidget *_nx, *_nz, *_dx, *_dz, *_npd, *_mm;
GtkWidget *_favg, *_fs, *_ds, *_ns, *_zs;
GtkWidget *_shot, *_snap;
gint nx, nz, npd, mm, fs, ds, ns, zs;
gfloat dx, dz, favg;
//###################################file OK return string#########################################
void file_OK( GtkWidget *w, gpointer *data )
{
     gtk_entry_set_text (GTK_ENTRY(data),gtk_file_selection_get_filename (GTK_FILE_SELECTION(file)));
     gtk_widget_destroy (file);
}
//###################################file select windows#########################################
void select_file(GtkWidget *w,gpointer *data)
{
     file= gtk_file_selection_new ("File selection");
     g_signal_connect (G_OBJECT(GTK_FILE_SELECTION(file)->ok_button), "clicked",G_CALLBACK (file_OK),data);
     g_signal_connect_swapped (G_OBJECT(GTK_FILE_SELECTION(file)->cancel_button),"clicked", G_CALLBACK(gtk_widget_destroy), NULL);
     gtk_widget_show (file);
}
//###################################function of label image box#########################################
GtkWidget *xpm_label_box( gchar *xpm_filename, gchar*label_text)
{
     GtkWidget *box;
     GtkWidget *label;
     GtkWidget *image;
     box = gtk_hbox_new (FALSE, 0);
     gtk_container_set_border_width (GTK_CONTAINER (box), 10);
     image = gtk_image_new_from_file (xpm_filename);
     label = gtk_label_new (label_text);
     gtk_box_pack_start (GTK_BOX (box), image, TRUE, TRUE,3);
     gtk_box_pack_start (GTK_BOX (box), label, TRUE,TRUE,3);
     gtk_widget_show (image);
     gtk_widget_show (label);
     return box;
}
//###################################connect to model#########################################
void write_txt(GtkWidget *widget,gpointer data)
{
     const char*str_nx=gtk_entry_get_text(GTK_ENTRY((GtkWidget*)_nx));
     const char*str_nz=gtk_entry_get_text(GTK_ENTRY((GtkWidget*)_nz));
     const char*str_dx=gtk_entry_get_text(GTK_ENTRY((GtkWidget*)_dx));
     const char*str_dz=gtk_entry_get_text(GTK_ENTRY((GtkWidget*)_dz));
     const char*str_npd=gtk_entry_get_text(GTK_ENTRY((GtkWidget*)_npd));
     const char*str_mm=gtk_entry_get_text(GTK_ENTRY((GtkWidget*)_mm));
     const char*str_vel=gtk_entry_get_text(GTK_ENTRY((GtkWidget*)_vel));
     const char*str_epsilu=gtk_entry_get_text(GTK_ENTRY((GtkWidget*)_epsilu));
     const char*str_deta=gtk_entry_get_text(GTK_ENTRY((GtkWidget*)_deta));
     const char*str_favg=gtk_entry_get_text(GTK_ENTRY((GtkWidget*)_favg));
     const char*str_fs=gtk_entry_get_text(GTK_ENTRY((GtkWidget*)_fs));
     const char*str_ns=gtk_entry_get_text(GTK_ENTRY((GtkWidget*)_ns));
     const char*str_ds=gtk_entry_get_text(GTK_ENTRY((GtkWidget*)_ds));
     const char*str_zs=gtk_entry_get_text(GTK_ENTRY((GtkWidget*)_zs));
     const char*str_shot=gtk_entry_get_text(GTK_ENTRY((GtkWidget*)_shot));
     const char*str_snap=gtk_entry_get_text(GTK_ENTRY((GtkWidget*)_snap));
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
     model(nx,nz,dx,dz,npd,mm,str_vel,str_epsilu,str_deta,favg,ns,fs,ds,zs,str_shot,str_snap);
}
//###################################VTI model#########################################
void VTI_model(GtkWidget *widget,gpointer data)
{
     GtkWidget *window;
     GtkWidget *button;
     GtkWidget *label;
     GtkWidget *table;
     GtkWidget *frame;
     GtkWidget *boxV,*boxH,*box;
     window=gtk_window_new(GTK_WINDOW_TOPLEVEL);
     gtk_window_set_title(GTK_WINDOW(window),"VTI model");
     gtk_window_set_position(GTK_WINDOW(window),GTK_WIN_POS_CENTER);
     g_signal_connect(GTK_OBJECT(window),"delete_event",GTK_SIGNAL_FUNC(gtk_widget_destroy),NULL);
     boxV = gtk_vbox_new (FALSE, 0);
     gtk_container_add (GTK_CONTAINER (window), boxV);
     gtk_widget_show (boxV);
     frame = gtk_frame_new (NULL);
     gtk_container_set_border_width (GTK_CONTAINER (frame), 10);
     gtk_box_pack_start (GTK_BOX (boxV), frame, TRUE, TRUE, 0);
     gtk_frame_set_label (GTK_FRAME (frame), "Input Filename:");
     gtk_frame_set_label_align (GTK_FRAME (frame), 0.0, 1.0);
     table=gtk_table_new(3,5,FALSE);
     gtk_container_set_border_width (GTK_CONTAINER (table), 5);
     gtk_container_add(GTK_CONTAINER(frame),table);
     label=gtk_label_new("vel:");
     gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,0,1);
     _vel=gtk_entry_new_with_max_length(100);
     gtk_entry_set_text (GTK_ENTRY (_vel), "vel.dat");
     gtk_table_attach_defaults(GTK_TABLE(table),_vel,1,5,0,1);
     button=gtk_button_new();
     box = xpm_label_box ("pic/seek.jpg", NULL);
     gtk_container_add (GTK_CONTAINER (button), box);
     g_signal_connect(G_OBJECT(button),"clicked",G_CALLBACK(select_file),_vel);
     gtk_table_attach_defaults(GTK_TABLE(table),button,5,6,0,1);
     label=gtk_label_new("epsilu:");
     gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,1,2);
     _epsilu=gtk_entry_new_with_max_length(100);gtk_entry_set_text(GTK_ENTRY (_epsilu), "epsilu.dat");
     gtk_table_attach_defaults(GTK_TABLE(table),_epsilu,1,5,1,2);
     button=gtk_button_new_with_label("...");
     button=gtk_button_new();
     box = xpm_label_box ("pic/seek.jpg", NULL);
     gtk_container_add (GTK_CONTAINER (button), box);
     g_signal_connect(G_OBJECT(button),"clicked",G_CALLBACK(select_file),_epsilu);
     gtk_table_attach_defaults(GTK_TABLE(table),button,5,6,1,2);
     label=gtk_label_new("deta:");gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,2,3);
     _deta=gtk_entry_new_with_max_length(100);
     gtk_entry_set_text (GTK_ENTRY (_deta), "deta.dat");
     gtk_table_attach_defaults(GTK_TABLE(table),_deta,1,5,2,3);
     button=gtk_button_new_with_label("...");
     button=gtk_button_new();
     box = xpm_label_box ("pic/seek.jpg", NULL);
     gtk_container_add (GTK_CONTAINER (button), box);
     g_signal_connect(G_OBJECT(button),"clicked",G_CALLBACK(select_file),_deta);
     gtk_table_attach_defaults(GTK_TABLE(table),button,5,6,2,3);
     boxH = gtk_hbox_new (FALSE, 0);
     gtk_box_pack_start (GTK_BOX (boxV), boxH, TRUE, TRUE, 0);
     frame = gtk_frame_new (NULL);
     gtk_container_set_border_width (GTK_CONTAINER (frame), 10);
     gtk_box_pack_start (GTK_BOX (boxH), frame, TRUE, TRUE, 0);
     gtk_frame_set_label (GTK_FRAME (frame), "Model Parameters:");
     gtk_frame_set_label_align (GTK_FRAME (frame), 0.0, 1.0);
     table=gtk_table_new(3,4,FALSE);
     gtk_container_set_border_width (GTK_CONTAINER (table), 10);
     gtk_container_add(GTK_CONTAINER(frame),table);
     label=gtk_label_new("nx=");
     gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,0,1);
     _nx=gtk_entry_new_with_max_length(10); gtk_entry_set_text(GTK_ENTRY(_nx), "100");
     gtk_table_attach_defaults(GTK_TABLE(table),_nx,1,2,0,1);
     label=gtk_label_new("dx=");
     gtk_table_attach_defaults(GTK_TABLE(table),label,2,3,0,1);
     _dx=gtk_entry_new_with_max_length(10);
     gtk_entry_set_text (GTK_ENTRY (_dx), "5");
     gtk_table_attach_defaults(GTK_TABLE(table),_dx,3,4,0,1);
     label=gtk_label_new("npd=");
     gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,2,3);
     _npd=gtk_entry_new_with_max_length(10);
     gtk_entry_set_text (GTK_ENTRY (_npd), "20");
     gtk_table_attach_defaults(GTK_TABLE(table),_npd,1,2,2,3);
     label=gtk_label_new("nz=");
     gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,1,2);
     _nz=gtk_entry_new_with_max_length(10); gtk_entry_set_text(GTK_ENTRY(_nz), "100");
     gtk_table_attach_defaults(GTK_TABLE(table),_nz,1,2,1,2);
     label=gtk_label_new("dz=");
     gtk_table_attach_defaults(GTK_TABLE(table),label,2,3,1,2);
     _dz=gtk_entry_new_with_max_length(10); gtk_entry_set_text(GTK_ENTRY(_dz), "5");
     gtk_table_attach_defaults(GTK_TABLE(table),_dz,3,4,1,2);
     label=gtk_label_new("mm=");
     gtk_table_attach_defaults(GTK_TABLE(table),label,2,3,2,3);
     _mm=gtk_entry_new_with_max_length(10); gtk_entry_set_text(GTK_ENTRY(_mm), "4");
     gtk_table_attach_defaults(GTK_TABLE(table),_mm,3,4,2,3);
     frame = gtk_frame_new (NULL);
     gtk_container_set_border_width (GTK_CONTAINER(frame),10);
     gtk_box_pack_start (GTK_BOX (boxH), frame, TRUE,TRUE,0);
     gtk_frame_set_label (GTK_FRAME (frame),"ShotParameters:");
     gtk_frame_set_label_align (GTK_FRAME (frame), 0.0, 1.0);
     table=gtk_table_new(3,4,FALSE);
     gtk_container_set_border_width (GTK_CONTAINER (table), 10);
     gtk_container_add(GTK_CONTAINER(frame),table);
     label=gtk_label_new("ns=");
     gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,0,1);
     _ns=gtk_entry_new_with_max_length(10); gtk_entry_set_text(GTK_ENTRY(_ns), "1");
     gtk_table_attach_defaults(GTK_TABLE(table),_ns,1,2,0,1);
     label=gtk_label_new("favg=");
     gtk_table_attach_defaults(GTK_TABLE(table),label,2,3,0,1);
     _favg=gtk_entry_new_with_max_length(10);gtk_entry_set_text(GTK_ENTRY (_favg), "20");
     gtk_table_attach_defaults(GTK_TABLE(table),_favg,3,4,0,1);
     label=gtk_label_new("fs=");
     gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,1,2);
     _fs=gtk_entry_new_with_max_length(10);gtk_entry_set_text (GTK_ENTRY(_fs), "1");
     gtk_table_attach_defaults(GTK_TABLE(table),_fs,1,2,1,2);
     label=gtk_label_new("ds=");
     gtk_table_attach_defaults(GTK_TABLE(table),label,2,3,1,2);
     _ds=gtk_entry_new_with_max_length(10); gtk_entry_set_text(GTK_ENTRY(_ds), "0");
     gtk_table_attach_defaults(GTK_TABLE(table),_ds,3,4,1,2);
     label=gtk_label_new("zs=");
     gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,2,3);
     _zs=gtk_entry_new_with_max_length(10); gtk_entry_set_text(GTK_ENTRY(_zs), "1");
     gtk_table_attach_defaults(GTK_TABLE(table),_zs,1,2,2,3);
     frame = gtk_frame_new (NULL);
     gtk_container_set_border_width (GTK_CONTAINER (frame), 10);
     gtk_box_pack_start (GTK_BOX (boxV), frame, TRUE, TRUE, 0);
     gtk_frame_set_label (GTK_FRAME (frame), "Output File:");
     gtk_frame_set_label_align (GTK_FRAME (frame), 0.0, 1.0);
     table=gtk_table_new(2,6,FALSE);
     gtk_container_set_border_width (GTK_CONTAINER (table), 5);
     gtk_container_add(GTK_CONTAINER(frame),table);
     label=gtk_label_new("shot=");
     gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,0,1);
     _shot=gtk_entry_new_with_max_length(10);gtk_entry_set_text(GTK_ENTRY (_shot), "shot.dat");
     gtk_table_attach_defaults(GTK_TABLE(table),_shot,1,4,0,1);
     label=gtk_label_new("snap=");
     gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,1,2);
     _snap=gtk_entry_new_with_max_length(10);
     gtk_entry_set_text (GTK_ENTRY (_snap), "snap.dat");
     gtk_table_attach_defaults(GTK_TABLE(table),_snap,1,4,1,2);
     button=gtk_button_new_with_label("RUN");
     table=gtk_table_new(3,6,TRUE);
     gtk_container_set_border_width (GTK_CONTAINER (table), 5);
     gtk_box_pack_start (GTK_BOX (boxV), table, TRUE, TRUE, 0);
     g_signal_connect(G_OBJECT(button),"clicked",G_CALLBACK(write_txt),NULL);
     gtk_table_attach_defaults(GTK_TABLE(table),button,2,4,0,1);
     label=gtk_label_new("Please ensure all parameters inputaccuratlly!");
     gtk_table_attach_defaults(GTK_TABLE(table),label,0,6,1,2);
     label=gtk_label_new("Rong Tao !");
     gtk_table_attach_defaults(GTK_TABLE(table),label,0,6,2,3);
     gtk_widget_show_all(window);
}
//###################################Main function#########################################
//###################################Main function#########################################
//###################################Main function#########################################
int main( int argc, char *argv[] )
{
     GtkWidget *window;
     GtkWidget *button;
     GtkWidget *label;
     GtkWidget *table;
     GtkWidget *frame;
     GtkWidget *boxV,*boxH,*box;


     gtk_init(&argc,&argv);
     window=gtk_window_new(GTK_WINDOW_TOPLEVEL);
     gtk_window_set_title(GTK_WINDOW(window),"RT's FD model");
     gtk_window_set_position(GTK_WINDOW(window),GTK_WIN_POS_CENTER);
     g_signal_connect(GTK_OBJECT(window),"delete_event",GTK_SIGNAL_FUNC(gtk_main_quit),NULL);

     boxV = gtk_vbox_new (FALSE, 0);
     gtk_container_add (GTK_CONTAINER (window), boxV);
     gtk_widget_show (boxV);
     frame = gtk_frame_new (NULL);
     gtk_container_set_border_width (GTK_CONTAINER (frame), 10);
     gtk_box_pack_start (GTK_BOX (boxV), frame, TRUE, TRUE, 0);
     gtk_frame_set_label (GTK_FRAME (frame), "FD kind of model:");
     gtk_frame_set_label_align (GTK_FRAME (frame), 0.0, 1.0);
     table=gtk_table_new(1,5,FALSE);
     gtk_container_set_border_width (GTK_CONTAINER (table), 5);
     gtk_container_add(GTK_CONTAINER(frame),table);

     button=gtk_button_new_with_label("VTI model");
     gtk_table_attach_defaults(GTK_TABLE(table),button,0,1,0,1);
     g_signal_connect(G_OBJECT(button),"clicked",G_CALLBACK(VTI_model),NULL);

     button=gtk_button_new_with_label("VTI RTM");
     gtk_table_attach_defaults(GTK_TABLE(table),button,0,1,1,2);
     g_signal_connect(G_OBJECT(button),"clicked",G_CALLBACK(VTI_model),NULL);


     gtk_widget_show_all(window);
     gtk_main ();
     return 0;
}
