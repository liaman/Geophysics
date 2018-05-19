//a############################################
//a##
//a## GUI of cpu vti 2dfd
//a##
//a##           codeing by Rong Tao
//a############################################

#include<gtk/gtk.h>
//#include"CPU_vti2dfd_kernels.c"
//###################################parameter of model#########################################
GtkWidget *file,*_vel,*_deta,*_epsilu;
GtkWidget *_nx, *_nz, *_dx, *_dz, *_npd, *_SV, *_nt, *_dt;
GtkWidget *_favg, *_fs, *_ds, *_ns, *_zs;
GtkWidget *_shot, *_snap;
gint nx, nz, npd, SV, fs, ds, ns, zs, dx, dz, favg, nt, dt;

//###################################connect to model#########################################
void GTK2CPU_vti2dfd(GtkWidget *widget,gpointer data)
{




     const char*str_nx=gtk_entry_get_text(GTK_ENTRY((GtkWidget*)_nx));
     const char*str_nz=gtk_entry_get_text(GTK_ENTRY((GtkWidget*)_nz));
     const char*str_dx=gtk_entry_get_text(GTK_ENTRY((GtkWidget*)_dx));
     const char*str_dz=gtk_entry_get_text(GTK_ENTRY((GtkWidget*)_dz));
     const char*str_npd=gtk_entry_get_text(GTK_ENTRY((GtkWidget*)_npd));
     const char*str_SV=gtk_entry_get_text(GTK_ENTRY((GtkWidget*)_SV));
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
     const char*str_nt=gtk_entry_get_text(GTK_ENTRY((GtkWidget*)_nt));
     const char*str_dt=gtk_entry_get_text(GTK_ENTRY((GtkWidget*)_dt));
     nx=strtod(str_nx,NULL);
     nz=strtod(str_nz,NULL);
     dx=strtod(str_dx,NULL);
     dz=strtod(str_dz,NULL); 
     npd=strtod(str_npd,NULL);
     SV=strtod(str_SV,NULL);
     favg=strtod(str_favg,NULL);
     fs=strtod(str_fs,NULL);
     ns=strtod(str_ns,NULL);
     ds=strtod(str_ds,NULL);
     zs=strtod(str_zs,NULL);
     nt=strtod(str_nt,NULL);
     dt=strtod(str_dt,NULL);

     /* vti2dfd CPU function */
     CPU_vti2dfd(nx,nz,dx,dz,npd,SV,str_vel,str_epsilu,str_deta,favg,ns,fs,ds,zs,str_shot,str_snap,nt,dt);

     gtk_widget_destroy(data);
}


//################################### button run or quit #########################################
void RunOrQuitCpu()
{
     GtkWidget *window;
     GtkWidget *button;
     GtkWidget *hbox;
     GtkWidget *vbox;
     GtkWidget *frame;
     GtkWidget *label;
     GtkWidget *align; 

     window=gtk_window_new(GTK_WINDOW_TOPLEVEL);
     gtk_window_set_title(GTK_WINDOW(window),"Run Window");
     //gtk_window_set_position(GTK_WINDOW(window),GTK_WIN_POS_CENTER);
     g_signal_connect(GTK_OBJECT(window),"delete_event",GTK_SIGNAL_FUNC(gtk_widget_destroy),NULL);

     vbox = gtk_vbox_new (FALSE, 1);
     gtk_container_add(GTK_CONTAINER(window),vbox);

     frame = gtk_frame_new ("Are you sure run this code ?");
    // gtk_widget_set_size_request (frame, 550, 250);
     gtk_box_pack_start (GTK_BOX (vbox), frame, TRUE, TRUE, 5);
     label=gtk_label_new("\n"
                         "   Please make sure the parameters\n"
                         "you input is accurate, and select Run\n"
                         "or Quit!\n"    );

     gtk_container_add(GTK_CONTAINER(frame),label);

     hbox = gtk_hbox_new (FALSE, 1);
     gtk_box_pack_start (GTK_BOX (vbox), hbox, TRUE, TRUE, 0);

      button = gtk_button_new_with_label ("Run");
      g_signal_connect (G_OBJECT (button), "clicked", G_CALLBACK (GTK2CPU_vti2dfd), window);
      align = gtk_alignment_new(0.5,0.5,0.5,0.5);  
      gtk_container_add(GTK_CONTAINER(align),button);  
      gtk_box_pack_start (GTK_BOX(hbox), align, TRUE, TRUE, 0);

      button = gtk_button_new_with_label ("Quit");
      g_signal_connect (G_OBJECT (button), "clicked", G_CALLBACK (window_destroy), window);
      align = gtk_alignment_new(0.5,0.5,0.5,0.5);  
      gtk_container_add(GTK_CONTAINER(align),button);  
      gtk_box_pack_start (GTK_BOX(hbox), align, TRUE, TRUE, 0);

      gtk_widget_show_all (window);

}
//###################################VTI model#########################################
void GTK_cpuvti2dfd(GtkWidget *widget,gpointer data)
{
     GtkWidget *window;
     GtkWidget *button;
     GtkWidget *label;
     GtkWidget *table;
     GtkWidget *frame;
     GtkWidget *boxV,*boxH,*box;
     GtkWidget *align; 


     window=gtk_window_new(GTK_WINDOW_TOPLEVEL);
     gtk_window_set_title(GTK_WINDOW(window),"2D CPU-VTI");
     //gtk_window_set_position(GTK_WINDOW(window),GTK_WIN_POS_CENTER);
     g_signal_connect(GTK_OBJECT(window),"delete_event",GTK_SIGNAL_FUNC(gtk_widget_destroy),NULL);

////////////////////////////////////a////////////////////////////////////////////
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
     label=gtk_label_new("vel[m/s]:");
     gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,0,1);
     _vel=gtk_entry_new_with_max_length(100);
     gtk_entry_set_text (GTK_ENTRY (_vel), "vel.dat");
     gtk_table_attach_defaults(GTK_TABLE(table),_vel,1,5,0,1);
     button=gtk_button_new_with_label("...");
     align = gtk_alignment_new(0,0.5,0,0);  
     gtk_container_add(GTK_CONTAINER(align),button);  
     g_signal_connect(G_OBJECT(button),"clicked",G_CALLBACK(select_file),_vel);
     gtk_table_attach_defaults(GTK_TABLE(table),align,5,6,0,1);
     label=gtk_label_new("epsilon:");
     gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,1,2);
     _epsilu=gtk_entry_new_with_max_length(100);gtk_entry_set_text(GTK_ENTRY (_epsilu), "epsilon.dat");
     gtk_table_attach_defaults(GTK_TABLE(table),_epsilu,1,5,1,2);
     button=gtk_button_new_with_label("...");
     align = gtk_alignment_new(0,0.5,0,0);  
     gtk_container_add(GTK_CONTAINER(align),button);  
     g_signal_connect(G_OBJECT(button),"clicked",G_CALLBACK(select_file),_epsilu);
     gtk_table_attach_defaults(GTK_TABLE(table),align,5,6,1,2);
     label=gtk_label_new("delta:");gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,2,3);
     _deta=gtk_entry_new_with_max_length(100);
     gtk_entry_set_text (GTK_ENTRY (_deta), "delta.dat");
     gtk_table_attach_defaults(GTK_TABLE(table),_deta,1,5,2,3);
     button=gtk_button_new_with_label("...");
     align = gtk_alignment_new(0,0.5,0,0);  
     gtk_container_add(GTK_CONTAINER(align),button);  
     g_signal_connect(G_OBJECT(button),"clicked",G_CALLBACK(select_file),_deta);
     gtk_table_attach_defaults(GTK_TABLE(table),align,5,6,2,3);

////////////////////////////////////a////////////////////////////////////////////
     boxH = gtk_hbox_new (FALSE, 0);
     gtk_box_pack_start (GTK_BOX (boxV), boxH, TRUE, TRUE, 0);

     frame = gtk_frame_new (NULL);
     gtk_container_set_border_width (GTK_CONTAINER (frame), 10);
     gtk_box_pack_start (GTK_BOX (boxH), frame, TRUE, TRUE, 0);

     gtk_frame_set_label (GTK_FRAME (frame), "Model Parameters:");
     gtk_frame_set_label_align (GTK_FRAME (frame), 0.0, 1.0);

     table=gtk_table_new(6,2,FALSE);
     gtk_container_set_border_width (GTK_CONTAINER (table), 10);
     gtk_container_add(GTK_CONTAINER(frame),table);

     label=gtk_label_new("nx[grid]=");
     gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,0,1);
     _nx=gtk_entry_new_with_max_length(10); 
     gtk_entry_set_text(GTK_ENTRY(_nx), "100");
     gtk_table_attach_defaults(GTK_TABLE(table),_nx,1,2,0,1);

     label=gtk_label_new("dx[m]=");
     gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,1,2);
     _dx=gtk_entry_new_with_max_length(10);
     gtk_entry_set_text (GTK_ENTRY (_dx), "5");
     gtk_table_attach_defaults(GTK_TABLE(table),_dx,1,2,1,2);

     label=gtk_label_new("nz[grid]=");
     gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,2,3);
     _nz=gtk_entry_new_with_max_length(10); 
     gtk_entry_set_text(GTK_ENTRY(_nz), "100");
     gtk_table_attach_defaults(GTK_TABLE(table),_nz,1,2,2,3);

     label=gtk_label_new("dz[m]=");
     gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,3,4);
     _dz=gtk_entry_new_with_max_length(10); 
     gtk_entry_set_text(GTK_ENTRY(_dz), "5");
     gtk_table_attach_defaults(GTK_TABLE(table),_dz,1,2,3,4);


     label=gtk_label_new("npd[grid]=");
     gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,4,5);
     _npd=gtk_entry_new_with_max_length(10);
     gtk_entry_set_text (GTK_ENTRY (_npd), "20");
     gtk_table_attach_defaults(GTK_TABLE(table),_npd,1,2,4,5);

     label=gtk_label_new("SV[0,1]=");
     gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,5,6);
     _SV=gtk_entry_new_with_max_length(10); 
     gtk_entry_set_text(GTK_ENTRY(_SV), "1");
     gtk_table_attach_defaults(GTK_TABLE(table),_SV,1,2,5,6);


////////////////////////////////////a////////////////////////////////////////////
     frame = gtk_frame_new (NULL);
     gtk_container_set_border_width (GTK_CONTAINER(frame),10);
     gtk_box_pack_start (GTK_BOX (boxH), frame, TRUE,TRUE,0);

     gtk_frame_set_label (GTK_FRAME (frame),"ShotParameters:");
     gtk_frame_set_label_align (GTK_FRAME (frame), 0.0, 1.0);

     table=gtk_table_new(6,2,FALSE);
     gtk_container_set_border_width (GTK_CONTAINER (table), 10);
     gtk_container_add(GTK_CONTAINER(frame),table);

     label=gtk_label_new("fm[Hz]=");
     gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,0,1);
     _favg=gtk_entry_new_with_max_length(10);
     gtk_entry_set_text(GTK_ENTRY (_favg), "20");
     gtk_table_attach_defaults(GTK_TABLE(table),_favg,1,2,0,1);

     label=gtk_label_new("ns=");
     gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,1,2);
     _ns=gtk_entry_new_with_max_length(10); 
     gtk_entry_set_text(GTK_ENTRY(_ns), "1");
     gtk_table_attach_defaults(GTK_TABLE(table),_ns,1,2,1,2);

     label=gtk_label_new("fs[grid]=");
     gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,2,3);
     _fs=gtk_entry_new_with_max_length(10);
     gtk_entry_set_text (GTK_ENTRY(_fs), "1");
     gtk_table_attach_defaults(GTK_TABLE(table),_fs,1,2,2,3);

     label=gtk_label_new("ds[grid]=");
     gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,3,4);
     _ds=gtk_entry_new_with_max_length(10); 
     gtk_entry_set_text(GTK_ENTRY(_ds), "0");
     gtk_table_attach_defaults(GTK_TABLE(table),_ds,1,2,3,4);

     label=gtk_label_new("zs[grid]=");
     gtk_table_attach_defaults(GTK_TABLE(table),label,0,1,4,5);
     _zs=gtk_entry_new_with_max_length(10); 
     gtk_entry_set_text(GTK_ENTRY(_zs), "1");
     gtk_table_attach_defaults(GTK_TABLE(table),_zs,1,2,4,5);


////////////////////////////////////a////////////////////////////////////////////
     boxH = gtk_hbox_new (FALSE, 0);
     gtk_box_pack_start (GTK_BOX (boxV), boxH, TRUE, TRUE, 0);

     frame = gtk_frame_new (NULL);
     gtk_container_set_border_width (GTK_CONTAINER (frame), 10);
     gtk_box_pack_start (GTK_BOX (boxH), frame, TRUE, TRUE, 0);

     gtk_frame_set_label (GTK_FRAME (frame), "Time Sampling:");
     gtk_frame_set_label_align (GTK_FRAME (frame), 0.0, 1.0);
     table=gtk_table_new(2,2,FALSE);
     gtk_container_set_border_width (GTK_CONTAINER (table), 10);
     gtk_container_add(GTK_CONTAINER(frame),table);
     label=gtk_label_new("nt=");
     align = gtk_alignment_new(0,0.5,0,0);  
     gtk_container_add(GTK_CONTAINER(align),label);  
     gtk_table_attach_defaults(GTK_TABLE(table),align,0,1,0,1);
     _nt=gtk_entry_new_with_max_length(10); 
     gtk_entry_set_text(GTK_ENTRY(_nt), "1001");
     gtk_table_attach_defaults(GTK_TABLE(table),_nt,1,2,0,1);
     label=gtk_label_new("dt(us)=");
     align = gtk_alignment_new(0,0.5,0,0);  
     gtk_container_add(GTK_CONTAINER(align),label);  
     gtk_table_attach_defaults(GTK_TABLE(table),align,0,1,1,2);
     _dt=gtk_entry_new_with_max_length(10);
     gtk_entry_set_text(GTK_ENTRY (_dt), "500");
     gtk_table_attach_defaults(GTK_TABLE(table),_dt,1,2,1,2);

////////////////////////////////////a////////////////////////////////////////////
     frame = gtk_frame_new (NULL);
     gtk_container_set_border_width (GTK_CONTAINER (frame), 10);
     gtk_box_pack_start (GTK_BOX (boxH), frame, TRUE, TRUE, 0);

     gtk_frame_set_label (GTK_FRAME (frame), "Output Filename:");
     gtk_frame_set_label_align (GTK_FRAME (frame), 0.0, 1.0);
     table=gtk_table_new(2,6,FALSE);
     gtk_container_set_border_width (GTK_CONTAINER (table), 5);
     gtk_container_add(GTK_CONTAINER(frame),table);
     label=gtk_label_new("shot=");
     align = gtk_alignment_new(0,0.5,0,0);  
     gtk_container_add(GTK_CONTAINER(align),label);  
     gtk_table_attach_defaults(GTK_TABLE(table),align,0,1,0,1);
     _shot=gtk_entry_new_with_max_length(50);gtk_entry_set_text(GTK_ENTRY (_shot), "shot.dat");
     gtk_table_attach_defaults(GTK_TABLE(table),_shot,1,4,0,1);
     label=gtk_label_new("snap=");
     align = gtk_alignment_new(0,0.5,0,0);  
     gtk_container_add(GTK_CONTAINER(align),label);  
     gtk_table_attach_defaults(GTK_TABLE(table),align,0,1,1,2);
     _snap=gtk_entry_new_with_max_length(50);
     gtk_entry_set_text (GTK_ENTRY (_snap), "snap.dat");
     gtk_table_attach_defaults(GTK_TABLE(table),_snap,1,4,1,2);


////////////////////////////////////a////////////////////////////////////////////
     button=gtk_button_new_with_label("Sure and Run");
     align = gtk_alignment_new(0,0.5,0,0);  
     gtk_container_add(GTK_CONTAINER(align),button);  

     table=gtk_table_new(2,6,TRUE);
     gtk_container_set_border_width (GTK_CONTAINER (table), 5);
     gtk_box_pack_start (GTK_BOX (boxV), table, TRUE, TRUE, 0);
     g_signal_connect(G_OBJECT(button),"clicked",G_CALLBACK(RunOrQuitCpu),NULL);
     gtk_table_attach_defaults(GTK_TABLE(table),align,5,6,1,2);
     label=gtk_label_new("Please ensure all parameters input accuratelly!");
     gtk_table_attach_defaults(GTK_TABLE(table),label,0,6,0,1);

     gtk_widget_show_all(window);
}

