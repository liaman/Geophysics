#include<gtk/gtk.h>
#include<stdio.h>
#include<malloc.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include"include/function/GTK_menu.c"
#include"GTK_gpuvti2dfd.c"
#include"GTK_gpuvti3dfd.c"
#include"GTK_cpuvti2dfd.c"

//###################################Main function#########################################
int main( int argc, char *argv[] )
{
     GtkWidget *window;
     GtkWidget *button;
     GtkWidget *label;
     GtkWidget *table;
     GtkWidget *frame;
     GtkWidget *boxV,*boxH,*boxbutton,*boxH2;
     GtkWidget *align; 


     gtk_init(&argc,&argv);
     window=gtk_window_new(GTK_WINDOW_TOPLEVEL);
     gtk_window_set_title(GTK_WINDOW(window),"VTI-FD-GTK");
     gtk_window_set_position(GTK_WINDOW(window),GTK_WIN_POS_CENTER);
     g_signal_connect(GTK_OBJECT(window),"delete_event",GTK_SIGNAL_FUNC(gtk_main_quit),NULL);

     boxV = gtk_vbox_new (FALSE, 0);
     gtk_container_add (GTK_CONTAINER (window), boxV);
     gtk_menu(window,boxV);

     boxH2 = gtk_hbox_new (FALSE, 0);
     gtk_container_add (GTK_CONTAINER (boxV), boxH2);
     gtk_widget_set_size_request (boxH2, 500, 100);

/*********************************frame**************************************/
          /****** CPU ******/
          frame = gtk_frame_new (NULL);
          gtk_container_set_border_width (GTK_CONTAINER (frame), 10);
          align = gtk_alignment_new(0.5,0,0,0);  
          gtk_container_add(GTK_CONTAINER(align),frame); 
          gtk_box_pack_start (GTK_BOX (boxH2), align, TRUE, TRUE, 0);
          gtk_frame_set_label (GTK_FRAME (frame), "CPU");
          gtk_frame_set_label_align (GTK_FRAME (frame), 0.0, 1.0);


               boxH = gtk_hbox_new (TRUE, 0);
               gtk_container_add (GTK_CONTAINER (frame), boxH);
               gtk_container_set_border_width (GTK_CONTAINER (boxH), 10);

                    boxbutton=xpm_label_box("include/picture/2d.png","2D");
                    button=gtk_button_new();
                    gtk_container_add(GTK_CONTAINER(button),boxbutton); 
                    align = gtk_alignment_new(0.5,0.5,0,0);  
                    gtk_container_add(GTK_CONTAINER(align),button);  
                    gtk_box_pack_start (GTK_BOX (boxH), align, FALSE, FALSE, 0);
                    g_signal_connect(G_OBJECT(button),"clicked",G_CALLBACK(GTK_cpuvti2dfd),NULL);

/*********************************frame**************************************/
          /****** GPU ******/
          frame = gtk_frame_new (NULL);
          gtk_container_set_border_width (GTK_CONTAINER (frame), 10);
          align = gtk_alignment_new(0.5,0,0,0);  
          gtk_container_add(GTK_CONTAINER(align),frame); 
          gtk_box_pack_start (GTK_BOX (boxH2), align, TRUE, TRUE, 0);
          gtk_frame_set_label (GTK_FRAME (frame), "CPU&GPU");
          gtk_frame_set_label_align (GTK_FRAME (frame), 0.0, 1.0);


               boxH = gtk_hbox_new (TRUE, 0);
               gtk_container_add (GTK_CONTAINER (frame), boxH);
               gtk_container_set_border_width (GTK_CONTAINER (boxH), 10);

                    boxbutton=xpm_label_box("include/picture/2d.png","2D");
                    button=gtk_button_new();
                    gtk_container_add(GTK_CONTAINER(button),boxbutton); 
                    align = gtk_alignment_new(0.5,0.5,0,0);  
                    gtk_container_add(GTK_CONTAINER(align),button);  
                    gtk_box_pack_start (GTK_BOX (boxH), align, FALSE, FALSE, 0);
                    g_signal_connect(G_OBJECT(button),"clicked",G_CALLBACK(GTK_gpuvti2dfd),NULL);


                    boxbutton=xpm_label_box("include/picture/3d.png","3D");
                    button=gtk_button_new();
                    gtk_container_add(GTK_CONTAINER(button),boxbutton); 
                    align = gtk_alignment_new(0.5,0.5,0,0);  
                    gtk_container_add(GTK_CONTAINER(align),button); 
                    gtk_box_pack_start (GTK_BOX (boxH), align, FALSE, FALSE, 0);
                    g_signal_connect(G_OBJECT(button),"clicked",G_CALLBACK(GTK_gpuvti3dfd),NULL);

/*********************************frame**************************************/
          /****** CPU&GPU ******/
          frame = gtk_frame_new (NULL);
          gtk_container_set_border_width (GTK_CONTAINER (frame), 10);
          align = gtk_alignment_new(0.5,0,0,0);  
          gtk_container_add(GTK_CONTAINER(align),frame); 
          gtk_box_pack_start (GTK_BOX (boxH2), align, TRUE, TRUE, 0);
          gtk_frame_set_label (GTK_FRAME (frame), "CPU&GPU&MPI");
          gtk_frame_set_label_align (GTK_FRAME (frame), 0.0, 1.0);


               boxH = gtk_hbox_new (TRUE, 0);
               gtk_container_add (GTK_CONTAINER (frame), boxH);
               gtk_container_set_border_width (GTK_CONTAINER (boxH), 10);

                    boxbutton=xpm_label_box("include/picture/2d.png","2D");
                    button=gtk_button_new();
                    gtk_container_add(GTK_CONTAINER(button),boxbutton); 
                    align = gtk_alignment_new(0.5,0.5,0,0);  
                    gtk_container_add(GTK_CONTAINER(align),button);  
                    gtk_box_pack_start (GTK_BOX (boxH), align, FALSE, FALSE, 0);
                    //g_signal_connect(G_OBJECT(button),"clicked",G_CALLBACK(GTK_vti2dfd),NULL);


                    boxbutton=xpm_label_box("include/picture/3d.png","3D");
                    button=gtk_button_new();
                    gtk_container_add(GTK_CONTAINER(button),boxbutton); 
                    align = gtk_alignment_new(0.5,0.5,0,0);  
                    gtk_container_add(GTK_CONTAINER(align),button); 
                    gtk_box_pack_start (GTK_BOX (boxH), align, FALSE, FALSE, 0);
                   // g_signal_connect(G_OBJECT(button),"clicked",G_CALLBACK(GTK_vti2dfd),NULL);


          frame = gtk_frame_new (NULL);
          gtk_widget_set_size_request (frame, 450, 370);
          boxbutton=xpm_label_box("include/picture/background3.png"," ");
          gtk_container_add (GTK_CONTAINER (frame), boxbutton);
          gtk_box_pack_start (GTK_BOX (boxV), frame, FALSE, FALSE, 0);
/*********************************frame**************************************/
          frame = gtk_frame_new (NULL);
          gtk_container_set_border_width (GTK_CONTAINER (frame), 10);
          gtk_box_pack_start (GTK_BOX (boxV), frame, FALSE, FALSE, 0);

               label=gtk_label_new("China University of Petroleum (East China)\n"
                                   "                                   Developer:Rong Tao");
               gtk_container_add(GTK_CONTAINER(frame),label);

     
     gtk_widget_show_all(window);
     gtk_main ();
     return 0;
}
