//a############################################
//a##
//a## GUI  Menu
//a##
//a##           codeing by Rong Tao
//a############################################
/* menu */
#include <gtk/gtk.h>
void on_menu_activate (GtkMenuItem* item,gpointer data)
{
   g_print("Menuitem %s is pressed.\n",(gchar*)data);
}
void on_button_clicked (GtkButton* button,gpointer data)
{
   g_print("Toolbaritem %s is pressed.\n",(gchar*)data);
}
void About_window()
{
     GtkWidget *window;
     GtkWidget *frame;
     GtkWidget *label;

     window=gtk_window_new(GTK_WINDOW_TOPLEVEL);
     gtk_window_set_title(GTK_WINDOW(window),"About");
     gtk_window_set_position(GTK_WINDOW(window),GTK_WIN_POS_CENTER);
     g_signal_connect(GTK_OBJECT(window),"delete_event",GTK_SIGNAL_FUNC(gtk_widget_destroy),NULL);

     frame = gtk_frame_new ("software introduction");
     gtk_widget_set_size_request (frame, 550, 250);
     gtk_container_add(GTK_CONTAINER(window),frame);

     label=gtk_label_new("**************************************************************************************\n"
                         "*\n"
                         "*          VTI media Finite difference forward software !\n"
                         "*\n"
                         "*   The software includes two-dimensional and three-dimensional\n"
                         "* VTI media GPU acceleration algorithm and CPU and GPU co-operation \n"
                         "* seismic finite difference forward modeling !\n"
                         "*\n"
                         "*\n"
                         "*                                            Developer: Rong Tao\n"
                         "*\n"
                         "**************************************************************************************\n");
   //  gtk_label_set_markup(GTK_LABEL(label),  
   //            "<span foreground='red' underline='double' underline_color='blue' font_desc='32'>test label!</span>"); 
     gtk_container_add(GTK_CONTAINER(frame),label);

     gtk_widget_show_all(window);
}

void gtk_menu(GtkWidget* window,GtkWidget* boxV)
{

   GtkWidget* box;
   GtkWidget* menubar;
   GtkWidget* menu;
   GtkWidget* editmenu;
   GtkWidget* helpmenu;
   GtkWidget* rootmenu;
   GtkWidget* menuitem;
   GtkAccelGroup* accel_group ;
   GtkWidget *align; 
   GtkWidget* toolbar ;
   GtkWidget* entry ;
   GtkWidget *label;

   accel_group = gtk_accel_group_new();

   gtk_window_add_accel_group(GTK_WINDOW(window),accel_group);// AccelGroup

   box = gtk_hbox_new(FALSE,0);
   align = gtk_alignment_new(0,0.5,0,0);  
   gtk_container_add(GTK_CONTAINER(align),box); 
   gtk_box_pack_start(GTK_BOX(boxV),align,FALSE,FALSE,0);


/********************************Menubar**********************************/
   menu = gtk_menu_new();
   menuitem = gtk_image_menu_item_new_from_stock(GTK_STOCK_NEW,accel_group);
   gtk_menu_shell_append(GTK_MENU_SHELL(menu),menuitem);
   g_signal_connect(G_OBJECT(menuitem),"activate",G_CALLBACK(on_menu_activate),(gpointer)("New"));
   menuitem = gtk_image_menu_item_new_from_stock(GTK_STOCK_OPEN,accel_group);
   gtk_menu_shell_append(GTK_MENU_SHELL(menu),menuitem);
   g_signal_connect(G_OBJECT(menuitem),"activate",G_CALLBACK(on_menu_activate),(gpointer)("Open"));
   menuitem = gtk_image_menu_item_new_from_stock(GTK_STOCK_SAVE,accel_group);
   gtk_menu_shell_append(GTK_MENU_SHELL(menu),menuitem);
   g_signal_connect(G_OBJECT(menuitem),"activate",G_CALLBACK(on_menu_activate),(gpointer)("Save"));
   menuitem = gtk_image_menu_item_new_from_stock(GTK_STOCK_SAVE_AS,accel_group);
   gtk_menu_shell_append(GTK_MENU_SHELL(menu),menuitem);
   g_signal_connect(G_OBJECT(menuitem),"activate",G_CALLBACK(on_menu_activate),(gpointer)("Save As"));
   menuitem = gtk_separator_menu_item_new();
   gtk_menu_shell_append(GTK_MENU_SHELL(menu),menuitem);
   menuitem = gtk_image_menu_item_new_from_stock(GTK_STOCK_QUIT,accel_group);
   gtk_menu_shell_append(GTK_MENU_SHELL(menu),menuitem);
   g_signal_connect(G_OBJECT(menuitem),"activate",G_CALLBACK(on_menu_activate),(gpointer)("Exit"));
   rootmenu = gtk_menu_item_new_with_label(" File ");
   gtk_menu_item_set_submenu(GTK_MENU_ITEM(rootmenu),menu);

   menubar = gtk_menu_bar_new();
   gtk_menu_shell_append(GTK_MENU_SHELL(menubar),rootmenu);
   rootmenu = gtk_menu_item_new_with_label(" Edit ");
   editmenu = gtk_menu_new();
   menuitem = gtk_image_menu_item_new_from_stock(GTK_STOCK_CUT,accel_group);
   gtk_menu_shell_append(GTK_MENU_SHELL(editmenu),menuitem);
   g_signal_connect(G_OBJECT(menuitem),"activate",G_CALLBACK(on_menu_activate),(gpointer)("Cut"));
   menuitem = gtk_image_menu_item_new_from_stock(GTK_STOCK_COPY,accel_group);
   gtk_menu_shell_append(GTK_MENU_SHELL(editmenu),menuitem);
   g_signal_connect(G_OBJECT(menuitem),"activate",G_CALLBACK(on_menu_activate),(gpointer)("Copy"));
   menuitem = gtk_image_menu_item_new_from_stock(GTK_STOCK_PASTE,accel_group);
   gtk_menu_shell_append(GTK_MENU_SHELL(editmenu),menuitem);
   g_signal_connect(G_OBJECT(menuitem),"activate",G_CALLBACK(on_menu_activate),(gpointer)("Paste"));
   menuitem = gtk_image_menu_item_new_from_stock(GTK_STOCK_FIND,accel_group);
   gtk_menu_shell_append(GTK_MENU_SHELL(editmenu),menuitem);
   g_signal_connect(G_OBJECT(menuitem),"activate",G_CALLBACK(on_menu_activate),(gpointer)("Search"));
   gtk_menu_item_set_submenu(GTK_MENU_ITEM(rootmenu),editmenu);
   gtk_menu_shell_append(GTK_MENU_SHELL(menubar),rootmenu);
   rootmenu = gtk_menu_item_new_with_label(" Help ");
   helpmenu = gtk_menu_new();
   menuitem = gtk_image_menu_item_new_from_stock(GTK_STOCK_HELP,accel_group);
   gtk_menu_shell_append(GTK_MENU_SHELL(helpmenu),menuitem);
   g_signal_connect(G_OBJECT(menuitem),"activate",G_CALLBACK(on_menu_activate),(gpointer)("Help"));
   menuitem = gtk_menu_item_new_with_label(" About... ");
   gtk_menu_shell_append(GTK_MENU_SHELL(helpmenu),menuitem);
   g_signal_connect(G_OBJECT(menuitem),"activate",G_CALLBACK(About_window),(gpointer)("About"));
   gtk_menu_item_set_submenu(GTK_MENU_ITEM(rootmenu),helpmenu);
   gtk_menu_shell_append(GTK_MENU_SHELL(menubar),rootmenu);



/********************************Toolbar**********************************/
   toolbar = gtk_toolbar_new();
   gtk_toolbar_insert_stock(GTK_TOOLBAR(toolbar),GTK_STOCK_NEW,
                "Build new file","New",GTK_SIGNAL_FUNC(on_button_clicked),(gpointer)("New"),-1);

   gtk_toolbar_insert_stock(GTK_TOOLBAR(toolbar),GTK_STOCK_OPEN,
                "Open file","Open", GTK_SIGNAL_FUNC(on_button_clicked), (gpointer)("Open"),-1);

   gtk_toolbar_insert_stock(GTK_TOOLBAR(toolbar),GTK_STOCK_SAVE,
                "Save file","Save", GTK_SIGNAL_FUNC(on_button_clicked), (gpointer)("Save"),-1);

   gtk_toolbar_append_space(GTK_TOOLBAR(toolbar));
   label = gtk_label_new(" Find：");
   gtk_toolbar_append_widget(GTK_TOOLBAR(toolbar),label, "Label","Label");
   entry = gtk_entry_new();
   gtk_toolbar_append_widget(GTK_TOOLBAR(toolbar),entry, "Entry","Entry");
   gtk_toolbar_append_space(GTK_TOOLBAR(toolbar));
   gtk_toolbar_insert_stock(GTK_TOOLBAR(toolbar),GTK_STOCK_CUT,
               "Cut","Cut", GTK_SIGNAL_FUNC(on_button_clicked), (gpointer)("Cut"),-1);
   gtk_toolbar_insert_stock(GTK_TOOLBAR(toolbar),GTK_STOCK_COPY,
               "Copy","Copy", GTK_SIGNAL_FUNC(on_button_clicked), (gpointer)("Copy"),-1);
   gtk_toolbar_insert_stock(GTK_TOOLBAR(toolbar),GTK_STOCK_PASTE,
               "Paste","Paste", GTK_SIGNAL_FUNC(on_button_clicked), (gpointer)("Paste"),-1);
   gtk_toolbar_append_space(GTK_TOOLBAR(toolbar));
   gtk_toolbar_insert_stock(GTK_TOOLBAR(toolbar),GTK_STOCK_QUIT,
               "Quit","Quit", GTK_SIGNAL_FUNC(on_button_clicked), (gpointer)("Quit"),-1);
   gtk_toolbar_set_style(GTK_TOOLBAR(toolbar),GTK_TOOLBAR_ICONS);
   gtk_toolbar_set_icon_size(GTK_TOOLBAR(toolbar),
      GTK_ICON_SIZE_SMALL_TOOLBAR);



   gtk_box_pack_start(GTK_BOX(box),menubar,FALSE,FALSE,0);
   gtk_box_pack_start(GTK_BOX(box),toolbar,FALSE,FALSE,0);


}
