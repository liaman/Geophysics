import java.io.*;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Insets;
import javax.swing.JPanel;
import javax.swing.JLabel;
import javax.swing.JFrame;
import javax.swing.JButton;
import java.util.Random;
import javax.swing.BorderFactory;
import java.io.*;
import java.awt.BorderLayout;
import java.awt.GridLayout;
import java.awt.GridBagLayout;
import java.awt.Toolkit;
import java.awt.Font;
import javax.swing.border.EmptyBorder;
import java.awt.Insets;
import javax.swing.Box;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.BoxLayout;
import javax.swing.JTextField;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseListener;
import java.awt.event.MouseAdapter;
import java.awt.event.WindowAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.Toolkit;
import java.awt.Dimension;
import java.awt.BorderLayout;
import java.awt.GridLayout;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Color;
import java.awt.FileDialog;
import java.awt.Color;
import java.awt.GraphicsEnvironment;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Insets;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.BasicStroke;
import java.awt.Image;
import java.awt.image.BufferedImage;


import javax.swing.JColorChooser;
import javax.swing.JComponent;
import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextPane;
import javax.swing.JEditorPane;
import javax.swing.JButton;
import javax.swing.JMenu;
import javax.swing.JPopupMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JProgressBar;
import javax.swing.JToolBar;
import javax.swing.JList;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JDialog;
import javax.swing.JComboBox;
import javax.swing.JCheckBox;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.UIManager;
import javax.swing.BorderFactory;
import javax.swing.border.EtchedBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.Timer;
import javax.swing.JOptionPane;
import javax.swing.JFileChooser;
import javax.swing.filechooser.FileFilter;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.text.Style;
import javax.swing.text.StyleConstants;
import javax.swing.text.StyledDocument;
import javax.swing.text.SimpleAttributeSet;
import javax.swing.text.MutableAttributeSet;
import javax.swing.text.AttributeSet;
import javax.swing.border.EmptyBorder;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;

import java.text.DateFormat;

import java.util.Calendar;
import java.util.Date;
import java.util.Locale;
import java.util.Scanner;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;

import java.nio.file.Path;

import javax.imageio.ImageIO;

import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JFileChooser;
import javax.swing.filechooser.FileFilter;
import javax.swing.filechooser.FileNameExtensionFilter;

import java.awt.GridBagConstraints;

import javax.swing.JPopupMenu;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;


import java.awt.FileDialog;


import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionAdapter;
import java.awt.event.MouseWheelListener;
import java.awt.event.MouseWheelEvent;
import java.awt.RenderingHints;
import java.awt.BasicStroke;


import java.awt.Point;
import java.nio.file.Path;

import javax.imageio.ImageIO;

import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;

import java.awt.Image;
import java.awt.image.BufferedImage;
import javax.swing.JOptionPane;



//a#############################################################################################
//a##
//a##                myLaganRayTracingParasJFrame -> the timage parameters basic frame
//a##
//a##                    myLaganRayTracingJFrame -> the frame include a image panel
//a##
//a##                        myLaganRayTracingJPanel -> the panel to paint
//a##
//a##                myLaganRayTracingAboutDialog -> About dialog
//a##
//a############################################################################################# 
//a########################################################
//a##
//a##   Lagan Raytracing
//a##
//a########################################################
class Lagan{

        private float pi = 3.141592653f;
        private int sx, sz, nx, nz, ngrid;
        private float dx,dz,s;
        private float[] v0;
        private int countbndr;

        Lagan(int sx, int sz, int nx, int nz, float dx, float dz, float[]v0, float s){
                this.sx = sx;
                this.sz = sz;
                this.nx = nx;
                this.nz = nz;
                this.dx = dx;
                this.dz = dz;
                this.v0 = v0;
                this.s = s;
                ngrid = nx * nz;
          }
        public float[][] makeRay(float angle)throws java.io.IOException{

                float p_x,p_z,n_x,n_z,l_x,l_z;
                float p_xend,p_zend;
                float n_xnew,n_znew;
                int ip_lux,ip_ldx,ip_rux,ip_rdx;
                int ip_luz,ip_ldz,ip_ruz,ip_rdz;
                float time,v;

                angle = angle*pi/180.0f;

                /** get the raypath length -> countbndr **/
                p_x = sx*dx;
                p_z = sz*dz;
                n_x = (float)Math.cos(angle);
                n_z = (float)Math.sin(angle);
                countbndr = 0;
                do{
                        /* cal_gridpoint */
                        ip_lux = (int)(p_x/dx);
                        ip_luz = (int)(p_z/dz);
                        ip_ldx = ip_lux;
                        ip_ldz = ip_luz+1;
                        ip_rux = ip_lux+1;
                        ip_ruz = ip_luz;
                        ip_rdx = ip_lux+1;
                        ip_rdz = ip_luz+1;

                        /* cal_gridvel */
                        float l_x1 = (v0[ip_rux*nz + ip_ruz]-v0[ip_lux*nz + ip_luz])/dx;
                        float l_x2 = (v0[ip_rdx*nz + ip_rdz]-v0[ip_ldx*nz + ip_ldz])/dx;
                        l_x = ( l_x1+l_x2 ) / 2.0f;
                        float l_z1 = (v0[ip_ldx*nz + ip_ldz]-v0[ip_lux*nz + ip_luz])/dz;
                        float l_z2 = (v0[ip_rdx*nz + ip_rdz]-v0[ip_rux*nz + ip_ruz])/dz;
                        l_z = (l_z1+l_z2)/2.0f;
                        v = v0[ip_rux*nz + ip_ruz]+v0[ip_lux*nz + ip_luz]+v0[ip_rdx*nz + ip_rdz]
                                        +v0[ip_ldx*nz + ip_ldz];
                        v = v/4.0f;

                        /* cal_path */
                        float dotmult_ln = l_x*n_x+l_z*n_z;
                        float dotmult_ll = l_x*l_x+l_z*l_z;
                        p_xend = p_x+n_x*s*(1f+dotmult_ln*0.5f*s/v)-0.5f*l_x*s*s/v
                                        -n_x*s*s*s*(dotmult_ll-dotmult_ln*dotmult_ln)/(6.0f*v*v);
                        p_zend = p_z+n_z*s*(1f+dotmult_ln*0.5f*s/v)-0.5f*l_z*s*s/v
                                        -n_z*s*s*s*(dotmult_ll-dotmult_ln*dotmult_ln)/(6.0f*v*v);
                        n_xnew = n_x*(1f+dotmult_ln*s/v)-l_x*s/v
                                        -n_x*s*s*(dotmult_ll-dotmult_ln*dotmult_ln)/(2.0f*v*v);
                        n_znew = n_z*(1f+dotmult_ln*s/v)-l_z*s/v
                                        -n_z*s*s*(dotmult_ll-dotmult_ln*dotmult_ln)/(2.0f*v*v);

                        /* buffer */
                        p_x = p_xend;
                        p_z = p_zend;
                        n_x = n_xnew;
                        n_z = n_znew;
                        
                        countbndr++;
                }while((p_xend>=0.0f)&&(p_xend<((nx-1f)*dx))&&(p_zend>=0.0f)&&(p_zend<(nz-1f)*dz));

                /** make raytracing -> save to raypath[2][..] **/
                p_x = sx*dx;
                p_z = sz*dz;
                n_x = (float)Math.cos(angle);
                n_z = (float)Math.sin(angle);
                time = 0.0f;
                /* initial raypath */
                float [][] raypath = new float[2][countbndr];

                for(int i=0;i<countbndr;i++){

                        /* cal_gridpoint */
                        ip_lux = (int)(p_x/dx);
                        ip_luz = (int)(p_z/dz);
                        ip_ldx = ip_lux;
                        ip_ldz = ip_luz+1;
                        ip_rux = ip_lux+1;
                        ip_ruz = ip_luz;
                        ip_rdx = ip_lux+1;
                        ip_rdz = ip_luz+1;

                        /* cal_gridvel */
                        float l_x1 = (v0[ip_rux*nz + ip_ruz]-v0[ip_lux*nz + ip_luz])/dx;
                        float l_x2 = (v0[ip_rdx*nz + ip_rdz]-v0[ip_ldx*nz + ip_ldz])/dx;
                        l_x = ( l_x1+l_x2 ) / 2.0f;
                        float l_z1 = (v0[ip_ldx*nz + ip_ldz]-v0[ip_lux*nz + ip_luz])/dz;
                        float l_z2 = (v0[ip_rdx*nz + ip_rdz]-v0[ip_rux*nz + ip_ruz])/dz;
                        l_z = (l_z1+l_z2)/2.0f;
                        v = v0[ip_rux*nz + ip_ruz]+v0[ip_lux*nz + ip_luz]+v0[ip_rdx*nz + ip_rdz]
                                        +v0[ip_ldx*nz + ip_ldz];
                        v = v/4.0f;
                        /* cal_path */
                        float dotmult_ln = l_x*n_x+l_z*n_z;
                        float dotmult_ll = l_x*l_x+l_z*l_z;
                        p_xend = p_x+n_x*s*(1f+dotmult_ln*0.5f*s/v)-0.5f*l_x*s*s/v
                                        -n_x*s*s*s*(dotmult_ll-dotmult_ln*dotmult_ln)/(6.0f*v*v);
                        p_zend = p_z+n_z*s*(1f+dotmult_ln*0.5f*s/v)-0.5f*l_z*s*s/v
                                        -n_z*s*s*s*(dotmult_ll-dotmult_ln*dotmult_ln)/(6.0f*v*v);
                        n_xnew = n_x*(1f+dotmult_ln*s/v)-l_x*s/v
                                        -n_x*s*s*(dotmult_ll-dotmult_ln*dotmult_ln)/(2.0f*v*v);
                        n_znew = n_z*(1f+dotmult_ln*s/v)-l_z*s/v
                                        -n_z*s*s*(dotmult_ll-dotmult_ln*dotmult_ln)/(2.0f*v*v);
                        time += s/v*(1f+s*s*(dotmult_ll+dotmult_ln*dotmult_ln)
                                        /(6f*v*v)-dotmult_ln*s/(2f*v));
                        /* buffer */
                        p_x = p_xend;
                        p_z = p_zend;
                        n_x = n_xnew;
                        n_z = n_znew;
                        raypath[0][i] = p_x;
                        raypath[1][i] = p_z;
                    }
                return raypath;
          }
}
//a########################################################
//a##
//a##  myLaganRayTracingJPanel to myLaganRayTracingJFrame
//a##
//a########################################################
class myLaganRayTracingJPanel extends JPanel{

        private int nx,nz,sx,sz;
        private float dx,dz,angle,s;
        private String FNvel;
        private float[] v0;
        private int nray;
        private float fangle, dangle;

        private int raypathcolor;
        private int raypathwidth;

        private String raypathtext;
        private Writer raypathtextwriter;
        private boolean writeraypath = true;

        private float imax, imin;

        private Lagan lagan;

        private boolean readvel;

        private float constvel;

        myLaganRayTracingJPanel(int nx,int nz, int sx,int sz,float dx,float dz,float s,String FNvel,
                                int nray, float fangle,float dangle,int raypathcolor,int raypathwidth,
                                String raypathtext, boolean readvel, float constvel)throws java.io.IOException{


                this.nx = nx;
                this.nz = nz;
                this.sx = sx;
                this.sz = sz;
                this.dx = dx;
                this.dz = dz;
                this.v0 = new float[nx*nz];
                this.s = s;
                this.FNvel = FNvel;
                this.nray = nray;
                this.fangle = fangle;
                this.dangle = dangle;
                this.raypathcolor = raypathcolor;
                this.raypathwidth = raypathwidth;
                this.raypathtext = raypathtext;

                this.readvel = readvel;

                this.constvel = constvel;

                if(sx>=nx-1)sx=nx-2;
                if(sx<0)sx=0;
                if(sz>=nz-1)sz=nz-2;
                if(sz<0)sz=0;

                raypathtextwriter = new FileWriter(raypathtext);

                readVelFile();
                getMaxMin();
                lagan = new Lagan(sx, sz, nx, nz, dx, dz, v0, s);
          }
        public void paintComponent(Graphics g) {


                super.paintComponent(g);

                Graphics2D g2d = (Graphics2D) g;

                Dimension size = getSize();
                Insets insets = getInsets();

                if(readvel){
                        int[] rgb = new int[3];
                        float radiox = (float)(getWidth())/(float)(nx-0);
                        float radioz = (float)(getHeight())/(float)(nz-0);
                        for (int ix = 0; ix<nx; ix ++)
                        for (int iz = 0; iz<nz; iz ++){

                                int i = ix*nz+iz;
                                /* gray */
                                rgb[0] = rgb[1] = rgb[2] = (int)( 1.0f*(255-(v0[i]-imin)*255.0f/(imax-imin)) );

                                int R = (int)( rgb[0] );
                                int G = (int)( rgb[1] );
                                int B = (int)( rgb[2] );

                                g2d.setColor(new Color(R, G, B) );

                                int drawx = (int)((float)(ix-0)*radiox);
                                int drawz = (int)((float)(iz-0)*radioz);

                                for(int j=0;j<=(int)(radiox);j++)
                                        for(int k=0;k<=(int)(radioz);k++)
                                                g2d.drawLine(drawx+j, drawz+k, drawx+j, drawz+k);
                             }
                    }

                float radioX = (float)(getWidth())/(float)(nx*dx);
                float radioZ = (float)(getHeight())/(float)(nz*dz);
                try{
                        g2d.setColor(getRaypathColor(raypathcolor));
                        g2d.setStroke(getRaypathWidth(raypathwidth));
                        for(int iray=0;iray<nray;iray++){
                                angle = fangle + iray*dangle;
                                float[][] raypath = lagan.makeRay(angle);

                                if(writeraypath){

                                        for(int ii=0;ii<raypath[0].length;ii++){

                                                raypathtextwriter.write(raypath[0][ii]+"         "+
                                                                        raypath[1][ii]+"        \n");
                                                  }
   
                                        }

                                for(int i=0;i<raypath[0].length-1;i++){

                                        int x0 = (int)(radioX*raypath[0][i]);
                                        int x1 = (int)(radioX*raypath[0][i+1]);
                                        int z0 = (int)(radioZ*raypath[1][i]);
                                        int z1 = (int)(radioZ*raypath[1][i+1]);
                                        g2d.drawLine(x0,z0,x1,z1);
                                        }
                                if(writeraypath){
                                        float M = -999999.9f;
                                        raypathtextwriter.write(M+"          "+M+"        \n");
                                        }
                              }
                        if(writeraypath)raypathtextwriter.close();
                        writeraypath = false;

                }catch(Exception ee){ee.printStackTrace();}


          }
        /* get raypath line width */
        public BasicStroke getRaypathWidth(int dim){

                BasicStroke tmp = new BasicStroke(dim);
                return tmp;
          }
        /* get raypath line color */
        public Color getRaypathColor(int dim){

                if(dim==1)return Color.RED;
                else if(dim==2)return Color.BLUE;
                else if(dim==3)return Color.GREEN;
                else if(dim==4)return Color.BLACK;
                else if(dim==5)return Color.GRAY;
                else if(dim==6)return Color.YELLOW;
                else if(dim==7)return Color.PINK;
                else if(dim==8)return Color.CYAN;
                else if(dim==9)return Color.MAGENTA;
                else if(dim==10)return Color.ORANGE;
                else return Color.BLACK;
          }
        /* get min & max */
        public void getMaxMin(){

                imax = imin = v0[0];
                for(int i=0;i<nx*nz;i++){

                        if(imin>v0[i])imin=v0[i];
                        if(imax<v0[i])imax=v0[i];
                    }
          }
        /* read file */
        public void readVelFile(){
                if(readvel){
                        DataInputStream fp = null;
                        try{    
                                if(!new File(FNvel).exists()){  
                                        System.out.println("The "+FNvel+" dont't exists");
                                        return;  
                                      } 
                                fp = new DataInputStream(new FileInputStream(new File(FNvel)));
                                int i = 0;
                                while(fp.available()>0&&i<nx*nz){   
                                        v0[i++] = swap(fp.readFloat());  
                                      } 
                        }catch(Exception e){  
                                e.printStackTrace();  
                              }
                }else{
                        for(int i=0;i<v0.length;i++)
                                v0[i] = constvel;
                    }
          }

        /* swap */
        public static float swap (float value){

                int intValue = Float.floatToRawIntBits (value);
                intValue = swap (intValue);
                return Float.intBitsToFloat (intValue);
          }
        public static int swap (int value){

                int b1 = (value >>  0) & 0xff;
                int b2 = (value >>  8) & 0xff;
                int b3 = (value >> 16) & 0xff;
                int b4 = (value >> 24) & 0xff;
                return b1 << 24 | b2 << 16 | b3 << 8 | b4 << 0;
          }

}
//a########################################################
//a##
//a##  myLaganRayTracingJFrame to myLaganRayTracingParasJFrame
//a##
//a########################################################
class myLaganRayTracingJFrame extends JFrame{

        private Toolkit toolkit;

        private myLaganRayTracingJPanel mylaganraytracingJPanel;

        myLaganRayTracingJFrame(int nx,int nz, int sx,int sz,float dx,float dz,float s,String FNvel,
                                int nray, float fangle,float dangle,int raypathcolor,int raypathwidth,
                                String raypathtext, boolean readvel, float constvel)throws java.io.IOException{

                setTitle("Lagan Raytracing use <"+FNvel+"> as background velocity.Copyright@Rong Tao");
                setSize(500, 600);
                setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

                toolkit = getToolkit();
                Dimension size = toolkit.getScreenSize();
                setLocation((size.width/2 - getWidth())/2, (size.height - getWidth())/2);

                /* menu: save, close */
                JMenuBar JMenuBarSave;
                final JMenu JMenuSave;
                JMenuItem JMenuItemSaveClose,JMenuItemSaveSave;

                JMenuBarSave=new JMenuBar();
                JMenuSave=new JMenu("File");
                JMenuItemSaveSave=new JMenuItem("Save as PNG Image");
                JMenuItemSaveClose=new JMenuItem("Close");

                JMenuSave.add(JMenuItemSaveSave);
                JMenuSave.add(JMenuItemSaveClose);
                JMenuBarSave.add(JMenuSave);

                setJMenuBar(JMenuBarSave);


                mylaganraytracingJPanel = new myLaganRayTracingJPanel(nx,nz, sx,sz,dx,dz,s,FNvel,nray,
                                                 fangle,dangle,raypathcolor,raypathwidth,raypathtext,
                                                 readvel,constvel);

                JMenuItemSaveClose.addActionListener(new ActionListener() {
                        public void actionPerformed(ActionEvent e) {
                                dispose();
                              }
                    });
                JMenuItemSaveSave.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                                saveJPanel.saveImage(mylaganraytracingJPanel, JMenuSave);
                              }
                   });
                add(mylaganraytracingJPanel);
                setVisible(true);
          }
}
//a########################################################
//a##
//a##  myLaganRayTracingParasJFrame
//a##
//a########################################################
class Demo extends JFrame{

        private Toolkit toolkit;
        private boolean flag_SelectVelFile = false;


        private int nx,nz,sx,sz,nray,raypathcolor,raypathwidth;
        private float dx,dz,s,fangle,dangle;
        private String FNvel,raypathtext;

        private float constvel;

        Demo(){

                setTitle("Toa ISO Lagan RayTracing. Copyright@ Rong Tao");
                setSize(700, 460);
                setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

                toolkit = getToolkit();
                Dimension size = toolkit.getScreenSize();
                setLocation((size.width/2 - getWidth())/2, (size.height - getWidth())/2);

                /* menu: About, close */
                JMenuBar menubar = new JMenuBar();
                JMenu menu=new JMenu("File");;
                JMenuItem about=new JMenuItem("About");
                menu.add(about);
                about.addActionListener(new ActionListener() {
                        public void actionPerformed(ActionEvent e) {
                                new myLaganRayTracingAboutDialog();
                              }
                    });
                JMenuItem close =new JMenuItem("Close");
                menu.add(close);
                close.addActionListener(new ActionListener() {
                        public void actionPerformed(ActionEvent e) {
                                dispose();
                              }
                    });
                menubar.add(menu);
                setJMenuBar(menubar);


                /* include  toolbars */
                final JPanel panel = new JPanel();
                panel.setLayout(new GridLayout(6,1));

                /* toolbar0 for velocity*/
                final JToolBar toolbar0 = new JToolBar();
                JCheckBox checkbox = new JCheckBox("Select File[Y/N]", false);
                checkbox.setFocusable(false);
                checkbox.setFont(new Font("Georgia", Font.BOLD, 20));
                final JLabel selectYN = new JLabel(" > NO");
                selectYN.setFont(new Font("Georgia", Font.BOLD, 20));

                toolbar0.add(checkbox);
                toolbar0.add(selectYN);

                ImageIcon folder = new ImageIcon("folder32.png");

                /* toolbar01 for FNvel*/
                final JToolBar toolbar01 = new JToolBar();
                toolbar01.add(new JLabel("Velocity:"));
                final JTextField JFselectvelocity = new JTextField(20);
                JFselectvelocity.setFont(new Font("Georgia", Font.BOLD, 20));
                JFselectvelocity.setText("/home/");
                FNvel = JFselectvelocity.getText();
                toolbar01.add(JFselectvelocity);

                final JButton buttonvelocity = new JButton(folder);
                buttonvelocity.setAlignmentX(1f);
                buttonvelocity.addActionListener(new ActionListener() {
                        public void actionPerformed(ActionEvent event) {

                                JFileChooser fileopen = new JFileChooser();
                                FileFilter filtertxt = new FileNameExtensionFilter("*.dat", "dat");
                                FileFilter filterpar = new FileNameExtensionFilter("*.bin", "bin");
                                fileopen.addChoosableFileFilter(filtertxt);
                                fileopen.addChoosableFileFilter(filterpar);

                                int ret = fileopen.showDialog(panel, "Select");

                                if (ret == JFileChooser.APPROVE_OPTION) {
                                        File file = fileopen.getSelectedFile();
                                        JFselectvelocity.setText(file.toString());
                                        FNvel = file.toString();
                                        System.out.println("Velocity file name : "+FNvel);
                                        }
                              }
                    });
                toolbar01.add(buttonvelocity);

                JFselectvelocity.setEnabled(false);
                buttonvelocity.setEnabled(false);


                /*constant velocity*/

                final JTextField JFconstvel = new JTextField(5);
                JFconstvel.setFont(new Font("Georgia", Font.BOLD, 20));
                JFconstvel.setText("2000");
                JFconstvel.getDocument().addDocumentListener (new DocumentListener() {  
                        @Override  
                        public void insertUpdate(DocumentEvent e) {  
                                if(!isNumeric(JFconstvel.getText()))
                                        JOptionPane.showMessageDialog(JFconstvel, 
                                                   "Input is not digital(0-9)!",
                                                   "Error", JOptionPane.ERROR_MESSAGE);
                                }  
                        @Override  
                        public void removeUpdate(DocumentEvent e) {}   
                        @Override  
                        public void changedUpdate(DocumentEvent e) {}        
                    });

                checkbox.addActionListener(new ActionListener() {
                        public void actionPerformed(ActionEvent e) {

                                flag_SelectVelFile = !flag_SelectVelFile;

                                if (flag_SelectVelFile) {
                                        selectYN.setText(">YES :select the binary file");
                                } else {
                                        selectYN.setText(">NO :just set constant velocity");
                                        }
                                JFselectvelocity.setEnabled(flag_SelectVelFile);
                                buttonvelocity.setEnabled(flag_SelectVelFile);

                                JFconstvel.setEnabled(!flag_SelectVelFile);
                              }
                    });


                /* toolbar1 */
                final JToolBar toolbar1 = new JToolBar();
                /*nx*/
                toolbar1.add(new JLabel(" nx:"));
                final JTextField JFnx = new JTextField(5);
                JFnx.setFont(new Font("Georgia", Font.BOLD, 20));
                JFnx.setText("200");
                final JTextField JFsx = new JTextField(5);
                JFsx.setFont(new Font("Georgia", Font.BOLD, 20));
                JFsx.setText("100");
                JFnx.getDocument().addDocumentListener (new DocumentListener() {  
                        @Override  
                        public void insertUpdate(DocumentEvent e) {  
                                if(!isNumeric(JFnx.getText()))
                                        JOptionPane.showMessageDialog(JFnx, 
                                                   "Input is not digital(0-9)!",
                                                   "Error", JOptionPane.ERROR_MESSAGE);
                                else JFsx.setText( Integer.toString(Integer.parseInt( JFnx.getText() )/2) );
                                }  
                        @Override  
                        public void removeUpdate(DocumentEvent e) {  
                                JFsx.setText( Integer.toString(Integer.parseInt( JFnx.getText() )/2) );
                                }   
                        @Override  
                        public void changedUpdate(DocumentEvent e) {}        
                    });
                toolbar1.add(JFnx);
                /*nz*/
                toolbar1.add(new JLabel(" nz:"));
                final JTextField JFnz = new JTextField(5);
                JFnz.setFont(new Font("Georgia", Font.BOLD, 20));
                JFnz.setText("200");
                final JTextField JFsz = new JTextField(5);
                JFsz.setFont(new Font("Georgia", Font.BOLD, 20));
                JFsz.setText("100");
                JFnz.getDocument().addDocumentListener (new DocumentListener() {  
                        @Override  
                        public void insertUpdate(DocumentEvent e) {  
                                if(!isNumeric(JFnz.getText()))
                                        JOptionPane.showMessageDialog(JFnz, 
                                                   "Input is not digital(0-9)!",
                                                   "Error", JOptionPane.ERROR_MESSAGE);
                                else JFsz.setText( Integer.toString(Integer.parseInt( JFnz.getText() )/2) );
                                }  
                        @Override  
                        public void removeUpdate(DocumentEvent e) {  
                                JFsz.setText( Integer.toString(Integer.parseInt( JFnz.getText() )/2) );
                                }   
                        @Override  
                        public void changedUpdate(DocumentEvent e) {}        
                    });
                toolbar1.add(JFnz);
                /*dx*/
                toolbar1.add(new JLabel(" dx[m]:"));
                final JTextField JFdx = new JTextField(5);
                JFdx.setFont(new Font("Georgia", Font.BOLD, 20));
                JFdx.setText("5");
                JFdx.getDocument().addDocumentListener (new DocumentListener() {  
                        @Override  
                        public void insertUpdate(DocumentEvent e) {  
                                if(!isNumeric(JFdx.getText()))
                                        JOptionPane.showMessageDialog(JFdx, 
                                                   "Input is not digital(0-9)!",
                                                   "Error", JOptionPane.ERROR_MESSAGE);
                                }  
                        @Override  
                        public void removeUpdate(DocumentEvent e) {}   
                        @Override  
                        public void changedUpdate(DocumentEvent e) {}        
                    });
                toolbar1.add(JFdx);
                /*dz*/
                toolbar1.add(new JLabel(" dz[m]:"));
                final JTextField JFdz = new JTextField(5);
                JFdz.setFont(new Font("Georgia", Font.BOLD, 20));
                JFdz.setText("5");
                JFdz.getDocument().addDocumentListener (new DocumentListener() {  
                        @Override  
                        public void insertUpdate(DocumentEvent e) {  
                                if(!isNumeric(JFdz.getText()))
                                        JOptionPane.showMessageDialog(JFdz, 
                                                   "Input is not digital(0-9)!",
                                                   "Error", JOptionPane.ERROR_MESSAGE);
                                }  
                        @Override  
                        public void removeUpdate(DocumentEvent e) {}   
                        @Override  
                        public void changedUpdate(DocumentEvent e) {}        
                    });
                toolbar1.add(JFdz);
                /*sx*/
                toolbar1.add(new JLabel(" sx:"));
                JFsx.getDocument().addDocumentListener (new DocumentListener() {  
                        @Override  
                        public void insertUpdate(DocumentEvent e) {  
                                if(!isNumeric(JFsx.getText())){
                                        JOptionPane.showMessageDialog(JFsx, 
                                                   "Input is not digital(0-9)!",
                                                   "Error", JOptionPane.ERROR_MESSAGE);
                                        return ;
                                        }
                                if(Integer.parseInt( JFnx.getText() )-4<Integer.parseInt( JFsx.getText() ))
                                        JOptionPane.showMessageDialog(JFsx, 
                                                   "You must set: sx <= nx-4",
                                                   "Error", JOptionPane.ERROR_MESSAGE);

                                }  
                        @Override  
                        public void removeUpdate(DocumentEvent e) {}   
                        @Override  
                        public void changedUpdate(DocumentEvent e) {}        
                    });
                toolbar1.add(JFsx);
                /*sz*/
                toolbar1.add(new JLabel(" sz:"));
                JFsz.getDocument().addDocumentListener (new DocumentListener() {  
                        @Override  
                        public void insertUpdate(DocumentEvent e) {  
                                if(!isNumeric(JFsz.getText())){
                                        JOptionPane.showMessageDialog(JFsz, 
                                                   "Input is not digital(0-9)!",
                                                   "Error", JOptionPane.ERROR_MESSAGE);
                                        return ;
                                        }
                                if(Integer.parseInt( JFnz.getText() )-4<Integer.parseInt( JFsz.getText() ))
                                        JOptionPane.showMessageDialog(JFsz, 
                                                   "You must set: sz <= nz-4",
                                                   "Error", JOptionPane.ERROR_MESSAGE);
                                }  
                        @Override  
                        public void removeUpdate(DocumentEvent e) {}   
                        @Override  
                        public void changedUpdate(DocumentEvent e) {}        
                    });
                toolbar1.add(JFsz);


                /* toolbar2 */
                JToolBar toolbar2 = new JToolBar();
                /*nray*/
                toolbar2.add(new JLabel(" nray:"));
                final JTextField JFnray = new JTextField(5);
                JFnray.setFont(new Font("Georgia", Font.BOLD, 20));
                JFnray.setText("90");
                JFnray.getDocument().addDocumentListener (new DocumentListener() {  
                        @Override  
                        public void insertUpdate(DocumentEvent e) {  
                                if(!isNumeric(JFnray.getText()))
                                        JOptionPane.showMessageDialog(JFnray, 
                                                   "Input is not digital(0-9)!",
                                                   "Error", JOptionPane.ERROR_MESSAGE);
                                }  
                        @Override  
                        public void removeUpdate(DocumentEvent e) {}   
                        @Override  
                        public void changedUpdate(DocumentEvent e) {}        
                    });
                toolbar2.add(JFnray);

                /*raypathwidth*/
                String[] LW = {"   1   ","   2   ","   3   ","   4   ",
                                        "   5   ","   6   ","   7   ","   8   ","   9   "};
                final JComboBox comboxLW = new JComboBox<String>(LW);
                comboxLW.setSelectedIndex(0);
                raypathwidth = 1;
                comboxLW.addItemListener(new ItemListener(){

                        public void itemStateChanged(ItemEvent e) {
                                if (e.getStateChange() == ItemEvent.SELECTED) {
                                        JComboBox combo = (JComboBox) e.getSource();
                                        int index = combo.getSelectedIndex();
                                        raypathwidth = index+1;
                                        System.out.println("raypath linewidth = "+raypathwidth);
                                        }
                              }
                    });
                toolbar2.add(new JLabel("Raypath LineWidth:"));
                toolbar2.add(comboxLW);

                /*raypathcolor*/
                String[] LC = {"  red  ","  blue  ","  green  ","  black  "
                              ,"  gray  ","  yellow ","  pink  ","  cyan  ","  magenta  ","  orange  "};
                final JComboBox comboxLC = new JComboBox<String>(LC);
                comboxLC.setSelectedIndex(0);
                raypathcolor = 1;
                comboxLC.addItemListener(new ItemListener(){

                        public void itemStateChanged(ItemEvent e) {
                                if (e.getStateChange() == ItemEvent.SELECTED) {
                                        JComboBox combo = (JComboBox) e.getSource();
                                        int index = combo.getSelectedIndex();
                                        raypathcolor = index+1;
                                        System.out.println("Raypath linecolor = "+raypathcolor);
                                        }
                              }
                    });
                toolbar2.add(new JLabel("Raypath LineColor:"));
                toolbar2.add(comboxLC);

                /*s*/
                toolbar2.add(new JLabel(" s[cm]:"));
                final JTextField JFs = new JTextField(5);
                JFs.setFont(new Font("Georgia", Font.BOLD, 20));
                JFs.setText("50");
                JFs.getDocument().addDocumentListener (new DocumentListener() {  
                        @Override  
                        public void insertUpdate(DocumentEvent e) {  
                                if(!isNumeric(JFs.getText()))
                                        JOptionPane.showMessageDialog(JFs, 
                                                   "Input is not digital(0-9)!",
                                                   "Error", JOptionPane.ERROR_MESSAGE);
                                }  
                        @Override  
                        public void removeUpdate(DocumentEvent e) {}   
                        @Override  
                        public void changedUpdate(DocumentEvent e) {}        
                    });
                toolbar2.add(JFs);


                /* toolbar3 for raypath.txt*/
                final JToolBar toolbar3 = new JToolBar();
                toolbar3.add(new JLabel("Output Raypath Text:"));
                final JTextField JFraypathtext = new JTextField(20);
                JFraypathtext.setFont(new Font("Georgia", Font.BOLD, 20));
                JFraypathtext.setText("raypath.txt");
                raypathtext = JFraypathtext.getText();
                toolbar3.add(JFraypathtext);

                final JButton buttonraypathtext = new JButton(folder);
                buttonraypathtext.setAlignmentX(1f);
                buttonraypathtext.addActionListener(new ActionListener() {
                        public void actionPerformed(ActionEvent event) {

                                JFileChooser fileopen = new JFileChooser();

                                int ret = fileopen.showDialog(panel, "Cover");

                                if (ret == JFileChooser.APPROVE_OPTION) {
                                        File file = fileopen.getSelectedFile();
                                        JFraypathtext.setText(file.toString());
                                        raypathtext = file.toString();
                                        System.out.println("Raypath text file name : "+raypathtext);
                                        }
                              }
                    });
                toolbar3.add(buttonraypathtext);



                /* toolbar4 */
                JToolBar toolbar4 = new JToolBar();
                /*fangle*/
                toolbar4.add(new JLabel(" fangle[o]:"));
                final JTextField JFfangle = new JTextField(5);
                JFfangle.setFont(new Font("Georgia", Font.BOLD, 20));
                JFfangle.setText("45");
                JFfangle.getDocument().addDocumentListener (new DocumentListener() {  
                        @Override  
                        public void insertUpdate(DocumentEvent e) {  
                                if(!isNumeric(JFfangle.getText()))
                                        JOptionPane.showMessageDialog(JFfangle, 
                                                   "Input is not digital(0-9)!",
                                                   "Error", JOptionPane.ERROR_MESSAGE);
                                }  
                        @Override  
                        public void removeUpdate(DocumentEvent e) {}   
                        @Override  
                        public void changedUpdate(DocumentEvent e) {}        
                    });
                toolbar4.add(JFfangle);

                /*dangle*/
                toolbar4.add(new JLabel(" dangle[o]:"));
                final JTextField JFdangle = new JTextField(5);
                JFdangle.setFont(new Font("Georgia", Font.BOLD, 20));
                JFdangle.setText("1");
                JFdangle.getDocument().addDocumentListener (new DocumentListener() {  
                        @Override  
                        public void insertUpdate(DocumentEvent e) {  
                                if(!isNumeric(JFdangle.getText()))
                                        JOptionPane.showMessageDialog(JFdangle, 
                                                   "Input is not digital(0-9)!",
                                                   "Error", JOptionPane.ERROR_MESSAGE);
                                }  
                        @Override  
                        public void removeUpdate(DocumentEvent e) {}   
                        @Override  
                        public void changedUpdate(DocumentEvent e) {}        
                    });
                toolbar4.add(JFdangle);
                toolbar4.add(new JLabel("constVel[m/s]:"));
                toolbar4.add(JFconstvel);

                /*OK button*/
                toolbar4.add(new JLabel("  "));
                ImageIcon iconOK = new ImageIcon("ok32.png");
                JButton buttonOK = new JButton(iconOK);
                buttonOK.setAlignmentX(1f);
                buttonOK.addActionListener(new ActionListener() {
                        public void actionPerformed(ActionEvent event) {

                                nx = Integer.parseInt( JFnx.getText() );
                                nz = Integer.parseInt( JFnz.getText() );
                                sx = Integer.parseInt( JFsx.getText() );
                                sz = Integer.parseInt( JFsz.getText() );
                                dx = (float)Integer.parseInt( JFdx.getText() );
                                dz = (float)Integer.parseInt( JFdz.getText() );

                                nray = Integer.parseInt( JFnray.getText() );
                                fangle = (float)Integer.parseInt( JFfangle.getText() );
                                dangle = (float)Integer.parseInt( JFdangle.getText() );
                                s = (float)(Integer.parseInt( JFs.getText() )/100f);

                                constvel = (float)Integer.parseInt( JFconstvel.getText() );

                                FNvel = JFselectvelocity.getText();
                                raypathtext = JFraypathtext.getText();

                                System.out.println("\nnx = "+nx+",  dx = "+dx+"\n");
                                System.out.println("nz = "+nz+",  dz = "+dz+"\n");
                                System.out.println("sz = "+sz+",  sz = "+sz+"\n");
                                System.out.println("fangle = "+fangle+",  dangle = "+dangle+"\n");
                                System.out.println("raypathcolor = "+raypathcolor
                                                  +",  raypathwidth = "+raypathwidth+"\n");
                                System.out.println("nray = "+nray+"\n");
                                System.out.println("s = "+s+"\n");
                                System.out.println("constant velocity = "+constvel+"\n");
                                System.out.println("velocity file path:"+FNvel+"\n");
                                System.out.println("Output raypath file path:"+raypathtext+"\n");



                                if(sx<0||sx>nx||sz<0||sz>nz){
                                        JOptionPane.showMessageDialog(JFdz, 
                                                 "shot location out of boundary!",
                                                 "Error", JOptionPane.ERROR_MESSAGE);
                                }else{

                                    try{
                                        if(flag_SelectVelFile){
                                                if(FNvel.equals("/home/")){
                                                        JOptionPane.showMessageDialog(JFdz, 
                                                             "Please select velocity file!",
                                                             "Error", JOptionPane.ERROR_MESSAGE);

                                                }else new myLaganRayTracingJFrame(nx,nz,sx,sz,dx,dz,s,FNvel,
                                                                        nray,fangle,dangle,
                                                                        raypathcolor,raypathwidth,
                                                                        raypathtext, true, constvel);
                                        }else{
                                                new myLaganRayTracingJFrame(nx,nz,sx,sz,dx,dz,s,FNvel,
                                                                        nray,fangle,dangle,
                                                                        raypathcolor,raypathwidth,
                                                                        raypathtext, false, constvel);
                                                  }
                                    }catch(Exception ee){ee.printStackTrace();}
                                        }
                              }
                    });
                toolbar4.add(buttonOK);


                panel.add(toolbar0);
                panel.add(toolbar01);
                panel.add(toolbar1);
                panel.add(toolbar2);
                panel.add(toolbar3);
                panel.add(toolbar4);

                add(panel, BorderLayout.NORTH);

                JLabel labelTitle = new JLabel("Click The >OK< Button to GO!",JLabel.CENTER);
                labelTitle.setFont(new Font("Georgia", Font.BOLD, 25));
                labelTitle.setForeground(new Color(50, 50, 25));
                add(labelTitle);

                setVisible(true);
          }
        public static boolean isNumeric(String str){
        /*Copyright: http://javapub.iteye.com/blog/666544*/
                for (int i = 0; i < str.length(); i++){
                        if (!Character.isDigit(str.charAt(i))){
                                return false;
                              }
                    }
                return true;
          }


        public static void main(String[]args)throws java.io.IOException{
                
                new Demo();
                //new myLaganRayTracingJFrame(600,301, 300,0,5f,5f,0.5f,"vel_883_301.dat",90, 45f,1f,1,1,"raypath.txt");
          }
}


//a#################################
//a##
//a##  myLaganRayTracingAboutDialog to myJMenuBar
//a##
//a#################################
class myLaganRayTracingAboutDialog extends JDialog {

        private Toolkit toolkit;

        private int itext =0 ;

        public myLaganRayTracingAboutDialog() {

                setTitle("About Information ");
                setSize(300, 350);
                toolkit = getToolkit();
                Dimension size = toolkit.getScreenSize();
                setLocation((size.width/2 - getWidth())/2, (size.height -getHeight())/2);

                JPanel basic = new JPanel();
                basic.setLayout(new BoxLayout(basic, BoxLayout.Y_AXIS));
                add(basic);

                JPanel topPanel = new JPanel(new BorderLayout(0, 0));
                topPanel.setMaximumSize(new Dimension(450, 0));

                JLabel hint = new JLabel("About Lagan RayTracing");
                hint.setBorder(BorderFactory.createEmptyBorder(0, 25, 0, 0));
                topPanel.add(hint);

                ImageIcon icon = new ImageIcon("raytracing32.png");
                JLabel label = new JLabel(icon);
                label.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
                topPanel.add(label, BorderLayout.EAST);

                JSeparator separator = new JSeparator();
                separator.setForeground(Color.gray);
                topPanel.add(separator, BorderLayout.SOUTH);

                basic.add(topPanel);

                JPanel textPanel = new JPanel(new BorderLayout());
                textPanel.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
                final JTextPane pane = new JTextPane();
                pane.setContentType("text/html");
                final String[] text ={ "<p><b>About The TLaganRayTracing</b></p>" +
                                "<p>. Author: Rong Tao</p>" +
                                "<p>. Time: 2017.7 </p>"+
                                "<p>. @Copyright Rong Tao </p>"+
                                "<p>. Location: @UPC </p>",

                                "<p><b></b></p>" +
                                "<p>.   You just change the velocity,</p>" +
                                "<p>. or use constant velocity.  </p>"+
                                "<p>.   This is a constant Paras. </p>"+
                                "<p>.   </p>"+
                                "<p>.   </p>",

                                "<p><b>Use It Stydy</b></p>" +
                                "<p>.   You can see the raytracing path.</p>" +
                                "<p>. after click OK button.</p>"+
                                "<p>.   You can change the shot location.</p>"+
                                "<p>.  </p>",

                                "<p><b>Hope you enjoy it.</b></p>" +
                                "<p>. Learning makes me happy.</p>" +
                                "<p>. I like fitness.</p>"+
                                "<p>. Wish me Luck.</p>"+
                                "<p>. Good Luck! </p>"};

                pane.setText(text[itext]);
                pane.setEditable(false);
                textPanel.add(pane);

                basic.add(textPanel);

                JPanel boxPanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 20, 0));

                basic.add(boxPanel);

                JPanel bottom = new JPanel(new FlowLayout(FlowLayout.RIGHT));
                JButton previous = new JButton("<<-Previous");

                previous.setMnemonic(KeyEvent.VK_N);
                previous.addActionListener(new ActionListener() {
                        public void actionPerformed(ActionEvent e) {
                                if(itext>0&&itext<=3)pane.setText(text[--itext]);
                              }
                    });
                JButton next = new JButton("Next->>");

                next.setMnemonic(KeyEvent.VK_N);
                next.addActionListener(new ActionListener() {
                        public void actionPerformed(ActionEvent e) {
                                if(itext>=0&&itext<3)pane.setText(text[++itext]);
                              }
                    });
                JButton close = new JButton("Close");
                close.setMnemonic(KeyEvent.VK_C);
                close.addActionListener(new ActionListener() {
                        public void actionPerformed(ActionEvent e) {
                                dispose();
                              }
                    });
                bottom.add(previous);
                bottom.add(next);
                bottom.add(close);

                basic.add(bottom);
                bottom.setMaximumSize(new Dimension(300, 0));

                setResizable(false);

                setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
                setVisible(true);
          }
}
class saveJPanel{

        static void saveImage(JPanel panel, JMenu parent) {
        /*https://stackoverflow.com/questions/19621105/save-image-from-jpanel-after-draw */
                BufferedImage img = new BufferedImage(panel.getWidth(), 
                                                      panel.getHeight(), BufferedImage.TYPE_INT_RGB);
                panel.paint(img.getGraphics());
                JFileChooser jFile = new JFileChooser();
                jFile.showSaveDialog(parent);
                Path pth = jFile.getSelectedFile().toPath();
                JOptionPane.showMessageDialog(parent, pth.toString());
                try {
                        ImageIO.write(img, "png", new File(pth.toString()));
                        System.out.println("panel saved as image");

                } catch (Exception e) {
                        System.out.println("panel not saved" + e.getMessage());
                    }
          }
}

