/**
 * Main Code
 * 
 *       Author: Rong Tao
 *     Location: UPC
 *         Time: 2017.04
 *    Modify by: Rong Tao
 */
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
import javax.swing.border.Border;
import javax.swing.border.TitledBorder;
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

import javax.swing.JInternalFrame;
import javax.swing.JTable;
import javax.swing.SwingUtilities;
import javax.swing.table.DefaultTableModel;


import java.text.DateFormat;

import java.util.Calendar;
import java.util.Date;
import java.util.Locale;
import java.util.Scanner;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.Writer;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;

import java.nio.file.Path;

import javax.imageio.ImageIO;

import net.java.dev.designgridlayout.DesignGridLayout;

import java.lang.reflect.Field;

import java.util.ArrayList;
import java.util.Vector;



/**
 *   Package of self
 *             Author:Rong Tao
 *
 */
import src.text.myJFrameTextOpenWithSaveAs;
import src.text.myJFrameTextNewWithEdit;
import src.text.myJFrameTextNew;

import src.about.myAuthorAboutDialog;

import src.seismic.myXmovie1OrderPMLParasJFrame;
import src.seismic.myXmovie2OrderAbsParasJFrame;
import src.seismic.myRungeKuttaRayTracingJPanel;
import src.seismic.saveJPanel;
import src.seismic.RungeKutta;
import src.seismic.myXraypathParasJFrame;
import src.seismic.cjbsegy;
import src.seismic.mySEGYJFrame;
import src.seismic.myRungeKuttaRayTracingJFrameParasJFrame;
import src.seismic.myLaganRayTracingParasJFrame;

import src.paint.myXimageParasJFrame;
import src.paint.myXcurveParasJFrame;
import src.paint.myXMutiCurveParasJFrame;

import src.swap.Swap;
//a##########################################
//a##
//a##              myJFrame -> the basic frame ,first show in screen
//a##
//a##########################################
//a#################################
//a##
//a##       myJFrame
//a##
//a#################################
/**
 * 
 *       Author: Rong Tao
 *     Location: UPC
 *         Time: 2017.04
 *    Modify by: Rong Tao
 */
class myJFrame extends JFrame{

    private Toolkit toolkit;

    myJFrame(String myJFrameName){

        setTitle(myJFrameName);
        setSize(850, 600);
				
				add(new BackgroundPanel((
                new ImageIcon("picture/BackGround/background5.png")).getImage()));
				
				
        toolkit = getToolkit();
        Dimension size = toolkit.getScreenSize();
        setLocation((size.width/2 - getWidth())/2, (size.height -getHeight())/2);

        setDefaultCloseOperation(EXIT_ON_CLOSE);          
    }//myJFrame
		class BackgroundPanel extends JPanel  
		{  
			Image image;  
			public BackgroundPanel(Image image)  
			{  
				this.image = image;  
				this.setOpaque(false);//false透明
			}  

			public void paintComponent(Graphics g)      
			{  
				super.paintComponents(g);  
				g.drawImage(image,0,0,this.getWidth(),this.getHeight(),this); 

			}  
		}
}
//a#########################################
//a##
//a##        myJMenuBar -> the menu in basic frame (in the head of basic frame)
//a##
//a##        myJPopupMenu -> the menu in basic frame (popup)
//a##
//a##        myAuthorAboutDialog -> aboutdialog for " menu->help->About "
//a##
//a#########################################
//a#################################
//a##
//a##       myJMenuBar
//a##
//a#################################
/**
 * 
 *       Author: Rong Tao
 *     Location: UPC
 *         Time: 2017.04
 *    Modify by: Rong Tao
 */
class myJMenuBar extends JMenuBar{

    private JLabel statusbar;

    private Toolkit toolkit;

    private static JLabel segycountlabel;

    myJMenuBar(final myJFrame frame){



/*******************File***********************{*/
/*******************File***********************{*/
        final JMenu file = new JMenu("File");
        file.setMnemonic(KeyEvent.VK_F);

        /* New */
        /* New */
        ImageIcon New = new ImageIcon("picture/ButtonImg/newdocument16.png");
        JMenuItem fileNew = new JMenuItem("New Text", New);
        fileNew.setMnemonic(KeyEvent.VK_N);
        fileNew.setToolTipText("New Text");
        fileNew.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                javax.swing.SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        new myJFrameTextNew();
                          }
                    });
                  }
            });


        /* New with edit */
        /* New with edit */
        ImageIcon Newwithedit = new ImageIcon("picture/ButtonImg/newdocumentwithedit16.png");
        JMenuItem fileNewwithedit = new JMenuItem("New Text With Edit", Newwithedit);
        fileNewwithedit.setMnemonic(KeyEvent.VK_N);
        fileNewwithedit.setToolTipText("New Text With Edit");
        fileNewwithedit.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                javax.swing.SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        new myJFrameTextNewWithEdit();

                          }
                    });
                  }
            });



        /* Open */
        /* Open */
        ImageIcon Open = new ImageIcon("picture/ButtonImg/opendocument16.png");
        JMenuItem fileOpen = new JMenuItem("Open Text", Open);
        fileOpen.setMnemonic(KeyEvent.VK_O);
        fileOpen.setToolTipText("Open Text");
        fileOpen.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                javax.swing.SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        new myJFrameTextOpenWithSaveAs(file);
                          }
                    });
                  }
            });



        /* Close */
        /* Close */
        ImageIcon Close = new ImageIcon("picture/ButtonImg/close16.png");
        JMenuItem fileClose = new JMenuItem("Close Software", Close);
        fileClose.setMnemonic(KeyEvent.VK_C);
        fileClose.setToolTipText("Exit/Close Software");
        fileClose.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                System.exit(0);
                  }
            });

        /* SubMenu */
        /* SubMenu */
        JMenu imp = new JMenu("Import");
        imp.setMnemonic(KeyEvent.VK_M);

        JMenuItem newsf = new JMenuItem("Import newsfeed list...");
        JMenuItem bookm = new JMenuItem("Import bookmarks...");
        JMenuItem mail = new JMenuItem("Import mail...");
        imp.add(newsf);
        imp.add(bookm);
        imp.add(mail);



        file.add(fileNew);
        file.add(fileNewwithedit);
        file.add(fileOpen);
        file.addSeparator();
        file.add(imp);
        file.addSeparator();
        file.add(fileClose);


/*******************File***********************}*/
/*******************File***********************}*/




/*******************Edit***********************{*/
/*******************Edit***********************{*/
        JMenu edit = new JMenu("Edit");
        edit.setMnemonic(KeyEvent.VK_E);

        /* Copy */
        /* Copy */
        ImageIcon Copy = new ImageIcon("picture/ButtonImg/copy16.png");
        JMenuItem editCopy = new JMenuItem("Copy", Copy);
        editCopy.setMnemonic(KeyEvent.VK_C);
        editCopy.setToolTipText("Copy");

        /* Paste */
        /* Paste */
        ImageIcon Paste = new ImageIcon("picture/ButtonImg/paste16.png");
        JMenuItem editPaste = new JMenuItem("Paste", Paste);
        editPaste.setMnemonic(KeyEvent.VK_V);
        editPaste.setToolTipText("Paste");

        /* Cut */
        /* Cut */
        ImageIcon Cut = new ImageIcon("picture/ButtonImg/cut16.png");
        JMenuItem editCut = new JMenuItem("Cut", Cut);
        editCut.setMnemonic(KeyEvent.VK_X);
        editCut.setToolTipText("Cut");

        edit.add(editCopy);
        edit.add(editPaste);
        edit.add(editCut);


/*******************Edit***********************}*/
/*******************Edit***********************}*/





/*******************View***********************{*/
/*******************View***********************{*/
        final JMenu view = new JMenu("View");
        view.setMnemonic(KeyEvent.VK_V);

        JMenu setMenu = new JMenu("Menu Font/Color");

        /* Set Menu Color */
        /* Set Menu Color */
        ImageIcon setMenuColor = new ImageIcon("picture/ButtonImg/color16.png");
        JMenuItem setmenucolor = new JMenuItem("Menu Color", setMenuColor);
        setmenucolor.setMnemonic(KeyEvent.VK_C);
        setmenucolor.setToolTipText("Set Menu Background Color");

        setmenucolor.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                JColorChooser clr = new JColorChooser();
                Color color = clr.showDialog(view, 
                                 "Choose Color",Color.white);
                setBackground(color);
                  }
            });
        setMenu.add(setmenucolor);


        /* Set StatusBar Color */
        /* Set StatusBar Color */
        ImageIcon setStatusBarColor = new ImageIcon("picture/ButtonImg/color16.png");
        JMenuItem setstatusbarcolor = new JMenuItem("StatusBar Color", 
                                 setStatusBarColor);
        setstatusbarcolor.setMnemonic(KeyEvent.VK_C);
        setstatusbarcolor.setToolTipText("StatusBar Color");

        setstatusbarcolor.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                JColorChooser clr = new JColorChooser();
                Color color = clr.showDialog(view, "Choose Color",Color.white);
                statusbar.setForeground(color);
                statusbar.setBackground(color);
                  }
            });

        setMenu.add(setstatusbarcolor);
        setMenu.addSeparator();


        /* Set Menu Font */
        /* Set Menu Font */
        JMenu SetMenuFont = new JMenu("Menu Font");
        SetMenuFont.setToolTipText("Menu Font");
        GraphicsEnvironment ge = GraphicsEnvironment.getLocalGraphicsEnvironment();
        String[] fonts = ge.getAvailableFontFamilyNames();
        final JList list = new JList<String>(fonts);
        list.addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent e) {
                if (!e.getValueIsAdjusting()) {
                    String name = (String) list.getSelectedValue();
                    Font font = new Font(name, Font.PLAIN, 20);


                    getComponent(0).setFont(font);
                    getComponent(1).setFont(font);
                    getComponent(2).setFont(font);
                    getComponent(3).setFont(font);
                    statusbar.setFont(font);
                    }
                  }
            });
        JScrollPane pane = new JScrollPane();
        pane.getViewport().add(list);
        pane.setPreferredSize(new Dimension(250, 200));
        SetMenuFont.add(pane);

        setMenu.add(SetMenuFont);


        /* Set Software Background */
        /* Set Software Background */
        JMenuItem Background = new JMenuItem("Background");
        Background.setMnemonic(KeyEvent.VK_C);
        Background.setToolTipText("Set Software Background");

        Background.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {

                new setBackgroundFrame(frame,"picture/BackGround/");
                  }
            });




        /* Show StatuBar */
        /* Show StatuBar */
        JCheckBoxMenuItem CheckBoxMenuItemStatuBar = new JCheckBoxMenuItem("Show StatuBar");
        CheckBoxMenuItemStatuBar.setState(true);

        CheckBoxMenuItemStatuBar.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                if (statusbar.isVisible()) {
                    statusbar.setVisible(false);
                } else {
                    statusbar.setVisible(true);
                    }
                  }
            });

        view.addSeparator();
        view.add(setMenu);
        view.add(Background);
        view.addSeparator();
        view.add(CheckBoxMenuItemStatuBar);




/*******************View***********************}*/
/*******************View***********************}*/






/*******************Help***********************{*/
/*******************Help***********************{*/

        final JMenu help = new JMenu("Help");
        help.setMnemonic(KeyEvent.VK_H);



        /* Error */
        /* Error */
        JMenuItem error = new JMenuItem("Error..");
        error.setToolTipText("Show Error");
        help.add(error);
        error.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                JOptionPane.showMessageDialog(help, "Error.",
                    "Error", JOptionPane.ERROR_MESSAGE);
                  }
            });

        /* Warning */
        /* Warning */
        JMenuItem warning = new JMenuItem("Warn..");
        warning.setToolTipText("Show Warning");
        help.add(warning);
        warning.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                JOptionPane.showMessageDialog(help, "Warning.",
                    "Warning", JOptionPane.WARNING_MESSAGE);
                  }
            });

        /* Question */
        /* Question */
        JMenuItem question = new JMenuItem("Ques..");
        question.setToolTipText("Show Question");
        help.add(question);
        question.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                JOptionPane.showMessageDialog(help, "Question.",
                    "Question", JOptionPane.QUESTION_MESSAGE);
                  }
            });

        /* Information */
        /* Information */
        JMenuItem information = new JMenuItem("Info..");
        information.setToolTipText("Show Information");
        help.add(information);
        information.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                JOptionPane.showMessageDialog(help, "Information.",
                    "Information", JOptionPane.INFORMATION_MESSAGE);
                  }
            });

        /* About */
        /* About */
        JMenuItem about = new JMenuItem("About...");
        about.setToolTipText("Show About Dialog");
        help.addSeparator();
        help.add(about);
        about.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                javax.swing.SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        new myAuthorAboutDialog();
                          }
                    });
                  }
            });

/*******************Help***********************}*/
/*******************Help***********************}*/




/*******************SEGY***********************{*/
/*******************SEGY***********************{*/


        final JTextField segytext = new JTextField("see the Segy header");
        segytext.setFont(new Font("Georgia",Font.BOLD,15));
        segytext.setBackground(Color.YELLOW);
        JButton segybutton 
                = new JButton(
                   new ImageIcon("picture/ButtonImg/folder16.png"));
        segybutton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {

                JFileChooser fileopen = new JFileChooser();
                fileopen.setFileFilter(
                        new FileNameExtensionFilter("*.su", "su")); 
                int ret = fileopen.showDialog(help, "select");

                if (ret == JFileChooser.APPROVE_OPTION) {
                    File file = fileopen.getSelectedFile();
                    segytext.setText(file.toString());
                    }
                  }
            });
        JButton segybuttonOK = new JButton(
                    new ImageIcon("picture/ButtonImg/ok16.png"));
        segybuttonOK.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {

                JOptionPane.showMessageDialog(help,
                     "file: "+segytext.getText()
                    +"\nIt may take some time."
                    +"\nPlease wait patiently!",
                    "Information", JOptionPane.INFORMATION_MESSAGE);

                try{
                    new mySEGYJFrame(segytext.getText(),segycountlabel);
                }catch(Exception ee){ee.printStackTrace();}
                  }
            });

        segycountlabel = new JLabel();
        segycountlabel.setFont(new Font("Georgia",Font.BOLD,15));





/*******************SEGY***********************}*/
/*******************SEGY***********************}*/






/*******************Statusbar***********************{*/
/*******************Statusbar***********************{*/


        statusbar = new JLabel(" Status: ");
        statusbar.setBorder(BorderFactory.createEtchedBorder(EtchedBorder.RAISED));
        statusbar.setFont(new Font("Georgia",Font.BOLD,20));
        statusbar.setForeground(Color.RED);
        /* Timer to show time */
        final Timer timer = new Timer(1000,null);
        timer.start();
        timer.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                Locale locale = Locale.getDefault();
                Date date = new Date();
                String s = DateFormat.getTimeInstance(DateFormat.SHORT,locale)
                             .format(date);

                statusbar.setText(s+"  "
                           +"[@Copyright: RongTao @UPC LEON_SWPI]");
                  }
            });




/*******************Statusbar***********************}*/
/*******************Statusbar***********************}*/



        add(file);
        add(edit);
        add(view);
        add(help);
        add(new JLabel("  "));
        add(segytext);
        add(new JLabel(" "));
        add(segybutton);
        add(new JLabel(" "));
        add(segybuttonOK);
        add(new JLabel(" "));
        add(segycountlabel);
        add(new JLabel("               "));

        getComponent(0).setFont(new Font("Georgia",Font.BOLD,20));
        getComponent(1).setFont(new Font("Georgia",Font.BOLD,20));
        getComponent(2).setFont(new Font("Georgia",Font.BOLD,20));
        getComponent(3).setFont(new Font("Georgia",Font.BOLD,20));
        statusbar.setFont(new Font("Georgia",Font.BOLD,20));



        frame.setJMenuBar(this);
        frame.add(statusbar, BorderLayout.SOUTH);

    }//myJMenuBar


    //a#################################
    //a##
    //a##  setBackgroundFrame  --> innerclass
    //a##
    //a#################################
    /**
     * 
     *       Author: Rong Tao
     *     Location: UPC
     *         Time: 2017.04
     *    Modify by: Rong Tao
     */
    class setBackgroundFrame extends JFrame{

        setBackgroundFrame(myJFrame myjframe,String dir){

            setTitle("Change background");
            setSize(500,400);

            toolkit = getToolkit();
            Dimension size = toolkit.getScreenSize();
            setLocation((size.width/2 - getWidth())/2, (size.height -getHeight())/2);


            File file = new File(dir);
            File [] files = file.listFiles();
            for (int i = 0; i < files.length; i++) {
                File file1 = files[i];
                file1.getName();   
                System.out.println(dir+file1.getName());
                  }

			     myjframe.add(new BackgroundPanel((
                new ImageIcon(dir+files[3].getName())).getImage()));


            setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE); 
            setVisible(true);
             }
        class BackgroundPanel extends JPanel   {  
			Image image;  
			public BackgroundPanel(Image image)  
			{  
				this.image = image;  
				this.setOpaque(false);//false
			}  

			public void paintComponent(Graphics g)      
			{  
				super.paintComponents(g);  
				g.drawImage(image,0,0,this.getWidth(),this.getHeight(),this); 

			}  
            }

    }//end of setBackgroundFrame







}//myJMenuBar


//a#################################
//a##
//a##      myJPopupMenu
//a##
//a#################################
/**
 * 
 *       Author: Rong Tao
 *     Location: UPC
 *         Time: 2017.04
 *    Modify by: Rong Tao
 */
class myJPopupMenu extends JPopupMenu{

    myJPopupMenu(final myJFrame frame){


        /*New Text*/
        /*New Text*/
        JMenuItem menuItemNew = new JMenuItem("New Text");
        menuItemNew.setToolTipText("build a text file");
        menuItemNew.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                javax.swing.SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        new myJFrameTextNew();
                          }
                    });
                  }
            });

        add(menuItemNew);


        /*New Text with Edit*/
        /*New Text with Edit*/
        JMenuItem menuItemNewwithEdit = new JMenuItem("New Text with Edit");
        menuItemNewwithEdit.setToolTipText("build a text file with edit");
        menuItemNewwithEdit.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                javax.swing.SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        new myJFrameTextNewWithEdit();
                          }
                    });
                  }
            });

        add(menuItemNewwithEdit);


        addSeparator();

        /*Set PopupMenu Color*/
        /*Set PopupMenu Color*/
        final JMenuItem setmenucolor = new JMenuItem("Set PopupMenu Color");
        setmenucolor.setMnemonic(KeyEvent.VK_C);
        setmenucolor.setToolTipText("Set Menu Background Color");

        setmenucolor.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                JColorChooser clr = new JColorChooser();
                Color color = clr.showDialog(frame, 
                               "Choose Color",Color.white);
                setBackground(color);
                  }
            });
        add(setmenucolor);

        addSeparator();


        /*quit*/
        /*quit*/
        JMenuItem menuItemClose = new JMenuItem("Close");
        menuItemClose.setToolTipText("colse the software");
        menuItemClose.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                System.exit(0);
                  }
            });

        add(menuItemClose);

        frame.addMouseListener(new MouseAdapter() {
            public void mouseReleased(MouseEvent e) {
                if (e.getButton() == e.BUTTON3) {
                    show(e.getComponent(), e.getX(), e.getY());
                    System.out.println("Popup Menu Located ("
                               +e.getX()+", "+e.getY()+")");
                    }
                  }
            });

    }//myJPopupMenu
}

    
//a#########################################
//a##
//a##        myJToolBarsHorizontal -> the Horizontal toolbar in the top of basic frame
//a##
//a##        myJToolBarsVertical -> the Vertical toolbar in the left of basic frame
//a##
//a#########################################   
//a#################################
//a##
//a##     myJToolBarsHorizontal
//a##
//a#################################
/**
 * 
 *       Author: Rong Tao
 *     Location: UPC
 *         Time: 2017.04
 *    Modify by: Rong Tao
 */
class myJToolBarsHorizontal extends JPanel{

    myJToolBarsHorizontal(myJFrame frame){

        final JToolBar toolbar1 = new JToolBar();
        final JToolBar toolbar2 = new JToolBar();
        final JToolBar toolbar3 = new JToolBar();
        final JToolBar toolbar4 = new JToolBar();

        toolbar1.setBackground(Color.RED);
        toolbar2.setBackground(Color.ORANGE);
        toolbar3.setBackground(Color.BLUE);
        toolbar4.setBackground(Color.GREEN);

        //setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));
        setLayout(new GridLayout(2,2));

        ImageIcon new1 = new ImageIcon(
                     getClass()
                     .getResource("picture/ButtonImg/newdocument32.png"));
        ImageIcon new2 = new ImageIcon( 
                     getClass()
                     .getResource("picture/ButtonImg/newdocumentwithedit32.png"));
        ImageIcon open = new ImageIcon( 
                     getClass()
                     .getResource("picture/ButtonImg/opendocument32.png"));
        ImageIcon save = new ImageIcon( 
                     getClass()
                     .getResource("picture/ButtonImg/save32.png"));
        ImageIcon exit = new ImageIcon( 
                     getClass()
                     .getResource("picture/ButtonImg/exit32.png"));
        ImageIcon color = new ImageIcon( 
                     getClass()
                     .getResource("picture/ButtonImg/color32.png"));
        ImageIcon image = new ImageIcon( 
                     getClass()
                     .getResource("picture/ButtonImg/timage32.png"));
        ImageIcon curve = new ImageIcon( 
                     getClass()
                     .getResource("picture/ButtonImg/tcurve32.png"));
        ImageIcon multicurve = new ImageIcon( 
                     getClass()
                     .getResource("picture/ButtonImg/tcurves32.png"));
        ImageIcon movie2order = new ImageIcon( 
                     getClass()
                     .getResource("picture/ButtonImg/tmovie2order32.png"));
        ImageIcon movie1order = new ImageIcon( 
                     getClass()
                     .getResource("picture/ButtonImg/tmovie1order32.png"));
        ImageIcon raypath = new ImageIcon( 
                     getClass()
                     .getResource("picture/ButtonImg/traypath32.png"));
        ImageIcon raytracingLG = new ImageIcon( 
                     getClass()
                     .getResource("picture/ButtonImg/traytracingLG32.png"));
        ImageIcon raytracingRK = new ImageIcon( 
                     getClass()
                     .getResource("picture/ButtonImg/traytracingRK32.png"));


        JButton newb1 = new JButton(new1);
        JButton newb2 = new JButton(new2);
        JButton openb = new JButton(open);
        JButton saveb = new JButton(save);
        JButton ximageb = new JButton(image);
        JButton xcurveb = new JButton(curve);
        JButton xmulticurveb = new JButton(multicurve);
        JButton xmovie2orderb = new JButton(movie2order);
        JButton xmovie1orderb = new JButton(movie1order);
        JButton xraypathb = new JButton(raypath);
        JButton xraytracingLGb = new JButton(raytracingLG);
        JButton xraytracingRKb = new JButton(raytracingRK);
        JButton exitb = new JButton(exit);
        JButton colorb1 = new JButton(color);
        JButton colorb2 = new JButton(color);
        JButton colorb3 = new JButton(color);
        JButton colorb4 = new JButton(color);



        newb1.setToolTipText("New Text");
        newb2.setToolTipText("New Text With Edit");
        openb.setToolTipText("Open Text");
        saveb.setToolTipText("Save");
        ximageb.setToolTipText("ToaXimage");
        xcurveb.setToolTipText("ToaXcurve");
        xmulticurveb.setToolTipText("ToaXmulticurve");
        xmovie2orderb.setToolTipText("ToaXmovie2order");
        xmovie1orderb.setToolTipText("ToaXmovie1order");
        xraypathb.setToolTipText("ToaXraypath");
        xraytracingLGb.setToolTipText("ToaXrayTracingLagan");
        xraytracingRKb.setToolTipText("Toa-VTI-rayTracingRunge-Kutta");
        exitb.setToolTipText("Exit Software");
        colorb1.setToolTipText("Set this Toolbar Color");
        colorb2.setToolTipText("Set this Toolbar Color");
        colorb3.setToolTipText("Set this Toolbar Color");
        colorb4.setToolTipText("Set this Toolbar Color");



        newb1.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                javax.swing.SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        new myJFrameTextNew();
                          }
                    });
                  }
            });
        newb2.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                javax.swing.SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        new myJFrameTextNewWithEdit();
                          }
                    });
                  }
            });

        openb.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                javax.swing.SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        new myJFrameTextOpenWithSaveAs(toolbar1);
                          }
                    });
                  }
            });

        ximageb.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                javax.swing.SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        new myXimageParasJFrame();
                          }
                    });
                  }
            });

        xcurveb.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                javax.swing.SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        new myXcurveParasJFrame();
                          }
                    });
                  }
            });
        xmulticurveb.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                javax.swing.SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        new myXMutiCurveParasJFrame();
                          }
                    });
                  }
            });

        xmovie2orderb.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                javax.swing.SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        new myXmovie2OrderAbsParasJFrame();
                          }
                    });
                  }
            });

        xmovie1orderb.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                javax.swing.SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        new myXmovie1OrderPMLParasJFrame();
                          }
                    });
                  }
            });

        xraypathb.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                javax.swing.SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        new myXraypathParasJFrame();
                          }
                    });
                  }
            });

        xraytracingLGb.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                javax.swing.SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        new myLaganRayTracingParasJFrame();
                          }
                    });
                  }
            });

        xraytracingRKb.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                javax.swing.SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        new myRungeKuttaRayTracingJFrameParasJFrame();
                          }
                    });
                  }
            });




        colorb1.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                JColorChooser clr = new JColorChooser();
                Color color = clr.showDialog(toolbar1, 
                                 "Choose Color",Color.white);
                toolbar1.setBackground(color);
                  }
            });

        colorb2.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                JColorChooser clr = new JColorChooser();
                Color color = clr.showDialog(toolbar1,  
                                 "Choose Color",Color.white);
                toolbar2.setBackground(color);
                  }
            });
        colorb3.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                JColorChooser clr = new JColorChooser();
                Color color = clr.showDialog(toolbar1,  
                                 "Choose Color",Color.white);
                toolbar3.setBackground(color);
                  }
            });

        colorb4.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                JColorChooser clr = new JColorChooser();
                Color color = clr.showDialog(toolbar1,  
                                 "Choose Color",Color.white);
                toolbar4.setBackground(color);
                  }
            });

        exitb.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                System.exit(0);
                  }
            });

        toolbar1.add(colorb1);
        toolbar1.add(newb1);
        toolbar1.add(newb2);
        toolbar1.add(openb);
        toolbar1.add(saveb);
        toolbar1.add(exitb);
        toolbar1.setAlignmentX(0);

        toolbar2.add(colorb2);
        toolbar2.add(ximageb);
        toolbar2.add(xcurveb);
        toolbar2.add(xmulticurveb);
        toolbar2.add(xmovie2orderb);
        toolbar2.add(xmovie1orderb);
        toolbar2.add(xraypathb);
        toolbar2.add(xraytracingLGb);
        toolbar2.add(xraytracingRKb);
        toolbar2.setAlignmentX(0);

        toolbar3.add(colorb3);
        toolbar3.setAlignmentX(0);

        toolbar4.add(colorb4);
        toolbar4.setAlignmentX(0);


        add(toolbar1);
        add(toolbar2);
        add(toolbar3);
        add(toolbar4);

        frame.add(this, BorderLayout.NORTH);

    }//myJToolBar
}
//a#################################
//a##
//a##     myJToolBarsVertical
//a##
//a#################################
/**
 * 
 *       Author: Rong Tao
 *     Location: UPC
 *         Time: 2017.04
 *    Modify by: Rong Tao
 */
class myJToolBarsVertical extends JPanel{

    myJToolBarsVertical(myJFrame frame){

        setLayout(new GridLayout(2,1));

        final JToolBar toolbar = new JToolBar(JToolBar.VERTICAL);

        toolbar.setBackground(Color.CYAN);

        ImageIcon i1 = new ImageIcon( 
                     getClass() 
                     .getResource("picture/ButtonImg/select32.png"));
        ImageIcon i2 = new ImageIcon( 
                     getClass() 
                     .getResource("picture/ButtonImg/person32.png"));
        ImageIcon i3 = new ImageIcon( 
                     getClass() 
                     .getResource("picture/ButtonImg/timage32.png"));
        ImageIcon i4 = new ImageIcon( 
                     getClass() 
                     .getResource("picture/ButtonImg/newdocument32.png"));
        ImageIcon i5 = new ImageIcon( 
                     getClass() 
                     .getResource("picture/ButtonImg/newdocumentwithedit32.png"));
        ImageIcon i6 = new ImageIcon( 
                     getClass() 
                     .getResource("picture/ButtonImg/ellipse32.png"));
        ImageIcon i7 = new ImageIcon( 
                     getClass() 
                     .getResource("picture/ButtonImg/text32.png"));
        ImageIcon i8 = new ImageIcon( 
                     getClass() 
                     .getResource("picture/ButtonImg/textwithedit32.png"));
        ImageIcon i9 = new ImageIcon( 
                     getClass() 
                     .getResource("picture/ButtonImg/color32.png"));

        JButton b1 = new JButton(i1);
        b1.setToolTipText("Select");
        JButton b2 = new JButton(i2);
        b2.setToolTipText("SoftWare About");
        JButton b3 = new JButton(i3);
        b3.setToolTipText("Image The Binary Data");
        JButton b4 = new JButton(i4);
        b4.setToolTipText("New Text");
        JButton b5 = new JButton(i5);
        b5.setToolTipText("Rectangle");
        JButton b6 = new JButton(i6);
        b6.setToolTipText("Ellipse");
        JButton b7 = new JButton(i7);
        b7.setToolTipText("PS");
        JButton b8 = new JButton(i8);
        b8.setToolTipText("New Text");
        JButton b9 = new JButton(i9);
        b9.setToolTipText("Set Toolbar Color");


        b2.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                javax.swing.SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        new myAuthorAboutDialog();
                          }
                    });
                  }
            });

        b3.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                javax.swing.SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        new myXimageParasJFrame();
                          }
                    });
                  }
            });

        b4.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                javax.swing.SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        new myJFrameTextNew();
                          }
                    });
                  }
            });
        b5.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                javax.swing.SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        new myJFrameTextNewWithEdit();
                          }
                    });
                  }
            });
        b7.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                javax.swing.SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        new myJFrameTextNew();
                          }
                    });
                  }
            });
        b8.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                javax.swing.SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        new myJFrameTextNewWithEdit();
                          }
                    });
                  }
            });

        b9.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent event) {
                JColorChooser clr = new JColorChooser();
                Color color = clr.showDialog(toolbar, 
                                 "Choose Color",Color.white);
                toolbar.setBackground(color);
                  }
            });


        //toolbar.add(b1);
        toolbar.add(b2);
        toolbar.add(b3);
        toolbar.add(b4);
        toolbar.add(b5);
        //toolbar.add(b6);
        toolbar.add(b7);
        toolbar.add(b8);
        toolbar.add(b9);

        add(toolbar);

        frame.add(toolbar, BorderLayout.WEST);

    }//myJToolBarsVertical
}

//a#########################################
//a##
//a##           myBorderJPanel
//a##
//a######################################### 
/**
 * 
 *       Author: Rong Tao
 *     Location: UPC
 *         Time: 2017.04
 *    Modify by: Rong Tao
 */
class myBorderJPanel extends JPanel{

    private TitledBorder titledborder;

    private JPanel panel = new JPanel();

    private GridBagConstraints gc = new GridBagConstraints();

    private GridBagLayout gridbag = new GridBagLayout();

    myBorderJPanel(String title){

        titledborder = BorderFactory.createTitledBorder(
                BorderFactory.createEtchedBorder(
                  EtchedBorder.LOWERED), title);
        titledborder.setTitleFont(new Font("Georgia",Font.BOLD,20));
        titledborder.setTitleJustification(TitledBorder.LEFT);
        titledborder.setTitlePosition(TitledBorder.DEFAULT_POSITION);

        setBorder(BorderFactory.createEmptyBorder(6,6,6,6));

        panel.setLayout(gridbag);
        panel.setBorder(titledborder);

        JPanel tmp = new JPanel();
        tmp.add(panel);
        add(tmp);
    }
    void addC(JComponent c, int gridx, int gridy,int gridwidth,int gridheight){

        gc.gridx = gridx; //x grid position
        gc.gridy = gridy; //y grid position
        gc.gridwidth = gridwidth;
        gc.gridheight = gridheight;
        gridbag.setConstraints(c, gc);
        c.setFont(new Font("Georgia",Font.BOLD,25));
        panel.add(c);
    }
}

//a#########################################
//a##
//a##           main class
//a##
//a######################################### 
/**
 * 
 *       Author: Rong Tao
 *     Location: UPC
 *         Time: 2017.04
 *    Modify by: Rong Tao
 */
public class MainWindow{

    public static void main(String[]args){

        myJFrame frame = new myJFrame("@Copyright: Rong Tao >> "
                        +"UPC.LEON_SWPI >> (QQ:2386499836@qq.com)");

        new myJMenuBar(frame);
        new myJPopupMenu(frame);
        new myJToolBarsHorizontal(frame);
        new myJToolBarsVertical(frame);


        frame.setVisible(true); 
        //frame.pack();
        //frame.setResizable(false);
    }//main
}
