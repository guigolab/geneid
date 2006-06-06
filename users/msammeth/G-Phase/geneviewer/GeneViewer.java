// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:07 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   GeneViewer.java

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.*;
import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

public class GeneViewer extends JFrame
    implements ListSelectionListener {

    static double version = 1.1000000000000001D;
    static JPanel framePanel = new JPanel();
    static LinkedList frameList = new LinkedList();
    static JFileChooser fChooser = new JFileChooser(".");
    JScrollPane sPane;
    HashSet cloneNames;
    private File gffFile;

    public GeneViewer() {
        sPane = new JScrollPane();
        cloneNames = new HashSet();
        gffFile = null;
        Container container = getContentPane();
        container.setLayout(new BorderLayout());
        JPanel jpanel = new JPanel();
        JLabel jlabel = new JLabel("GeneView", 0);
        jlabel.setFont(new Font("Times New Roman", 0, 20));
        jlabel.setForeground(Color.red);
        jpanel.add(jlabel);
        jpanel.add(new JLabel("    Version " + version));
        container.add(jpanel, "North");
        jpanel = new JPanel();
        sPane.setPreferredSize(new Dimension(100, 100));
        jpanel.add(sPane);
        container.add(jpanel, "Center");
        framePanel.setLayout(new BoxLayout(framePanel, 1));
        container.add(framePanel, "East");
        jpanel = new JPanel();
        jpanel.add(new JLabel("search:"));
        JTextField jtextfield = new JTextField(10);
        jpanel.add(jtextfield);
        jtextfield.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent actionevent) {
                String s = ((JTextField)actionevent.getSource()).getText().trim().toLowerCase();
                Vector vector = new Vector();
                Object obj = cloneNames.iterator();
                do {
                    if(!((Iterator) (obj)).hasNext())
                        break;
                    String s1 = (String)((Iterator) (obj)).next();
                    if(s.length() == 0 || s1.toLowerCase().indexOf(s) >= 0)
                        vector.add(s1);
                } while(true);
                obj = new JList(vector);
                ((JList) (obj)).setSelectionMode(0);
                ((JList) (obj)).addListSelectionListener(GeneViewer.this);
                sPane.setViewportView(((java.awt.Component) (obj)));
            }

             {
                super();
            }
        });
        container.add(jpanel, "South");
        JMenu jmenu = new JMenu("File");
        jmenu.setMnemonic(70);
        JMenuItem jmenuitem = new JMenuItem("Open", 79);
        jmenuitem.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent actionevent) {
                if(GeneViewer.fChooser.showOpenDialog(GeneViewer.this) != 0) {
                    return;
                } else {
                    gffFile = GeneViewer.fChooser.getSelectedFile();
                    readDataFile();
                    return;
                }
            }

             {
                super();
            }
        });
        jmenu.add(jmenuitem);
        JMenuItem jmenuitem1 = new JMenuItem("Quit", 81);
        jmenuitem1.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent actionevent) {
                System.exit(0);
            }

             {
                super();
            }
        });
        jmenu.add(jmenuitem1);
        JMenuBar jmenubar = new JMenuBar();
        jmenubar.add(jmenu);
        setJMenuBar(jmenubar);
        setTitle("Gene Viewer");
    }

    public void valueChanged(ListSelectionEvent listselectionevent) {
        if(listselectionevent.getValueIsAdjusting())
            return;
        JList jlist = (JList)sPane.getViewport().getView();
        if(jlist.getSelectedIndex() != -1) {
            String s = jlist.getSelectedValue().toString();
            createGeneFrame();
        }
    }

    private void removeGeneFrame(String s) {
        java.util.ListIterator listiterator = frameList.listIterator();
        do {
            if(!listiterator.hasNext())
                break;
            JFrame jframe = (JFrame)listiterator.next();
            if(jframe.getTitle() != s)
                continue;
            listiterator.remove();
            break;
        } while(true);
        java.awt.Component acomponent[] = framePanel.getComponents();
        int i = 0;
        do {
            if(i >= acomponent.length)
                break;
            JButton jbutton = (JButton)acomponent[i];
            if(jbutton.getActionCommand() == s) {
                framePanel.remove(i);
                break;
            }
            i++;
        } while(true);
    }

    public void createGeneFrame() {
        GeneFrame geneframe = new GeneFrame();
        try {
            geneframe.addWindowListener(new WindowAdapter() {

                public void windowClosing(WindowEvent windowevent) {
                    String s = ((JFrame)windowevent.getSource()).getTitle();
                    removeGeneFrame(s);
                    pack();
                }

                public void windowClosed(WindowEvent windowevent) {
                    String s = ((JFrame)windowevent.getSource()).getTitle();
                    removeGeneFrame(s);
                    pack();
                }

             {
                super();
            }
            });
            geneframe.pack();
            geneframe.setVisible(true);
            frameList.addLast(geneframe);
            JButton jbutton = new JButton();
            jbutton.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent actionevent) {
                    String s = actionevent.getActionCommand();
                    java.util.ListIterator listiterator = GeneViewer.frameList.listIterator();
                    do {
                        if(!listiterator.hasNext())
                            break;
                        JFrame jframe = (JFrame)listiterator.next();
                        if(jframe.getTitle() == s) {
                            int i = jframe.getExtendedState();
                            i &= -2;
                            jframe.setExtendedState(i);
                            jframe.requestFocus();
                        }
                    } while(true);
                }

             {
                super();
            }
            });
            framePanel.add(jbutton);
            pack();
        }
        catch(Exception exception) {
            System.out.println(exception);
        }
    }

    public void readDataFile() {
        readDataFile(gffFile.getAbsolutePath());
    }

    public void readDataFile(String s) {
        try {
            FileInputStream fileinputstream = new FileInputStream(s);
            BufferedReader bufferedreader = new BufferedReader(new InputStreamReader(fileinputstream));
            Object obj = null;
            do {
                String s1;
                if((s1 = bufferedreader.readLine()) == null)
                    break;
                s1 = s1.trim();
                if(!s1.startsWith("#")) {
                    StringTokenizer stringtokenizer = new StringTokenizer(s1);
                    cloneNames.add(stringtokenizer.nextToken());
                }
            } while(true);
            Vector vector = new Vector();
            String s2;
            for(Iterator iterator = cloneNames.iterator(); iterator.hasNext(); vector.add(s2))
                s2 = (String)iterator.next();

            JList jlist = new JList(vector);
            jlist.setSelectionMode(0);
            jlist.addListSelectionListener(this);
            sPane.setPreferredSize(null);
            sPane.setViewportView(jlist);
        }
        catch(IOException ioexception) {
            System.err.println(ioexception);
        }
    }

    public static void main(String args[]) {
        GeneViewer geneviewer = new GeneViewer();
        geneviewer.addWindowListener(new WindowAdapter() {

            public void windowClosing(WindowEvent windowevent) {
                System.exit(0);
            }

        });
        geneviewer.pack();
        geneviewer.setVisible(true);
        if(args.length >= 1)
            geneviewer.readDataFile(args[0]);
    }



}