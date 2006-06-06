// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:04 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   GeneFrame.java

import java.awt.Container;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.net.URL;
import javax.swing.*;

public class GeneFrame extends JFrame {

    File gffFile;
    String species;
    GenePane gp;

    public GeneFrame() {
        gffFile = null;
        species = "hs";
        gp = null;
    }

    public void init(File file) {
        try {
            JMenu jmenu = new JMenu("File");
            JMenuItem jmenuitem = new JMenuItem("Close", 67);
            jmenuitem.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent actionevent) {
                    dispose();
                }

             {
                super();
            }
            });
            jmenu.add(jmenuitem);
            JMenuBar jmenubar = new JMenuBar();
            jmenubar.add(jmenu);
            setJMenuBar(jmenubar);
            URL url = file.toURL();
            int i = url.toString().lastIndexOf('.');
            URL url1 = new URL(url.toString().substring(0, i) + ".xml");
            gp = new GenePane();
            gp.createImageIcons(new ApplicationImageIconMaker());
            getContentPane().add(gp);
            gp.init(species, url, url1);
        }
        catch(Exception exception) {
            exception.printStackTrace();
        }
    }
}