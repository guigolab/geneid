// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:06 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   GeneViewer.java

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JFileChooser;

class GeneViewer$2
    implements ActionListener {

    public void actionPerformed(ActionEvent actionevent) {
        if(GeneViewer.fChooser.showOpenDialog(GeneViewer.this) != 0) {
            return;
        } else {
            GeneViewer.access$002(GeneViewer.this, GeneViewer.fChooser.getSelectedFile());
            readDataFile();
            return;
        }
    }

    GeneViewer$2() {
        super();
    }
}