// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:06 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   GeneViewer.java

import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import javax.swing.JFrame;

class GeneViewer$4 extends WindowAdapter {

    public void windowClosing(WindowEvent windowevent) {
        String s = ((JFrame)windowevent.getSource()).getTitle();
        GeneViewer.access$100(GeneViewer.this, s);
        pack();
    }

    public void windowClosed(WindowEvent windowevent) {
        String s = ((JFrame)windowevent.getSource()).getTitle();
        GeneViewer.access$100(GeneViewer.this, s);
        pack();
    }

    GeneViewer$4() {
        super();
    }
}