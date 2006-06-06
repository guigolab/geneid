// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:06 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   GeneViewer.java

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Iterator;
import java.util.LinkedList;
import javax.swing.JFrame;

class GeneViewer$5
    implements ActionListener {

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

    GeneViewer$5() {
        super();
    }
}