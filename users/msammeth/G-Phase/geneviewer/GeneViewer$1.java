// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:06 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   GeneViewer.java

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.*;
import javax.swing.*;

class GeneViewer$1
    implements ActionListener {

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

    GeneViewer$1() {
        super();
    }
}