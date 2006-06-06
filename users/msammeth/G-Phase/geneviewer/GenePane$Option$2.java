// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:05 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   GenePane.java

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JColorChooser;
import javax.swing.JFrame;

class GenePane$Option$2
    implements ActionListener {

    public void actionPerformed(ActionEvent actionevent) {
        java.awt.Color color = JColorChooser.showDialog(null, "Choose UTR Color", GenePane.UTRCOLOR1);
        if(color == null) {
            return;
        } else {
            GenePane.UTRCOLOR1 = color;
            tionPane.setVisible(false);
            cess._mth600(GenePane$Option.this).repaint();
            return;
        }
    }

    GenePane$Option$2() {
        super();
    }
}