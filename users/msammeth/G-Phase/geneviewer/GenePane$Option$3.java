// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:05 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   GenePane.java

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JColorChooser;
import javax.swing.JFrame;

class GenePane$Option$3
    implements ActionListener {

    public void actionPerformed(ActionEvent actionevent) {
        Color color = JColorChooser.showDialog(null, "Choose Background Color", GenePane.BGCOLOR);
        if(color == null) {
            return;
        } else {
            GenePane.BGCOLOR = color;
            tionPane.setVisible(false);
            GenePane.access$700(cess._mth600(GenePane$Option.this)).setBackground(GenePane.BGCOLOR);
            GenePane.access$800(cess._mth600(this._cls1.this)).setBackground(GenePane.BGCOLOR);
            GenePane.access$800(cess._mth600(this._cls1.this)).resetAllAlternativeContainersColor();
            GenePane.access$900(cess._mth600(GenePane$Option.this)).setBackground(GenePane.BGCOLOR.darker());
            cess._mth600(this._cls1.this).repaint();
            return;
        }
    }

    GenePane$Option$3() {
        super();
    }
}