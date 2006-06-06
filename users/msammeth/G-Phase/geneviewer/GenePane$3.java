// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:05 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   GenePane.java

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JColorChooser;
import javax.swing.JFrame;

private final class GenePane$3
    implements ActionListener {

    public void actionPerformed(ActionEvent e) {
        Color newColor = JColorChooser.showDialog(null, "Choose Background Color", GenePane.BGCOLOR);
        if(newColor == null) {
            return;
        } else {
            GenePane.BGCOLOR = newColor;
            optionPane.setVisible(false);
            GenePane.access$6(tion.access._mth0(tion.this)).setBackground(GenePane.BGCOLOR);
            GenePane.access$7(tion.access._mth0(tion.this)).setBackground(GenePane.BGCOLOR);
            GenePane.access$7(tion.access._mth0(tion.this)).resetAllAlternativeContainersColor();
            GenePane.access$8(tion.access._mth0(tion.this)).setBackground(GenePane.BGCOLOR.darker());
            tion.access._mth0(tion.this).repaint();
            return;
        }
    }

    GenePane$3() {
    }
}