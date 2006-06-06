// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:05 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   GenePane.java

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JButton;

class GenePane$ChangeView
    implements ActionListener {

    public void actionPerformed(ActionEvent actionevent) {
        GenePane.access$1002(GenePane.this, !GenePane.access$1000(GenePane.this));
        JButton jbutton = (JButton)actionevent.getSource();
        if(GenePane.access$1000(GenePane.this))
            jbutton.setText("Show proportional view");
        else
            jbutton.setText("Show monospaced View");
        GenePane.access$1100(GenePane.this);
        GenePane.access$1200(GenePane.this).repaint();
    }

    GenePane$ChangeView() {
        super();
    }
}