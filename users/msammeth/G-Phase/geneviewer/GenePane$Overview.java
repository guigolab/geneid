// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:05 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   GenePane.java

import java.awt.event.ActionEvent;
import javax.swing.JFrame;

class GenePane$Overview extends GenePane$Action {

    public void actionPerformed(ActionEvent actionevent) {
        GenePane.access$100(GenePane.this).setVisible(true);
    }

    public GenePane$Overview() {
        super(GenePane.this, "Overview", GenePane.access$000());
        mnemonic = 'o';
        toolTipText = "Show Overview";
    }
}