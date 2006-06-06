// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:06 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   GenePane.java

import java.awt.event.ActionEvent;

class GenePane$ScrollRight extends GenePane$Action {

    public void actionPerformed(ActionEvent actionevent) {
        GenePane.access$300(GenePane.this, -50);
    }

    public GenePane$ScrollRight() {
        super(GenePane.this, "Right", GenePane.access$400());
        mnemonic = 'f';
        toolTipText = "Scroll Right";
    }
}