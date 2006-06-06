// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:05 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   GenePane.java

import javax.swing.AbstractAction;
import javax.swing.ImageIcon;

abstract class GenePane$Action extends AbstractAction {

    protected char mnemonic;
    protected String toolTipText;

    char getMnemonic() {
        return mnemonic;
    }

    String getToolTipText() {
        return toolTipText;
    }

    public GenePane$Action(String s, ImageIcon imageicon) {
        super(s, imageicon);
    }
}