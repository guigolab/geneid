// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:06 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   GenePane.java

import java.awt.Dimension;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

class GenePane$ZoomSpinner extends JSpinner
    implements ChangeListener {

    public void stateChanged(ChangeEvent changeevent) {
        int i = ((Integer)getValue()).intValue();
        if(i <= 100)
            GenePane.access$1500(GenePane.this, (GenePane.access$1300(GenePane.this).getMinGenomePos() + GenePane.access$1300(GenePane.this).getMaxGenomePos()) / 2, i - GenePane.access$1400(GenePane.this));
        else
            GenePane.access$1500(GenePane.this, -1, i - GenePane.access$1400(GenePane.this));
    }

    public GenePane$ZoomSpinner(SpinnerNumberModel spinnernumbermodel) {
        super(spinnernumbermodel);
        setMaximumSize(new Dimension(50, 20));
        setPreferredSize(new Dimension(50, 20));
        addChangeListener(this);
    }
}