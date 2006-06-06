// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:06 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   GenePane.java

import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.ArrayList;
import java.util.Iterator;
import javax.swing.JPanel;

class GenePane$SplicingPatternPane extends JPanel
    implements MouseListener {

    static final int LANE_HEIGHT = 20;
    static final int TYPE_COLUMN_WIDTH = 50;
    int noticedAlternativeId;
    JPanel alternativesContainer;
    Component alternativePanels[];

    public void layoutComponents() {
        setBackground(GenePane.BGCOLOR);
        setLayout(new FlowLayout(0));
        alternativesContainer = new JPanel();
        alternativesContainer.setBackground(GenePane.BGCOLOR);
        packAllAlternativePanels();
        add(alternativesContainer);
    }

    public void resetAllAlternativeContainersColor() {
        alternativesContainer.setBackground(GenePane.BGCOLOR);
        for(int i = 0; i < alternativePanels.length; i++)
            resetAlternativeContainerColor(i);

    }

    public void resetAlternativeContainerColor() {
        resetAlternativeContainerColor(noticedAlternativeId);
    }

    public void resetAlternativeContainerColor(int i) {
        if(noticedAlternativeId < 0) {
            return;
        } else {
            alternativePanels[noticedAlternativeId].setBackground(GenePane.BGCOLOR);
            return;
        }
    }

    private void packAllAlternativePanels() {
        alternativePanels = new Component[GenePane.access$1300(GenePane.this).getSplicingPatternList().size()];
        alternativesContainer.setLayout(new GridBagLayout());
        int i = 0;
        for(Iterator iterator = GenePane.access$1300(GenePane.this).getSplicingPatternList().iterator(); iterator.hasNext();) {
            Alternative alternative = (Alternative)iterator.next();
              = new (GenePane.this, alternative);
            .setAlternativeId(i);
            alternativePanels[i] = ;
            packOneAlternativePanel(, i);
            i++;
        }

    }

    private void packOneAlternativePanel(Component component, int i) {
        GridBagConstraints gridbagconstraints = new GridBagConstraints();
        gridbagconstraints.fill = 0;
        gridbagconstraints.gridx = 0;
        gridbagconstraints.gridy = i;
        gridbagconstraints.anchor = 17;
        alternativesContainer.add(component, gridbagconstraints);
    }

    public void setHighlightedAlternativeId(int i) {
        noticedAlternativeId = i;
    }

    public void mouseClicked(MouseEvent mouseevent) {
        resetAllAlternativeContainersColor();
        GenePane.access$700(GenePane.this).ticedAlternative(null);
        GenePane.access$700(GenePane.this).nt();
    }

    public void mouseEntered(MouseEvent mouseevent) {
    }

    public void mouseExited(MouseEvent mouseevent) {
    }

    public void mousePressed(MouseEvent mouseevent) {
    }

    public void mouseReleased(MouseEvent mouseevent) {
    }

    GenePane$SplicingPatternPane() {
        super();
        noticedAlternativeId = -1;
        alternativesContainer = null;
        alternativePanels = null;
        addMouseListener(this);
    }
}