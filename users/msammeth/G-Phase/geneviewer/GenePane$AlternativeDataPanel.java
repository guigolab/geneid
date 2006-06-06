// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:05 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   GenePane.java

import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.PrintStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import javax.swing.*;

class GenePane$AlternativeDataPanel extends JPanel
    implements MouseListener {

    private static final String EXON_SEPARATOR = "/ ";
    private JLabel patternImageLabel;
    private JLabel alternativeTypeLabel;
    private int alternativeId;

    private void addComponent(Component component, int i, int j, int k, int l) {
        GridBagConstraints gridbagconstraints = new GridBagConstraints();
        gridbagconstraints.fill = 0;
        gridbagconstraints.gridx = i;
        gridbagconstraints.gridy = j;
        gridbagconstraints.gridwidth = k;
        gridbagconstraints.gridheight = l;
        gridbagconstraints.anchor = 17;
        gridbagconstraints.insets = new Insets(1, 1, 1, 1);
        add(component, gridbagconstraints);
    }

    private JLabel createPatternImageLabel(Alternative alternative) {
        JLabel jlabel = null;
        try {
            URL url = new URL(GenePane.access$4200(GenePane.this), alternative.getPatternImageFileName());
            jlabel = new JLabel(new ImageIcon(url));
        }
        catch(MalformedURLException malformedurlexception) {
            System.err.println(malformedurlexception);
        }
        return jlabel;
    }

    private JLabel createNagnagImageLabel() {
        JLabel jlabel = null;
        try {
            URL url = new URL(GenePane.access$4200(GenePane.this), "images/nagnag.jpg");
            jlabel = new JLabel(new ImageIcon(url));
        }
        catch(MalformedURLException malformedurlexception) {
            System.err.println(malformedurlexception);
        }
        return jlabel;
    }

    private JLabel createAlternativeTypeLabel(Alternative alternative) {
        StringBuffer stringbuffer = new StringBuffer(alternative.getCaption());
        if(!alternative.getType().equals("others") && !alternative.getType().equals("5end")) {
            stringbuffer.append(" [");
            stringbuffer.append(alternative.getType());
            stringbuffer.append(']');
        }
        stringbuffer.append(' ');
        stringbuffer.append(alternative.getGenomeLeft());
        stringbuffer.append(",");
        stringbuffer.append(alternative.getGenomeRight());
        stringbuffer.append(" (length ");
        stringbuffer.append(alternative.getAlternaLength());
        stringbuffer.append(')');
        return new JLabel(stringbuffer.toString());
    }

    public int getAlternativeId() {
        return alternativeId;
    }

    public void setAlternativeId(int i) {
        alternativeId = i;
    }

    public void mouseClicked(MouseEvent mouseevent) {
        ArrayList arraylist = GenePane.access$1300(GenePane.this).getSplicingPatternList();
        GenePane.access$800(GenePane.this).esetAlternativeContainerColor();
        Alternative alternative = (Alternative)arraylist.get(alternativeId);
        GenePane.access$1300(GenePane.this).sortGeneListBySplicingPattern(alternativeId);
        int i = -1;
        int j = -1;
        if(!GenePane.access$1000(GenePane.this)) {
            i = alternative.getGenomeLeft();
            j = alternative.getGenomeRight();
            if(getStrand() == -1) {
                i = complementPos2directPos(i);
                j = complementPos2directPos(j);
            }
        } else {
            int k = GenePane.access$3400(GenePane.this, alternative.getGenomeLeft());
            int i1 = GenePane.access$3400(GenePane.this, alternative.getGenomeRight());
            if(k != -1 && i1 == -1)
                i1 = k + 1;
            if(k != -1)
                if(i1 != -1);
            i = GenePane.access$3600(GenePane.this, k * GenePane.access$3500(GenePane.this));
            j = GenePane.access$3600(GenePane.this, i1 * GenePane.access$3500(GenePane.this));
            if(getStrand() == -1) {
                i = complementPos2directPos(i);
                j = complementPos2directPos(j);
            }
        }
        int l = GenePane.access$1300(GenePane.this).getMinGenomePos();
        if(GenePane.access$1700(GenePane.this) > i || j > GenePane.access$1700(GenePane.this) + (GenePane.access$1800(GenePane.this) - 20) * GenePane.access$1600(GenePane.this)) {
            if(l > i || j > l + (GenePane.access$1800(GenePane.this) - 20) * GenePane.access$1600(GenePane.this))
                l = (i + j) / 2 - (GenePane.access$1800(GenePane.this) * GenePane.access$1600(GenePane.this)) / 2;
            if(l < GenePane.access$1300(GenePane.this).getMinGenomePos() - 20 * GenePane.access$1600(GenePane.this))
                l = GenePane.access$1300(GenePane.this).getMinGenomePos() - 20 * GenePane.access$1600(GenePane.this);
            else
            if(l > GenePane.access$1300(GenePane.this).getMaxGenomePos() - (GenePane.access$1800(GenePane.this) - 20) * GenePane.access$1600(GenePane.this))
                l = GenePane.access$1300(GenePane.this).getMaxGenomePos() - (GenePane.access$1800(GenePane.this) - 20) * GenePane.access$1600(GenePane.this);
            GenePane.access$4300(GenePane.this, l);
        }
        GenePane.access$700(GenePane.this).icedAlternative(alternative);
        GenePane.access$700(GenePane.this).t();
        GenePane.access$1200(GenePane.this).repaint();
        setBackground(GenePane.ALTERNATIVERANGE_COLOR);
        GenePane.access$800(GenePane.this).etHighlightedAlternativeId(alternativeId);
    }

    public void mouseEntered(MouseEvent mouseevent) {
    }

    public void mouseExited(MouseEvent mouseevent) {
    }

    public void mousePressed(MouseEvent mouseevent) {
    }

    public void mouseReleased(MouseEvent mouseevent) {
    }

    public GenePane$AlternativeDataPanel(Alternative alternative) {
        super();
        patternImageLabel = null;
        alternativeTypeLabel = null;
        alternativeId = -1;
        addMouseListener(this);
        setLayout(new GridBagLayout());
        setBackground(GenePane.BGCOLOR);
        patternImageLabel = createPatternImageLabel(alternative);
        alternativeTypeLabel = createAlternativeTypeLabel(alternative);
        ArrayList arraylist = alternative.getExonList();
        int i = arraylist.size();
        if(alternative.isNagnag())
            addComponent(createNagnagImageLabel(), 0, 0, 1, 3);
        addComponent(patternImageLabel, 1, 0, 1, 3);
        addComponent(alternativeTypeLabel, 2, 0, i, 1);
        for(int j = 0; j < i; j++) {
            String s = alternative.getExonRangeString(j);
            if(j != 0)
                s = "/ " + s;
            addComponent(new JLabel(s), j + 2, 1, 1, 1);
            s = alternative.getLengthSubtypeString(j);
            if(j != 0)
                s = "/ " + s;
            addComponent(new JLabel(s), j + 2, 2, 1, 1);
        }

    }
}