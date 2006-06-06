// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:05 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   GenePane.java

import java.awt.*;
import javax.swing.JPanel;

class GenePane$rulerPane extends JPanel {

    static final int PREFERRED_HEIGHT = 18;
    final GradientPaint rulerPaint = new GradientPaint(0.0F, 0.0F, new Color(238, 238, 238), 0.0F, 10F, new Color(128, 128, 128));
    final Font font = new Font("Helvetica", 0, 8);

    public void paintComponent(Graphics g) {
        Graphics2D graphics2d = (Graphics2D)g;
        int i = GenePane.access$1600(GenePane.this);
        if(i == 0)
            return;
        int j = GenePane.access$1700(GenePane.this);
        if(getStrand() == -1)
            j = directPos2complementPos(GenePane.access$1700(GenePane.this));
        int k = GenePane.access$900(GenePane.this).getWidth();
        graphics2d.setPaint(rulerPaint);
        graphics2d.fillRect(0, 0, getWidth(), getHeight());
        g.setColor(Color.black);
        g.drawString(GenePane.access$1300(GenePane.this).getChnoAndStrandString(), 2, 11);
        if(!GenePane.access$1000(GenePane.this)) {
            g.setFont(font);
            int l = i * 50;
            int i1 = (1 + j / l) * l;
            int j1 = (i1 - j) / i;
            FontMetrics fontmetrics = g.getFontMetrics();
            for(; j1 < GenePane.access$1800(GenePane.this); j1 += 50) {
                g.drawLine(j1 + k, 0, j1 + k, 4);
                for(int k1 = 1; k1 < 5; k1++)
                    g.drawLine(j1 + k + k1 * 10, 0, j1 + k + k1 * 10, 2);

                String s = String.valueOf(i1);
                g.drawString(s, (j1 + k) - fontmetrics.stringWidth(s) / 2, 12);
                if(getStrand() == 1) {
                    i1 += l;
                    continue;
                }
                if(getStrand() == -1)
                    i1 -= l;
            }

        }
    }

    GenePane$rulerPane() {
        super();
        setPreferredSize(new Dimension(0, 18));
        setMaximumSize(new Dimension(1200, 18));
    }
}