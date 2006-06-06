// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:05 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   GenePane.java

import java.awt.*;
import java.awt.event.*;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.HashMap;
import javax.swing.JPanel;

class GenePane$overviewPane extends JPanel
    implements MouseListener, MouseMotionListener {

    float scale;
    java.awt.geom. rect;
    Point mouseFrom;

    void init() {
        addMouseListener(this);
        addMouseMotionListener(this);
        scale = 200F / (float)GenePane.access$1300(GenePane.this).getGenomeLen();
        float f = (float)(GenePane.access$1300(GenePane.this).getGeneList().size() + 1) * 5F;
        setPreferredSize(new Dimension(220, (int)f));
    }

    void drawView() {
        int i = (int)(((float)(GenePane.access$1900(GenePane.this) * -1 + 7) * 5F) / 30F + 3F);
        int j = (int)((float)(GenePane.access$1700(GenePane.this) - GenePane.access$1300(GenePane.this).getMinGenomePos()) * scale + 10F);
        drawView(j, i);
    }

    void drawView(float f, float f1) {
        float f2 = ((float)GenePane.access$700(GenePane.this).getHeight() * 5F) / 30F;
        float f3 = (float)(GenePane.access$700(GenePane.this).getWidth() * GenePane.access$1600(GenePane.this)) * scale;
        Graphics2D graphics2d = (Graphics2D)getGraphics();
        graphics2d.setXORMode(Color.white);
        graphics2d.setColor(Color.yellow);
        graphics2d.fill(rect);
        rect.ect(f, f1, f3, f2);
        graphics2d.fill(rect);
        graphics2d.setPaintMode();
    }

    public void paintComponent(Graphics g) {
        super.paintComponent(g);
        Graphics2D graphics2d = (Graphics2D)g;
        ArrayList arraylist = GenePane.access$1300(GenePane.this).getGeneList();
        for(int i = 0; i < arraylist.size(); i++) {
            Gene gene = (Gene)arraylist.get(i);
            ArrayList arraylist1 = (ArrayList)GenePane.access$2500(GenePane.this).get(gene.getCdnaId());
            for(int j = 0; j < arraylist1.size(); j++) {
                GffRecord gffrecord = (GffRecord)arraylist1.get(j);
                int k = -1;
                int l = -1;
                if(!GenePane.access$1000(GenePane.this)) {
                    k = gffrecord.getStart() - GenePane.access$1300(GenePane.this).getMinGenomePos();
                    l = gffrecord.getEnd() - GenePane.access$1300(GenePane.this).getMinGenomePos();
                } else {
                    int i1 = GenePane.access$3800(GenePane.this, gene.getId(), j);
                    int j1 = GenePane.access$3700(GenePane.this, gene.getId(), j);
                    k = GenePane.access$3600(GenePane.this, i1 * GenePane.access$3500(GenePane.this));
                    l = GenePane.access$3600(GenePane.this, j1 * GenePane.access$3500(GenePane.this));
                    if(getStrand() == -1) {
                        k = complementPos2directPos(k);
                        l = complementPos2directPos(l);
                    }
                    int k1 = GenePane.access$3600(GenePane.this, 0 * GenePane.access$3500(GenePane.this));
                    if(getStrand() == -1)
                        k1 = complementPos2directPos(k1);
                    k -= k1;
                    l -= k1;
                }
                if(gene == GenePane.access$3300(GenePane.this))
                    graphics2d.setColor(GenePane.access$700(GenePane.this).getHighlightColor1());
                else
                if(gffrecord.getFeature().equals("ORF"))
                    graphics2d.setColor(GenePane.access$700(GenePane.this).getORFColor1());
                else
                if(gffrecord.getFeature().equals("UTR"))
                    graphics2d.setColor(GenePane.access$700(GenePane.this).getUTRColor1());
                else
                    graphics2d.setColor(Color.BLACK);
                graphics2d.drawLine((int)((float)k * scale) + 10, (i + 1) * 5, (int)((float)l * scale) + 10, (i + 1) * 5);
            }

        }

        graphics2d.setXORMode(Color.white);
        graphics2d.setColor(Color.yellow);
        graphics2d.fill(rect);
        graphics2d.setPaintMode();
    }

    public void mouseClicked(MouseEvent mouseevent) {
    }

    public void mouseMoved(MouseEvent mouseevent) {
        if(rect == null)
            return;
        Point point = mouseevent.getPoint();
        if(rect.ains(point))
            setCursor(Cursor.getPredefinedCursor(13));
        else
            setCursor(Cursor.getPredefinedCursor(0));
    }

    public void mousePressed(MouseEvent mouseevent) {
        if(getCursor().getType() == 13)
            mouseFrom = mouseevent.getPoint();
        else
            mouseFrom = null;
    }

    public void mouseDragged(MouseEvent mouseevent) {
        if(mouseFrom == null)
            return;
        Point point = mouseevent.getPoint();
        int i = (point.x - mouseFrom.x) + (int)rect.();
        int j = (point.y - mouseFrom.y) + (int)rect.();
        if(i < 0)
            i = 0;
        if(j < 0)
            j = 0;
        if((double)i > (double)getWidth() - rect.idth())
            i = getWidth() - (int)rect.idth();
        if((double)j > (double)getHeight() - rect.eight())
            j = getHeight() - (int)rect.eight();
        drawView(i, j);
        mouseFrom = point;
    }

    public void mouseReleased(MouseEvent mouseevent) {
        if(rect == null) {
            return;
        } else {
            mouseFrom = null;
            int i = (int)rect.();
            int j = (int)rect.();
            GenePane.access$1902(GenePane.this, (int)((float)((-j + 3) * 30) / 5F + 7F));
            GenePane.access$1702(GenePane.this, (int)((float)(i - 10) / scale) + GenePane.access$1300(GenePane.this).getMinGenomePos());
            GenePane.access$2800(GenePane.this, 0);
            GenePane.access$300(GenePane.this, 0);
            return;
        }
    }

    public void mouseEntered(MouseEvent mouseevent) {
    }

    public void mouseExited(MouseEvent mouseevent) {
    }

    GenePane$overviewPane() {
        super();
        scale = 0.0F;
        rect = new java.awt.geom.Rectangle2D$Float();
        mouseFrom = null;
    }
}