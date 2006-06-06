// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:05 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   GenePane.java

import java.awt.*;
import java.awt.event.*;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import javax.imageio.ImageIO;
import javax.swing.*;

class GenePane$labelPane extends JPanel
    implements MouseMotionListener, MouseListener {

    static final int PREFERRED_WIDTH = 120;
    static final int PREFERRED_HEIGHT = 300;
    BufferedImage nmdImage;
    int highlight;
    ArrayList labels;

    public void paintComponent(Graphics g) {
        super.paintComponent(g);
        Graphics2D graphics2d = (Graphics2D)g;
        Font font = graphics2d.getFont();
        graphics2d.setFont(font.deriveFont(font.getStyle() ^ 1));
        if(nmdImage == null)
            try {
                nmdImage = ImageIO.read(new URL(getBaseURL(), "images/nmd.jpg"));
            }
            catch(IOException ioexception) {
                ioexception.printStackTrace();
            }
        labels.clear();
        ArrayList arraylist = GenePane.access$1300(GenePane.this).getGeneList();
        FontMetrics fontmetrics = g.getFontMetrics();
        int i = fontmetrics.getHeight() - 4;
        int j = GenePane.access$1900(GenePane.this);
        for(int k = 0; k < arraylist.size(); k++) {
            Gene gene = (Gene)arraylist.get(k);
            String s = gene.getCdnaId();
            if(gene.isNmd())
                graphics2d.drawImage(nmdImage, 0, j - 5, this);
            graphics2d.setColor(highlight != k ? Color.black : Color.red);
            graphics2d.drawString(s, 10, i + j);
            Rectangle2D rectangle2d = fontmetrics.getStringBounds(s, graphics2d);
            rectangle2d.setRect(8D, j, rectangle2d.getWidth(), rectangle2d.getHeight());
            if(GenePane.access$3300(GenePane.this) == gene) {
                graphics2d.setColor(Color.yellow);
                graphics2d.draw(rectangle2d);
            }
            labels.add(rectangle2d);
            j += 30;
        }

    }

    public void mouseMoved(MouseEvent mouseevent) {
        java.awt.Point point = mouseevent.getPoint();
        for(int i = 0; i < labels.size(); i++)
            if(((Rectangle2D)labels.get(i)).contains(point)) {
                if(highlight != i) {
                    highlight = i;
                    String s = ((Gene)GenePane.access$1300(GenePane.this).getGeneList().get(i)).getGeneName();
                    setToolTipText("Press mouse-right button for details.");
                    GenePane.access$2200(GenePane.this).setText(s);
                    repaint();
                }
                return;
            }

        if(highlight != -1) {
            highlight = -1;
            repaint();
            setToolTipText(null);
        }
    }

    public void mouseClicked(MouseEvent mouseevent) {
        GenePane.access$3302(GenePane.this, null);
        java.awt.Point point = mouseevent.getPoint();
        for(int i = 0; i < labels.size(); i++)
            if(((Rectangle2D)labels.get(i)).contains(point))
                GenePane.access$3302(GenePane.this, (Gene)GenePane.access$1300(GenePane.this).getGeneList().get(i));

        if(GenePane.access$3300(GenePane.this) == null)
            return;
        if(SwingUtilities.isRightMouseButton(mouseevent))
            showAnnotationFrame();
        int j = mouseevent.getClickCount();
        if(j == 1) {
            repaint();
            GenePane.access$1200(GenePane.this).repaint();
        } else {
            showAnnotationFrame();
        }
    }

    private void showAnnotationFrame() {
        GenePane.access$2700(GenePane.this).init();
        GenePane.access$2700(GenePane.this).setNoticedGene(GenePane.access$3300(GenePane.this));
        GenePane.access$2700(GenePane.this).setVisible(true);
    }

    public void mouseDragged(MouseEvent mouseevent) {
    }

    public void mouseEntered(MouseEvent mouseevent) {
    }

    public void mouseExited(MouseEvent mouseevent) {
    }

    public void mousePressed(MouseEvent mouseevent) {
    }

    public void mouseReleased(MouseEvent mouseevent) {
    }

    GenePane$labelPane() {
        super();
        nmdImage = null;
        highlight = -1;
        labels = new ArrayList();
        setBackground(GenePane.BGCOLOR.darker());
        addMouseListener(this);
        addMouseMotionListener(this);
        setPreferredSize(new Dimension(120, 300));
        setMinimumSize(new Dimension(120, 0));
        setMaximumSize(new Dimension(120, 10000));
    }
}