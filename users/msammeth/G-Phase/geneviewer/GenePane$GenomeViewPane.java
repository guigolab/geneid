// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:05 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   GenePane.java

import java.awt.*;
import java.awt.event.*;
import java.awt.geom.Rectangle2D;
import java.util.*;
import javax.swing.*;

class GenePane$GenomeViewPane extends JPanel
    implements MouseListener, MouseMotionListener, KeyListener {

    private static final int REGION_INTERVAL = 15;
    private static final int EXON_SHAPE_HEIGHT = 15;
    private static final int INTRON_SHAPE_HEIGHT = 4;
    private static final int REGION_HEIGHT = 30;
    private Alternative noticedAlternative;

    public int genomePosition2Pixel(int i) {
        int j = GenePane.access$1600(GenePane.this);
        int k = GenePane.access$1700(GenePane.this);
        int l = GenePane.access$1700(GenePane.this) + GenePane.access$1800(GenePane.this) * j;
        if(i < k)
            i = k;
        else
        if(i > l)
            i = l;
        return (i - GenePane.access$1700(GenePane.this)) / j;
    }

    private GeneShape getClickedGeneShape(Point point) {
        if(point.y < GenePane.access$1900(GenePane.this) || point.y > 30 * GenePane.access$1300(GenePane.this).getGeneList().size())
            return null;
        for(Iterator iterator = GenePane.access$2000(GenePane.this).iterator(); iterator.hasNext();) {
            GeneShape geneshape = (GeneShape)iterator.next();
            if(geneshape.getSegmentRect().contains(point))
                return geneshape;
        }

        return null;
    }

    private GffRecord getGffData(MouseEvent mouseevent) {
        GeneShape geneshape = getClickedGeneShape(mouseevent.getPoint());
        if(geneshape != null)
            return geneshape.getGffRecord();
        else
            return null;
    }

    private Exon getExonData(MouseEvent mouseevent) {
        GffRecord gffrecord = getGffData(mouseevent);
        if(gffrecord == null)
            return null;
        int i = gffrecord.getStart();
        int j = gffrecord.getEnd();
        if(getStrand() == -1) {
            i = directPos2complementPos(i);
            j = directPos2complementPos(j);
        }
        ArrayList arraylist = GenePane.access$1300(GenePane.this).getGeneList();
        int k = 0;
        do {
            if(k >= arraylist.size())
                break;
            Gene gene = (Gene)arraylist.get(k);
            if(gene.getCdnaId().equals(gffrecord.getLabel())) {
                Exon aexon[] = gene.getExons();
                for(int l = 0; l < aexon.length; l++)
                    if(aexon[l].getGenomeLeft() == i || aexon[l].getGenomeRight() == j)
                        return aexon[l];

                break;
            }
            k++;
        } while(true);
        return null;
    }

    public void mouseEntered(MouseEvent mouseevent) {
    }

    public void mouseExited(MouseEvent mouseevent) {
    }

    public void mouseMoved(MouseEvent mouseevent) {
        GenePane.access$2102(GenePane.this, getGffData(mouseevent));
        if(GenePane.access$2100(GenePane.this) == null) {
            setCursor(Cursor.getPredefinedCursor(0));
            setToolTipText(null);
        } else {
            setCursor(Cursor.getPredefinedCursor(12));
            GffRecord gffrecord;
            if(getStrand() == -1) {
                int i = GenePane.access$2100(GenePane.this).getStart();
                int j = GenePane.access$2100(GenePane.this).getEnd();
                gffrecord = (GffRecord)GenePane.access$2100(GenePane.this).clone();
                gffrecord.setStart(directPos2complementPos(i));
                gffrecord.setEnd(directPos2complementPos(j));
            } else {
                gffrecord = GenePane.access$2100(GenePane.this);
            }
            String s = gffrecord.toString();
            setToolTipText("Press mouse-right button for details.");
            GenePane.access$2200(GenePane.this).setText(s);
        }
    }

    public void mousePressed(MouseEvent mouseevent) {
        if(SwingUtilities.isRightMouseButton(mouseevent)) {
            showAnnotationFrame(mouseevent);
        } else {
            GenePane.access$2302(GenePane.this, mouseevent.getPoint());
            GeneShape geneshape = getClickedGeneShape(mouseevent.getPoint());
            if(geneshape == null) {
                setCursor(Cursor.getPredefinedCursor(13));
            } else {
                setCursor(Cursor.getPredefinedCursor(8));
                GenePane.access$2402(GenePane.this, (GenePane.access$1900(GenePane.this) * -1 + GenePane.access$2300(GenePane.this).y + 7) / 30);
                GenePane.access$2302(GenePane.this, null);
            }
        }
    }

    public void mouseReleased(MouseEvent mouseevent) {
        if(getCursor().getType() == 8) {
            Point point = mouseevent.getPoint();
            int i = (GenePane.access$1900(GenePane.this) * -1 + point.y + 15 + 7) / 30;
            if(i > GenePane.access$2500(GenePane.this).size())
                i = GenePane.access$2500(GenePane.this).size();
            ArrayList arraylist = GenePane.access$1300(GenePane.this).getGeneList();
            Object obj = arraylist.remove(GenePane.access$2400(GenePane.this));
            if(GenePane.access$2400(GenePane.this) >= i)
                arraylist.add(i, obj);
            else
                arraylist.add(i - 1, obj);
            GenePane.access$2602(GenePane.this, GenePane.access$2102(GenePane.this, null));
            GenePane.access$1200(GenePane.this).paint();
            GenePane.access$900(GenePane.this).nt();
            GenePane.access$1100(GenePane.this);
        }
        setCursor(Cursor.getPredefinedCursor(0));
    }

    public void mouseClicked(MouseEvent mouseevent) {
        if(mouseevent.getClickCount() > 1)
            showAnnotationFrame(mouseevent);
    }

    private void showAnnotationFrame(MouseEvent mouseevent) {
        GeneShape geneshape = getClickedGeneShape(mouseevent.getPoint());
        GenePane.access$2700(GenePane.this).init();
        String s = geneshape.getGffRecord().getLabel();
        Gene gene = getGeneByCdnaId(s);
        GenePane.access$2700(GenePane.this).setNoticedGene(gene);
        GenePane.access$2700(GenePane.this).setHighlightGeneShape(geneshape);
        if(geneshape.getSegmentType() == 1) {
            Exon exon = getExonData(mouseevent);
            if(exon == null)
                return;
            geneshape.setExon(exon);
            if(isAlternativeExonClicked(mouseevent.getX(), mouseevent.getY()))
                GenePane.access$2700(GenePane.this).setHighlightAlternativeExon(true);
            GenePane.access$2700(GenePane.this).setAlternativeExon(findAlternativeExon(gene, exon));
        }
        GenePane.access$2700(GenePane.this).showSplicingSequence();
        GenePane.access$2700(GenePane.this).setVisible(true);
    }

    private String numberSuffix(int i) {
        if(i == 1)
            return "st";
        if(i == 2)
            return "nd";
        if(i == 3)
            return "rd";
        else
            return "th";
    }

    public void mouseDragged(MouseEvent mouseevent) {
        Point point = mouseevent.getPoint();
        if(point == null)
            return;
        if(getCursor().getType() == 8) {
            if(GenePane.access$2300(GenePane.this) != null)
                paintGeneTemporarily();
            int i = (GenePane.access$1900(GenePane.this) + point.y + 7) / 30;
            GenePane.access$2200(GenePane.this).setText("Moving the " + (GenePane.access$2400(GenePane.this) + 1) + numberSuffix(GenePane.access$2400(GenePane.this) + 1) + " row to the " + i + numberSuffix(i) + ".");
        } else {
            if(point.distance(GenePane.access$2300(GenePane.this)) < 5D)
                return;
            if(mouseevent.isShiftDown()) {
                GenePane.access$1500(GenePane.this, GenePane.access$1700(GenePane.this) + point.x * GenePane.access$1600(GenePane.this), GenePane.access$2300(GenePane.this).x - point.x);
            } else {
                GenePane.access$300(GenePane.this, point.x - GenePane.access$2300(GenePane.this).x);
                GenePane.access$2800(GenePane.this, point.y - GenePane.access$2300(GenePane.this).y);
            }
        }
        GenePane.access$2302(GenePane.this, point);
        if(getCursor().getType() == 8)
            paintGeneTemporarily();
    }

    public void keyPressed(KeyEvent keyevent) {
        doKeyPressed(keyevent);
    }

    public void keyTyped(KeyEvent keyevent) {
    }

    public void keyReleased(KeyEvent keyevent) {
    }

    public void paintComponent(Graphics g) {
        GenePane.access$900(GenePane.this).nt();
        super.paintComponent(g);
        GenePane.access$1802(GenePane.this, getWidth());
        GenePane.access$2900(GenePane.this);
        if(GenePane.access$1400(GenePane.this) <= 100)
            GenePane.access$1702(GenePane.this, GenePane.access$1300(GenePane.this).getMinGenomePos() - GenePane.access$1600(GenePane.this) * 5);
        GenePane.access$3000(GenePane.this).nt();
        GenePane.access$3100(GenePane.this).nt();
        GenePane.access$1200(GenePane.this).awView();
        Graphics2D graphics2d = (Graphics2D)g;
        paintAlternativeRegion(graphics2d);
        GenePane.access$2000(GenePane.this).clear();
        GenePane.access$3200(GenePane.this).clear();
        int i = GenePane.access$1900(GenePane.this);
        for(Iterator iterator = GenePane.access$1300(GenePane.this).getGeneList().iterator(); iterator.hasNext();) {
            Gene gene = (Gene)iterator.next();
            ArrayList arraylist = (ArrayList)GenePane.access$2500(GenePane.this).get(gene.getCdnaId());
            paintIntron(graphics2d, i, arraylist, gene.getId());
            paintGene(graphics2d, i, arraylist, false, gene == GenePane.access$3300(GenePane.this), gene.getId());
            i += 30;
        }

        if(noticedAlternative != null)
            highlightAlternativeExons(graphics2d, GenePane.access$1900(GenePane.this));
    }

    private void paintAlternativeRegion(Graphics2D graphics2d) {
        if(noticedAlternative == null)
            return;
        int i = getAlternativeRegionHeight();
        if(i <= Math.abs(GenePane.access$1900(GenePane.this)))
            return;
        int j = noticedAlternative.getGenomeLeft();
        int k = noticedAlternative.getGenomeRight();
        if(GenePane.access$1000(GenePane.this)) {
            int l = GenePane.access$3400(GenePane.this, j);
            int j1 = GenePane.access$3400(GenePane.this, k);
            j = GenePane.access$3600(GenePane.this, l * GenePane.access$3500(GenePane.this));
            k = GenePane.access$3600(GenePane.this, j1 * GenePane.access$3500(GenePane.this));
        }
        if(getStrand() == -1) {
            j = complementPos2directPos(j);
            k = complementPos2directPos(k);
        }
        int i1 = genomePosition2Pixel(j);
        int k1 = genomePosition2Pixel(k);
        java.awt.geom.Rectangle2D$Float rectangle2d$float = new java.awt.geom.(i1, GenePane.access$1900(GenePane.this) != 20 ? 0.0F : GenePane.access$1900(GenePane.this) - 1, k1 - i1, i - Math.abs(GenePane.access$1900(GenePane.this)));
        graphics2d.setPaint(GenePane.ALTERNATIVERANGE_COLOR);
        graphics2d.fill(rectangle2d$float);
    }

    private void paintIntron(Graphics2D graphics2d, int i, ArrayList arraylist, int j) {
        GffRecord gffrecord = (GffRecord)arraylist.get(0);
        GffRecord gffrecord1 = (GffRecord)arraylist.get(arraylist.size() - 1);
        int k = -1;
        int l = -1;
        if(!GenePane.access$1000(GenePane.this)) {
            k = gffrecord.getEnd();
            l = gffrecord1.getStart();
        } else {
            int i1 = GenePane.access$3700(GenePane.this, j, 0);
            int j1 = GenePane.access$3800(GenePane.this, j, arraylist.size() - 1);
            k = GenePane.access$3600(GenePane.this, i1 * GenePane.access$3500(GenePane.this));
            l = GenePane.access$3600(GenePane.this, j1 * GenePane.access$3500(GenePane.this));
            if(getStrand() == -1) {
                k = complementPos2directPos(k);
                l = complementPos2directPos(l);
            }
        }
        paintIntron(graphics2d, k + 1, l - 1, i);
    }

    private void paintIntron(Graphics2D graphics2d, int i, int j, int k) {
        int l = genomePosition2Pixel(i);
        int i1 = genomePosition2Pixel(j);
        int j1 = calcIntronYpos(k);
        Rectangle rectangle = new Rectangle(l, j1, (i1 - l) + 1, 4);
        GradientPaint gradientpaint = new GradientPaint(l, j1, GenePane.INTRONCOLOR1, l, j1 + 4, GenePane.INTRONCOLOR2);
        graphics2d.setPaint(gradientpaint);
        graphics2d.fill(rectangle);
    }

    private void paintGeneTemporarily() {
        Graphics2D graphics2d = (Graphics2D)getGraphics();
        graphics2d.setXORMode(Color.white);
        String s = ((Gene)GenePane.access$1300(GenePane.this).getGeneList().get(GenePane.access$2400(GenePane.this))).getCdnaId();
        ArrayList arraylist = (ArrayList)GenePane.access$2500(GenePane.this).get(s);
        paintGene(graphics2d, GenePane.access$2300(GenePane.this).y, arraylist, true, false, ((Gene)GenePane.access$1300(GenePane.this).getGeneList().get(GenePane.access$2400(GenePane.this))).getId());
        graphics2d.setPaintMode();
    }

    private void paintGene(Graphics2D graphics2d, Gene gene, int i) {
        ActualDirectDrawer actualdirectdrawer = new ActualDirectDrawer();
        actualdirectdrawer.setGraphics2D(graphics2d);
        actualdirectdrawer.setVOffset(i);
        Exon aexon[] = gene.getExons();
        for(int j = 0; j < aexon.length; j++)
            aexon[j].acceptDrawer(actualdirectdrawer);

    }

    private void paintGene(Graphics2D graphics2d, int i, ArrayList arraylist, boolean flag, boolean flag1, int j) {
        GffRecord gffrecord = null;
        int k = GenePane.access$1600(GenePane.this);
        int l = calcIntronYpos(i);
        int i1 = GenePane.access$1300(GenePane.this).getMinGenomePos();
        int k1 = -1;
        for(int l1 = 0; l1 < arraylist.size(); l1++) {
            GffRecord gffrecord1 = (GffRecord)arraylist.get(l1);
            int i2 = -1;
            int j2 = -1;
            if(!GenePane.access$1000(GenePane.this)) {
                i2 = gffrecord1.getStart();
                j2 = gffrecord1.getEnd();
            } else {
                int k2 = GenePane.access$3800(GenePane.this, j, l1);
                int k3 = GenePane.access$3700(GenePane.this, j, l1);
                i2 = GenePane.access$3600(GenePane.this, k2 * GenePane.access$3500(GenePane.this));
                j2 = GenePane.access$3600(GenePane.this, k3 * GenePane.access$3500(GenePane.this));
                if(getStrand() == -1) {
                    i2 = complementPos2directPos(i2);
                    j2 = complementPos2directPos(j2);
                }
            }
            if(gffrecord != null)
                i1 = gffrecord.getEnd() + 1;
            int j1 = gffrecord1.getStart() - 1;
            if(j2 < GenePane.access$1700(GenePane.this)) {
                i1 = gffrecord1.getEnd() + 1;
                continue;
            }
            if(i2 > GenePane.access$1700(GenePane.this) + GenePane.access$1800(GenePane.this) * k) {
                if(i1 < GenePane.access$1700(GenePane.this) + GenePane.access$1800(GenePane.this) * k) {
                    int l2 = genomePosition2Pixel(i1);
                    int l3 = genomePosition2Pixel(j1);
                    GeneShape geneshape = new GeneShape();
                    geneshape.setSegmentType(2);
                    geneshape.setSegmentRect(new Rectangle(l2, l, (l3 - l2) + 1, 4));
                    geneshape.setGffRecord(getIntronGffRecord(i1, j1, gffrecord1.getLabel()));
                    GenePane.access$2000(GenePane.this).add(geneshape);
                }
                break;
            }
            paintGff(graphics2d, gffrecord1, i, flag, j, l1);
            if(gffrecord == null) {
                if(l1 != 0) {
                    int i3 = genomePosition2Pixel(i2);
                    if(i3 > 0) {
                        j1 = i2 - 1;
                        int i4 = 0;
                        int k4 = genomePosition2Pixel(i2 - 1);
                        GeneShape geneshape2 = new GeneShape();
                        geneshape2.setSegmentRect(new Rectangle(i4, l, k4 - 1, 4));
                        geneshape2.setSegmentType(2);
                        geneshape2.setGffRecord(getIntronGffRecord(i1, j1, gffrecord1.getLabel()));
                        GenePane.access$2000(GenePane.this).add(geneshape2);
                    }
                }
            } else {
                int j3 = -1;
                int j4 = -1;
                if(!GenePane.access$1000(GenePane.this)) {
                    j3 = genomePosition2Pixel(gffrecord.getEnd() + 1);
                    j4 = genomePosition2Pixel(gffrecord1.getStart() - 1);
                } else {
                    j3 = genomePosition2Pixel(k1 + 1);
                    j4 = genomePosition2Pixel(i2 - 1);
                }
                GeneShape geneshape1 = new GeneShape();
                geneshape1.setSegmentRect(new Rectangle(j3, l, (j4 - j3) + 1, 4));
                geneshape1.setSegmentType(2);
                geneshape1.setGffRecord(getIntronGffRecord(i1, j1, gffrecord1.getLabel()));
                GenePane.access$2000(GenePane.this).add(geneshape1);
            }
            gffrecord = gffrecord1;
            k1 = j2;
        }

    }

    private void paintGff(Graphics2D graphics2d, GffRecord gffrecord, int i, boolean flag, int j, int k) {
        int l = -1;
        int i1 = -1;
        int j1 = -1;
        int k1 = -1;
        if(!GenePane.access$1000(GenePane.this)) {
            l = gffrecord.getStart();
            i1 = gffrecord.getEnd();
            if(getStrand() == -1) {
                l = directPos2complementPos(l);
                i1 = directPos2complementPos(i1);
            }
        } else {
            int l1 = GenePane.access$3800(GenePane.this, j, k);
            int i2 = GenePane.access$3700(GenePane.this, j, k);
            j1 = l = GenePane.access$3600(GenePane.this, l1 * GenePane.access$3500(GenePane.this));
            k1 = i1 = GenePane.access$3600(GenePane.this, i2 * GenePane.access$3500(GenePane.this));
        }
        GeneShape geneshape = new GeneShape();
        Rectangle rectangle = calcExonPixelRect(l, i1, i);
        geneshape.setSegmentRect(rectangle);
        geneshape.setSegmentType(1);
        geneshape.setGffRecord(gffrecord);
        TreeMap treemap = gffrecord.getAttribute();
        String s = (String)treemap.get("ORFLEFT");
        String s1 = (String)treemap.get("ORFRIGHT");
        String s2 = (String)treemap.get("UTRLEFT");
        String s3 = (String)treemap.get("UTRRIGHT");
        int j2 = s != null ? Integer.parseInt(s) : -1;
        int k2 = s1 != null ? Integer.parseInt(s1) : -1;
        int l2 = s2 != null ? Integer.parseInt(s2) : -1;
        int i3 = s3 != null ? Integer.parseInt(s3) : -1;
        if(!GenePane.access$1000(GenePane.this)) {
            if(s != null && s1 != null && s2 != null && s3 != null) {
                if(getStrand() == 1 && j2 < l2 || getStrand() == -1 && j2 > l2) {
                    paintExon(graphics2d, j2, k2, i, GenePane.ORFCOLOR1, GenePane.ORFCOLOR2, flag);
                    paintExon(graphics2d, l2, i3, i, GenePane.UTRCOLOR1, GenePane.UTRCOLOR2, flag);
                } else {
                    paintExon(graphics2d, l2, i3, i, GenePane.UTRCOLOR1, GenePane.UTRCOLOR2, flag);
                    paintExon(graphics2d, j2, k2, i, GenePane.ORFCOLOR1, GenePane.ORFCOLOR2, flag);
                }
            } else
            if(s != null && s1 != null)
                paintExon(graphics2d, j2, k2, i, GenePane.ORFCOLOR1, GenePane.ORFCOLOR2, flag);
            else
            if(s2 != null && s3 != null)
                paintExon(graphics2d, l2, i3, i, GenePane.UTRCOLOR1, GenePane.UTRCOLOR2, flag);
        } else
        if(s != null && s1 != null && s2 != null && s3 != null) {
            int j3 = Math.abs(j2 - k2) + 1;
            int k3 = Math.abs(l2 - i3) + 1;
            if(getStrand() == 1 && j2 < l2 || getStrand() == -1 && j2 > l2) {
                float f = (float)j3 / (float)(j3 + k3);
                int l3 = j1 + Math.round((float)((k1 - j1) + 1) * f);
                paintExon(graphics2d, j1, l3, i, GenePane.ORFCOLOR1, GenePane.ORFCOLOR2, flag);
                paintExon(graphics2d, l3 + 1, k1, i, GenePane.UTRCOLOR1, GenePane.UTRCOLOR2, flag);
            } else {
                float f1 = (float)k3 / (float)(j3 + k3);
                int i4 = j1 + Math.round((float)((k1 - j1) + 1) * f1);
                paintExon(graphics2d, j1, i4, i, GenePane.UTRCOLOR1, GenePane.UTRCOLOR2, flag);
                paintExon(graphics2d, i4 + 1, k1, i, GenePane.ORFCOLOR1, GenePane.ORFCOLOR2, flag);
            }
        } else
        if(s != null && s1 != null)
            paintExon(graphics2d, j1, k1, i, GenePane.ORFCOLOR1, GenePane.ORFCOLOR2, flag);
        else
        if(s2 != null && s3 != null)
            paintExon(graphics2d, j1, k1, i, GenePane.UTRCOLOR1, GenePane.UTRCOLOR2, flag);
        GenePane.access$2000(GenePane.this).add(geneshape);
    }

    private void paintExon(Graphics2D graphics2d, int i, int j, int k, Color color, Color color1, boolean flag) {
        Rectangle rectangle = calcExonPixelRect(i, j, k);
        GradientPaint gradientpaint = new GradientPaint(rectangle.x, rectangle.y, color, rectangle.x, (rectangle.y + rectangle.height) - 1, color1);
        paintExon(graphics2d, rectangle, ((Paint) (gradientpaint)), flag);
    }

    private void paintExon(Graphics2D graphics2d, Rectangle rectangle, Paint paint, boolean flag) {
        graphics2d.setPaint(paint);
        if(flag)
            graphics2d.draw(rectangle);
        else
            graphics2d.fill(rectangle);
    }

    private void highlightAlternativeExons(Graphics2D graphics2d, int i) {
        for(Iterator iterator = noticedAlternative.getExonList().iterator(); iterator.hasNext(); highlightAlternativeExon(graphics2d, i, (Exon)iterator.next()));
    }

    private void highlightAlternativeExon(Graphics2D graphics2d, int i, Exon exon) {
        int ai[] = getAlternativeGeneIds(exon);
        if(ai != null) {
            int l = ai[0];
            int j;
            int k;
            if(!GenePane.access$1000(GenePane.this)) {
                j = exon.getGenomeLeft();
                k = exon.getGenomeRight();
            } else {
                Range range = getAlternativeExonColumnRange(exon);
                int i1 = range.getLeft();
                int k1 = range.getRight();
                j = GenePane.access$3600(GenePane.this, i1 * GenePane.access$3500(GenePane.this));
                k = GenePane.access$3600(GenePane.this, k1 * GenePane.access$3500(GenePane.this));
            }
            Rectangle rectangle = calcExonPixelRect(j, k, i);
            graphics2d.setColor(GenePane.HIGHLIGHTCOLOR1);
            int j1 = rectangle.x;
            int l1 = rectangle.y;
            int i2 = rectangle.width;
            int j2 = rectangle.height;
            int ai1[] = noticedAlternative.getGeneIds1();
            if(ai1[0] != ai[0])
                l1 += ai1.length * 30;
            for(int k2 = 0; k2 < ai.length; k2++) {
                graphics2d.fill(new Rectangle(j1, l1, i2, 2));
                graphics2d.fill(new Rectangle(j1, (l1 + j2) - 2, i2, 2));
                graphics2d.fill(new Rectangle(j1, l1, 2, j2));
                graphics2d.fill(new Rectangle((j1 + i2) - 2, l1, 2, j2));
                GenePane.access$3200(GenePane.this).add(new Rectangle(j1, l1, i2, j2));
                l1 += 30;
            }

        }
    }

    private int[] getAlternativeGeneIds(Exon exon) {
        int ai[] = null;
        Gene gene = (Gene)GenePane.access$1300(GenePane.this).getGeneList().get(0);
        if(gene.hasExon(exon.getGenomeLeft(), exon.getGenomeRight(), false)) {
            ai = noticedAlternative.getGeneIds1();
        } else {
            int ai1[] = noticedAlternative.getGeneIds1();
            Gene gene1 = (Gene)GenePane.access$1300(GenePane.this).getGeneList().get(ai1.length);
            if(gene1.hasExon(exon.getGenomeLeft(), exon.getGenomeRight(), false))
                ai = noticedAlternative.getGeneIds2();
        }
        return ai;
    }

    private Range getAlternativeExonColumnRange(Exon exon) {
        Range range = new Range();
        int i = exon.getGenomeLeft();
        int j = exon.getGenomeRight();
        for(Iterator iterator = GenePane.access$1300(GenePane.this).getGeneList().iterator(); iterator.hasNext();) {
            Gene gene = (Gene)iterator.next();
            Exon exon1 = gene.getExon(i, j);
            if(exon1 != null) {
                range.setLeft(GenePane.access$3900(GenePane.this, gene.getId(), exon1));
                range.setRight(GenePane.access$4000(GenePane.this, gene.getId(), exon1));
                return range;
            }
        }

        int k = GenePane.access$3400(GenePane.this, i);
        int l = GenePane.access$3400(GenePane.this, j);
        if(k == l)
            if(k == 0)
                l = 1;
            else
                k--;
        range.setLeft(k);
        range.setRight(l);
        return range;
    }

    private Rectangle calcExonPixelRect(int i, int j, int k) {
        int l = 0;
        int i1 = 0;
        if(getStrand() == 1) {
            l = genomePosition2Pixel(i);
            i1 = genomePosition2Pixel(j);
        } else
        if(getStrand() == -1) {
            l = genomePosition2Pixel(complementPos2directPos(i));
            i1 = genomePosition2Pixel(complementPos2directPos(j));
        }
        return new Rectangle(l, k, (i1 - l) + 1, 15);
    }

    private int calcIntronYpos(int i) {
        return i + 5;
    }

    private GffRecord getIntronGffRecord(int i, int j, String s) {
        GffRecord gffrecord = new GffRecord();
        gffrecord.setFeature("intron");
        gffrecord.setSource("ALN");
        gffrecord.setStart(i);
        gffrecord.setEnd(j);
        gffrecord.setLabel(s);
        return gffrecord;
    }

    private int getAlternativeRegionHeight() {
        int i = 0;
        if(noticedAlternative != null)
            i = noticedAlternative.getNumAlternativeGenes() * 30 - 7;
        return i;
    }

    private Exon findAlternativeExon(Gene gene, Exon exon) {
        Exon exon1;
label0:
        {
            exon1 = null;
            if(noticedAlternative == null || !noticedAlternative.isAlternativeGene(gene.getId()))
                break label0;
            Iterator iterator = noticedAlternative.getExonList().iterator();
            Exon exon2;
            Range range;
            Range range1;
            do {
                if(!iterator.hasNext())
                    break label0;
                exon2 = (Exon)iterator.next();
                range = new Range(exon2.getGenomeLeft(), exon2.getGenomeRight());
                range1 = new Range(exon.getGenomeLeft(), exon.getGenomeRight());
            } while(!range1.contains(range));
            exon1 = exon2;
        }
        return exon1;
    }

    private boolean isAlternativeExonClicked(int i, int j) {
        boolean flag = false;
        Iterator iterator = GenePane.access$3200(GenePane.this).iterator();
        do {
            if(!iterator.hasNext())
                break;
            Rectangle rectangle = (Rectangle)iterator.next();
            if(rectangle.contains(i, j))
                flag = true;
        } while(true);
        return flag;
    }

    public void doLayout() {
        GenePane.access$4102(GenePane.this, Math.max(GenePane.access$2500(GenePane.this).size() * 30 + 15, getHeight()));
        GenePane.access$1802(GenePane.this, getWidth());
        setPreferredSize(new Dimension(GenePane.access$1800(GenePane.this), GenePane.access$4100(GenePane.this)));
        super.doLayout();
    }

    public Color getORFColor1() {
        return GenePane.ORFCOLOR1;
    }

    public Color getORFColor2() {
        return GenePane.ORFCOLOR2;
    }

    public Color getUTRColor1() {
        return GenePane.UTRCOLOR1;
    }

    public Color getUTRColor2() {
        return GenePane.UTRCOLOR2;
    }

    public Color getHighlightColor1() {
        return GenePane.HIGHLIGHTCOLOR1;
    }

    public Color getHighlightColor2() {
        return GenePane.HIGHLIGHTCOLOR2;
    }

    public Alternative getNoticedAlternative() {
        return noticedAlternative;
    }

    public void setNoticedAlternative(Alternative alternative) {
        noticedAlternative = alternative;
    }

    GenePane$GenomeViewPane() {
        super();
        noticedAlternative = null;
        addKeyListener(this);
        addMouseListener(this);
        addMouseMotionListener(this);
        setBackground(GenePane.BGCOLOR);
        setPreferredSize(new Dimension(600, 300));
    }
}