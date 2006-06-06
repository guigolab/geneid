// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:06 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   GenePane.java

import java.applet.AppletContext;
import java.awt.*;
import java.awt.event.*;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.*;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.*;
import javax.imageio.ImageIO;
import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.xml.parsers.*;
import org.xml.sax.SAXException;

public class GenePane extends JPanel {
    class overviewPane extends JPanel
        implements MouseListener, MouseMotionListener {

        float scale;
        java.awt.geom.Rectangle2D.Float rect;
        Point mouseFrom;

        void init() {
            addMouseListener(this);
            addMouseMotionListener(this);
            scale = 200F / (float)locus.getGenomeLen();
            float f = (float)(locus.getGeneList().size() + 1) * 5F;
            setPreferredSize(new Dimension(220, (int)f));
        }

        void drawView() {
            int i = (int)(((float)(vOffset * -1 + 7) * 5F) / 30F + 3F);
            int j = (int)((float)(hOffset - locus.getMinGenomePos()) * scale + 10F);
            drawView(j, i);
        }

        void drawView(float f, float f1) {
            float f2 = ((float)regionPanel.getHeight() * 5F) / 30F;
            float f3 = (float)(regionPanel.getWidth() * getBasePerPixel()) * scale;
            Graphics2D graphics2d = (Graphics2D)getGraphics();
            graphics2d.setXORMode(Color.white);
            graphics2d.setColor(Color.yellow);
            graphics2d.fill(rect);
            rect.setRect(f, f1, f3, f2);
            graphics2d.fill(rect);
            graphics2d.setPaintMode();
        }

        public void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D graphics2d = (Graphics2D)g;
            ArrayList arraylist = locus.getGeneList();
            for(int i = 0; i < arraylist.size(); i++) {
                Gene gene = (Gene)arraylist.get(i);
                ArrayList arraylist1 = (ArrayList)gffRecordsMap.get(gene.getCdnaId());
                for(int j = 0; j < arraylist1.size(); j++) {
                    GffRecord gffrecord = (GffRecord)arraylist1.get(j);
                    int k = -1;
                    int l = -1;
                    if(!equalScaleView) {
                        k = gffrecord.getStart() - locus.getMinGenomePos();
                        l = gffrecord.getEnd() - locus.getMinGenomePos();
                    } else {
                        int i1 = getLeftLinePos(gene.getId(), j);
                        int j1 = getRightLinePos(gene.getId(), j);
                        k = translateGenomePos(i1 * getOneExonLenEqualView());
                        l = translateGenomePos(j1 * getOneExonLenEqualView());
                        if(getStrand() == -1) {
                            k = complementPos2directPos(k);
                            l = complementPos2directPos(l);
                        }
                        int k1 = translateGenomePos(0 * getOneExonLenEqualView());
                        if(getStrand() == -1)
                            k1 = complementPos2directPos(k1);
                        k -= k1;
                        l -= k1;
                    }
                    if(gene == focusedGene)
                        graphics2d.setColor(regionPanel.getHighlightColor1());
                    else
                    if(gffrecord.getFeature().equals("ORF"))
                        graphics2d.setColor(regionPanel.getORFColor1());
                    else
                    if(gffrecord.getFeature().equals("UTR"))
                        graphics2d.setColor(regionPanel.getUTRColor1());
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
            if(rect.contains(point))
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
            int i = (point.x - mouseFrom.x) + (int)rect.getX();
            int j = (point.y - mouseFrom.y) + (int)rect.getY();
            if(i < 0)
                i = 0;
            if(j < 0)
                j = 0;
            if((double)i > (double)getWidth() - rect.getWidth())
                i = getWidth() - (int)rect.getWidth();
            if((double)j > (double)getHeight() - rect.getHeight())
                j = getHeight() - (int)rect.getHeight();
            drawView(i, j);
            mouseFrom = point;
        }

        public void mouseReleased(MouseEvent mouseevent) {
            if(rect == null) {
                return;
            } else {
                mouseFrom = null;
                int i = (int)rect.getX();
                int j = (int)rect.getY();
                vOffset = (int)((float)((-j + 3) * 30) / 5F + 7F);
                hOffset = (int)((float)(i - 10) / scale) + locus.getMinGenomePos();
                vScroll(0);
                hScroll(0);
                return;
            }
        }

        public void mouseEntered(MouseEvent mouseevent) {
        }

        public void mouseExited(MouseEvent mouseevent) {
        }

        overviewPane() {
            super();
            scale = 0.0F;
            rect = new java.awt.geom.Rectangle2D.Float();
            mouseFrom = null;
        }
    }

    class AlternativeDataPanel extends JPanel
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
                URL url = new URL(baseURL, alternative.getPatternImageFileName());
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
                URL url = new URL(baseURL, "images/nagnag.jpg");
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
            ArrayList arraylist = locus.getSplicingPatternList();
            alternativePatternPanel.resetAlternativeContainerColor();
            Alternative alternative = (Alternative)arraylist.get(alternativeId);
            locus.sortGeneListBySplicingPattern(alternativeId);
            int i = -1;
            int j = -1;
            if(!equalScaleView) {
                i = alternative.getGenomeLeft();
                j = alternative.getGenomeRight();
                if(getStrand() == -1) {
                    i = complementPos2directPos(i);
                    j = complementPos2directPos(j);
                }
            } else {
                int k = genomePos2equalViewColumnNo(alternative.getGenomeLeft());
                int i1 = genomePos2equalViewColumnNo(alternative.getGenomeRight());
                if(k != -1 && i1 == -1)
                    i1 = k + 1;
                if(k != -1)
                    if(i1 != -1);
                i = translateGenomePos(k * getOneExonLenEqualView());
                j = translateGenomePos(i1 * getOneExonLenEqualView());
                if(getStrand() == -1) {
                    i = complementPos2directPos(i);
                    j = complementPos2directPos(j);
                }
            }
            int l = locus.getMinGenomePos();
            if(hOffset > i || j > hOffset + (panelWidth - 20) * getBasePerPixel()) {
                if(l > i || j > l + (panelWidth - 20) * getBasePerPixel())
                    l = (i + j) / 2 - (panelWidth * getBasePerPixel()) / 2;
                if(l < locus.getMinGenomePos() - 20 * getBasePerPixel())
                    l = locus.getMinGenomePos() - 20 * getBasePerPixel();
                else
                if(l > locus.getMaxGenomePos() - (panelWidth - 20) * getBasePerPixel())
                    l = locus.getMaxGenomePos() - (panelWidth - 20) * getBasePerPixel();
                setHOffset(l);
            }
            regionPanel.setNoticedAlternative(alternative);
            regionPanel.repaint();
            ovPanel.repaint();
            setBackground(GenePane.ALTERNATIVERANGE_COLOR);
            alternativePatternPanel.setHighlightedAlternativeId(alternativeId);
        }

        public void mouseEntered(MouseEvent mouseevent) {
        }

        public void mouseExited(MouseEvent mouseevent) {
        }

        public void mousePressed(MouseEvent mouseevent) {
        }

        public void mouseReleased(MouseEvent mouseevent) {
        }

        public AlternativeDataPanel(Alternative alternative) {
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

    class SplicingPatternPane extends JPanel
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
            alternativePanels = new Component[locus.getSplicingPatternList().size()];
            alternativesContainer.setLayout(new GridBagLayout());
            int i = 0;
            for(Iterator iterator = locus.getSplicingPatternList().iterator(); iterator.hasNext();) {
                Alternative alternative = (Alternative)iterator.next();
                AlternativeDataPanel alternativedatapanel = new AlternativeDataPanel(alternative);
                alternativedatapanel.setAlternativeId(i);
                alternativePanels[i] = alternativedatapanel;
                packOneAlternativePanel(alternativedatapanel, i);
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
            regionPanel.setNoticedAlternative(null);
            regionPanel.repaint();
        }

        public void mouseEntered(MouseEvent mouseevent) {
        }

        public void mouseExited(MouseEvent mouseevent) {
        }

        public void mousePressed(MouseEvent mouseevent) {
        }

        public void mouseReleased(MouseEvent mouseevent) {
        }

        SplicingPatternPane() {
            super();
            noticedAlternativeId = -1;
            alternativesContainer = null;
            alternativePanels = null;
            addMouseListener(this);
        }
    }

    class labelPane extends JPanel
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
            ArrayList arraylist = locus.getGeneList();
            FontMetrics fontmetrics = g.getFontMetrics();
            int i = fontmetrics.getHeight() - 4;
            int j = vOffset;
            for(int k = 0; k < arraylist.size(); k++) {
                Gene gene = (Gene)arraylist.get(k);
                String s = gene.getCdnaId();
                if(gene.isNmd())
                    graphics2d.drawImage(nmdImage, 0, j - 5, this);
                graphics2d.setColor(highlight != k ? Color.black : Color.red);
                graphics2d.drawString(s, 10, i + j);
                Rectangle2D rectangle2d = fontmetrics.getStringBounds(s, graphics2d);
                rectangle2d.setRect(8D, j, rectangle2d.getWidth(), rectangle2d.getHeight());
                if(focusedGene == gene) {
                    graphics2d.setColor(Color.yellow);
                    graphics2d.draw(rectangle2d);
                }
                labels.add(rectangle2d);
                j += 30;
            }

        }

        public void mouseMoved(MouseEvent mouseevent) {
            Point point = mouseevent.getPoint();
            for(int i = 0; i < labels.size(); i++)
                if(((Rectangle2D)labels.get(i)).contains(point)) {
                    if(highlight != i) {
                        highlight = i;
                        String s = ((Gene)locus.getGeneList().get(i)).getGeneName();
                        setToolTipText("Press mouse-right button for details.");
                        messagePanel.setText(s);
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
            focusedGene = null;
            Point point = mouseevent.getPoint();
            for(int i = 0; i < labels.size(); i++)
                if(((Rectangle2D)labels.get(i)).contains(point))
                    focusedGene = (Gene)locus.getGeneList().get(i);

            if(focusedGene == null)
                return;
            if(SwingUtilities.isRightMouseButton(mouseevent))
                showAnnotationFrame();
            int j = mouseevent.getClickCount();
            if(j == 1) {
                repaint();
                ovPanel.repaint();
            } else {
                showAnnotationFrame();
            }
        }

        private void showAnnotationFrame() {
            annotationFrame.init();
            annotationFrame.setNoticedGene(focusedGene);
            annotationFrame.setVisible(true);
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

        labelPane() {
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

    class rulerPane extends JPanel {

        static final int PREFERRED_HEIGHT = 18;
        final GradientPaint rulerPaint = new GradientPaint(0.0F, 0.0F, new Color(238, 238, 238), 0.0F, 10F, new Color(128, 128, 128));
        final Font font = new Font("Helvetica", 0, 8);

        public void paintComponent(Graphics g) {
            Graphics2D graphics2d = (Graphics2D)g;
            int i = getBasePerPixel();
            if(i == 0)
                return;
            int j = hOffset;
            if(getStrand() == -1)
                j = directPos2complementPos(hOffset);
            int k = labelPanel.getWidth();
            graphics2d.setPaint(rulerPaint);
            graphics2d.fillRect(0, 0, getWidth(), getHeight());
            g.setColor(Color.black);
            g.drawString(locus.getChnoAndStrandString(), 2, 11);
            if(!equalScaleView) {
                g.setFont(font);
                int l = i * 50;
                int i1 = (1 + j / l) * l;
                int j1 = (i1 - j) / i;
                FontMetrics fontmetrics = g.getFontMetrics();
                for(; j1 < panelWidth; j1 += 50) {
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

        rulerPane() {
            super();
            setPreferredSize(new Dimension(0, 18));
            setMaximumSize(new Dimension(1200, 18));
        }
    }

    class GenomeViewPane extends JPanel
        implements MouseListener, MouseMotionListener, KeyListener {

        private static final int REGION_INTERVAL = 15;
        private static final int EXON_SHAPE_HEIGHT = 15;
        private static final int INTRON_SHAPE_HEIGHT = 4;
        private static final int REGION_HEIGHT = 30;
        private Alternative noticedAlternative;

        public int genomePosition2Pixel(int i) {
            int j = getBasePerPixel();
            int k = hOffset;
            int l = hOffset + panelWidth * j;
            if(i < k)
                i = k;
            else
            if(i > l)
                i = l;
            return (i - hOffset) / j;
        }

        private GeneShape getClickedGeneShape(Point point) {
            if(point.y < vOffset || point.y > 30 * locus.getGeneList().size())
                return null;
            for(Iterator iterator = drawnGeneSegmentList.iterator(); iterator.hasNext();) {
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
            ArrayList arraylist = locus.getGeneList();
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
            pointedGff = getGffData(mouseevent);
            if(pointedGff == null) {
                setCursor(Cursor.getPredefinedCursor(0));
                setToolTipText(null);
            } else {
                setCursor(Cursor.getPredefinedCursor(12));
                GffRecord gffrecord;
                if(getStrand() == -1) {
                    int i = pointedGff.getStart();
                    int j = pointedGff.getEnd();
                    gffrecord = (GffRecord)pointedGff.clone();
                    gffrecord.setStart(directPos2complementPos(i));
                    gffrecord.setEnd(directPos2complementPos(j));
                } else {
                    gffrecord = pointedGff;
                }
                String s = gffrecord.toString();
                setToolTipText("Press mouse-right button for details.");
                messagePanel.setText(s);
            }
        }

        public void mousePressed(MouseEvent mouseevent) {
            if(SwingUtilities.isRightMouseButton(mouseevent)) {
                showAnnotationFrame(mouseevent);
            } else {
                mouseFrom = mouseevent.getPoint();
                GeneShape geneshape = getClickedGeneShape(mouseevent.getPoint());
                if(geneshape == null) {
                    setCursor(Cursor.getPredefinedCursor(13));
                } else {
                    setCursor(Cursor.getPredefinedCursor(8));
                    focusedRow = (vOffset * -1 + mouseFrom.y + 7) / 30;
                    mouseFrom = null;
                }
            }
        }

        public void mouseReleased(MouseEvent mouseevent) {
            if(getCursor().getType() == 8) {
                Point point = mouseevent.getPoint();
                int i = (vOffset * -1 + point.y + 15 + 7) / 30;
                if(i > gffRecordsMap.size())
                    i = gffRecordsMap.size();
                ArrayList arraylist = locus.getGeneList();
                Object obj = arraylist.remove(focusedRow);
                if(focusedRow >= i)
                    arraylist.add(i, obj);
                else
                    arraylist.add(i - 1, obj);
                focusedGff = pointedGff = null;
                ovPanel.repaint();
                labelPanel.repaint();
                repaintPanel();
            }
            setCursor(Cursor.getPredefinedCursor(0));
        }

        public void mouseClicked(MouseEvent mouseevent) {
            if(mouseevent.getClickCount() > 1)
                showAnnotationFrame(mouseevent);
        }

        private void showAnnotationFrame(MouseEvent mouseevent) {
            GeneShape geneshape = getClickedGeneShape(mouseevent.getPoint());
            annotationFrame.init();
            String s = geneshape.getGffRecord().getLabel();
            Gene gene = getGeneByCdnaId(s);
            annotationFrame.setNoticedGene(gene);
            annotationFrame.setHighlightGeneShape(geneshape);
            if(geneshape.getSegmentType() == 1) {
                Exon exon = getExonData(mouseevent);
                if(exon == null)
                    return;
                geneshape.setExon(exon);
                if(isAlternativeExonClicked(mouseevent.getX(), mouseevent.getY()))
                    annotationFrame.setHighlightAlternativeExon(true);
                annotationFrame.setAlternativeExon(findAlternativeExon(gene, exon));
            }
            annotationFrame.showSplicingSequence();
            annotationFrame.setVisible(true);
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
                if(mouseFrom != null)
                    paintGeneTemporarily();
                int i = (vOffset + point.y + 7) / 30;
                messagePanel.setText("Moving the " + (focusedRow + 1) + numberSuffix(focusedRow + 1) + " row to the " + i + numberSuffix(i) + ".");
            } else {
                if(point.distance(mouseFrom) < 5D)
                    return;
                if(mouseevent.isShiftDown()) {
                    zoom(hOffset + point.x * getBasePerPixel(), mouseFrom.x - point.x);
                } else {
                    hScroll(point.x - mouseFrom.x);
                    vScroll(point.y - mouseFrom.y);
                }
            }
            mouseFrom = point;
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
            labelPanel.repaint();
            super.paintComponent(g);
            panelWidth = getWidth();
            initZoomScale();
            if(magPower <= 100)
                hOffset = locus.getMinGenomePos() - getBasePerPixel() * 5;
            rulerPanelA.repaint();
            rulerPanelB.repaint();
            ovPanel.drawView();
            Graphics2D graphics2d = (Graphics2D)g;
            paintAlternativeRegion(graphics2d);
            drawnGeneSegmentList.clear();
            alternativeExonPixelRectList.clear();
            int i = vOffset;
            for(Iterator iterator = locus.getGeneList().iterator(); iterator.hasNext();) {
                Gene gene = (Gene)iterator.next();
                ArrayList arraylist = (ArrayList)gffRecordsMap.get(gene.getCdnaId());
                paintIntron(graphics2d, i, arraylist, gene.getId());
                paintGene(graphics2d, i, arraylist, false, gene == focusedGene, gene.getId());
                i += 30;
            }

            if(noticedAlternative != null)
                highlightAlternativeExons(graphics2d, vOffset);
        }

        private void paintAlternativeRegion(Graphics2D graphics2d) {
            if(noticedAlternative == null)
                return;
            int i = getAlternativeRegionHeight();
            if(i <= Math.abs(vOffset))
                return;
            int j = noticedAlternative.getGenomeLeft();
            int k = noticedAlternative.getGenomeRight();
            if(equalScaleView) {
                int l = genomePos2equalViewColumnNo(j);
                int j1 = genomePos2equalViewColumnNo(k);
                j = translateGenomePos(l * getOneExonLenEqualView());
                k = translateGenomePos(j1 * getOneExonLenEqualView());
            }
            if(getStrand() == -1) {
                j = complementPos2directPos(j);
                k = complementPos2directPos(k);
            }
            int i1 = genomePosition2Pixel(j);
            int k1 = genomePosition2Pixel(k);
            java.awt.geom.Rectangle2D.Float float1 = new java.awt.geom.Rectangle2D.Float(i1, vOffset != 20 ? 0.0F : vOffset - 1, k1 - i1, i - Math.abs(vOffset));
            graphics2d.setPaint(GenePane.ALTERNATIVERANGE_COLOR);
            graphics2d.fill(float1);
        }

        private void paintIntron(Graphics2D graphics2d, int i, ArrayList arraylist, int j) {
            GffRecord gffrecord = (GffRecord)arraylist.get(0);
            GffRecord gffrecord1 = (GffRecord)arraylist.get(arraylist.size() - 1);
            int k = -1;
            int l = -1;
            if(!equalScaleView) {
                k = gffrecord.getEnd();
                l = gffrecord1.getStart();
            } else {
                int i1 = getRightLinePos(j, 0);
                int j1 = getLeftLinePos(j, arraylist.size() - 1);
                k = translateGenomePos(i1 * getOneExonLenEqualView());
                l = translateGenomePos(j1 * getOneExonLenEqualView());
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
            String s = ((Gene)locus.getGeneList().get(focusedRow)).getCdnaId();
            ArrayList arraylist = (ArrayList)gffRecordsMap.get(s);
            paintGene(graphics2d, mouseFrom.y, arraylist, true, false, ((Gene)locus.getGeneList().get(focusedRow)).getId());
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
            int k = getBasePerPixel();
            int l = calcIntronYpos(i);
            int i1 = locus.getMinGenomePos();
            int k1 = -1;
            for(int l1 = 0; l1 < arraylist.size(); l1++) {
                GffRecord gffrecord1 = (GffRecord)arraylist.get(l1);
                int i2 = -1;
                int j2 = -1;
                if(!equalScaleView) {
                    i2 = gffrecord1.getStart();
                    j2 = gffrecord1.getEnd();
                } else {
                    int k2 = getLeftLinePos(j, l1);
                    int k3 = getRightLinePos(j, l1);
                    i2 = translateGenomePos(k2 * getOneExonLenEqualView());
                    j2 = translateGenomePos(k3 * getOneExonLenEqualView());
                    if(getStrand() == -1) {
                        i2 = complementPos2directPos(i2);
                        j2 = complementPos2directPos(j2);
                    }
                }
                if(gffrecord != null)
                    i1 = gffrecord.getEnd() + 1;
                int j1 = gffrecord1.getStart() - 1;
                if(j2 < hOffset) {
                    i1 = gffrecord1.getEnd() + 1;
                    continue;
                }
                if(i2 > hOffset + panelWidth * k) {
                    if(i1 < hOffset + panelWidth * k) {
                        int l2 = genomePosition2Pixel(i1);
                        int l3 = genomePosition2Pixel(j1);
                        GeneShape geneshape = new GeneShape();
                        geneshape.setSegmentType(2);
                        geneshape.setSegmentRect(new Rectangle(l2, l, (l3 - l2) + 1, 4));
                        geneshape.setGffRecord(getIntronGffRecord(i1, j1, gffrecord1.getLabel()));
                        drawnGeneSegmentList.add(geneshape);
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
                            drawnGeneSegmentList.add(geneshape2);
                        }
                    }
                } else {
                    int j3 = -1;
                    int j4 = -1;
                    if(!equalScaleView) {
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
                    drawnGeneSegmentList.add(geneshape1);
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
            if(!equalScaleView) {
                l = gffrecord.getStart();
                i1 = gffrecord.getEnd();
                if(getStrand() == -1) {
                    l = directPos2complementPos(l);
                    i1 = directPos2complementPos(i1);
                }
            } else {
                int l1 = getLeftLinePos(j, k);
                int i2 = getRightLinePos(j, k);
                j1 = l = translateGenomePos(l1 * getOneExonLenEqualView());
                k1 = i1 = translateGenomePos(i2 * getOneExonLenEqualView());
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
            if(!equalScaleView) {
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
            drawnGeneSegmentList.add(geneshape);
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
                if(!equalScaleView) {
                    j = exon.getGenomeLeft();
                    k = exon.getGenomeRight();
                } else {
                    Range range = getAlternativeExonColumnRange(exon);
                    int i1 = range.getLeft();
                    int k1 = range.getRight();
                    j = translateGenomePos(i1 * getOneExonLenEqualView());
                    k = translateGenomePos(k1 * getOneExonLenEqualView());
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
                    alternativeExonPixelRectList.add(new Rectangle(j1, l1, i2, j2));
                    l1 += 30;
                }

            }
        }

        private int[] getAlternativeGeneIds(Exon exon) {
            int ai[] = null;
            Gene gene = (Gene)locus.getGeneList().get(0);
            if(gene.hasExon(exon.getGenomeLeft(), exon.getGenomeRight(), false)) {
                ai = noticedAlternative.getGeneIds1();
            } else {
                int ai1[] = noticedAlternative.getGeneIds1();
                Gene gene1 = (Gene)locus.getGeneList().get(ai1.length);
                if(gene1.hasExon(exon.getGenomeLeft(), exon.getGenomeRight(), false))
                    ai = noticedAlternative.getGeneIds2();
            }
            return ai;
        }

        private Range getAlternativeExonColumnRange(Exon exon) {
            Range range = new Range();
            int i = exon.getGenomeLeft();
            int j = exon.getGenomeRight();
            for(Iterator iterator = locus.getGeneList().iterator(); iterator.hasNext();) {
                Gene gene = (Gene)iterator.next();
                Exon exon1 = gene.getExon(i, j);
                if(exon1 != null) {
                    range.setLeft(getLeftLinePos(gene.getId(), exon1));
                    range.setRight(getRightLinePos(gene.getId(), exon1));
                    return range;
                }
            }

            int k = genomePos2equalViewColumnNo(i);
            int l = genomePos2equalViewColumnNo(j);
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
            Iterator iterator = alternativeExonPixelRectList.iterator();
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
            panelHeight = Math.max(gffRecordsMap.size() * 30 + 15, getHeight());
            panelWidth = getWidth();
            setPreferredSize(new Dimension(panelWidth, panelHeight));
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

        GenomeViewPane() {
            super();
            noticedAlternative = null;
            addKeyListener(this);
            addMouseListener(this);
            addMouseMotionListener(this);
            setBackground(GenePane.BGCOLOR);
            setPreferredSize(new Dimension(600, 300));
        }
    }

    class HomoloGene
        implements ActionListener {

        public void actionPerformed(ActionEvent actionevent) {
            showHomoloGene();
        }

        private void showHomoloGene() {
            try {
                HashMap hashmap = locus.getHomoloGeneMap();
                Iterator iterator = hashmap.keySet().iterator();
                do {
                    if(!iterator.hasNext())
                        break;
                    String s = (String)iterator.next();
                    Integer integer = (Integer)hashmap.get(s);
                    if(!s.equals(getSpecies())) {
                        String s1 = "viewer.php?sp=" + s + "&lociId=" + integer.toString();
                        getAppletContext().showDocument(new URL(getBaseURL(), s1), "_blank");
                    }
                } while(true);
            }
            catch(MalformedURLException malformedurlexception) { }
        }

        HomoloGene() {
            super();
        }
    }

    class ZoomSpinner extends JSpinner
        implements ChangeListener {

        public void stateChanged(ChangeEvent changeevent) {
            int i = ((Integer)getValue()).intValue();
            if(i <= 100)
                zoom((locus.getMinGenomePos() + locus.getMaxGenomePos()) / 2, i - magPower);
            else
                zoom(-1, i - magPower);
        }

        public ZoomSpinner(SpinnerNumberModel spinnernumbermodel) {
            super(spinnernumbermodel);
            setMaximumSize(new Dimension(50, 20));
            setPreferredSize(new Dimension(50, 20));
            addChangeListener(this);
        }
    }

    class ChangeView
        implements ActionListener {

        public void actionPerformed(ActionEvent actionevent) {
            equalScaleView = !equalScaleView;
            JButton jbutton = (JButton)actionevent.getSource();
            if(equalScaleView)
                jbutton.setText("Show proportional view");
            else
                jbutton.setText("Show monospaced View");
            repaintPanel();
            ovPanel.repaint();
        }

        ChangeView() {
            super();
        }
    }

    class Option extends Action {

        JFrame optionPane;

        public void actionPerformed(ActionEvent actionevent) {
            showOptionPane();
        }

        void init() {
            JPanel jpanel = new JPanel();
            jpanel.setLayout(new GridLayout(0, 2));
            jpanel.add(new JLabel("ORF Color"));
            JButton jbutton = new JButton("change");
            jbutton.setBackground(GenePane.ORFCOLOR1);
            jbutton.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent actionevent) {
                    Color color = JColorChooser.showDialog(null, "Choose ORF Color", GenePane.ORFCOLOR1);
                    if(color == null) {
                        return;
                    } else {
                        GenePane.ORFCOLOR1 = color;
                        optionPane.setVisible(false);
                        repaint();
                        return;
                    }
                }

                 {
                    super();
                }
            });
            jpanel.add(jbutton);
            jpanel.add(new JLabel("UTR Color"));
            jbutton = new JButton("change");
            jbutton.setBackground(GenePane.UTRCOLOR1);
            jbutton.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent actionevent) {
                    Color color = JColorChooser.showDialog(null, "Choose UTR Color", GenePane.UTRCOLOR1);
                    if(color == null) {
                        return;
                    } else {
                        GenePane.UTRCOLOR1 = color;
                        optionPane.setVisible(false);
                        repaint();
                        return;
                    }
                }

                 {
                    super();
                }
            });
            jpanel.add(jbutton);
            jpanel.add(new JLabel("Background Color"));
            jbutton = new JButton("change");
            jbutton.setBackground(GenePane.BGCOLOR);
            jbutton.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent actionevent) {
                    Color color = JColorChooser.showDialog(null, "Choose Background Color", GenePane.BGCOLOR);
                    if(color == null) {
                        return;
                    } else {
                        GenePane.BGCOLOR = color;
                        optionPane.setVisible(false);
                        regionPanel.setBackground(GenePane.BGCOLOR);
                        alternativePatternPanel.setBackground(GenePane.BGCOLOR);
                        alternativePatternPanel.resetAllAlternativeContainersColor();
                        labelPanel.setBackground(GenePane.BGCOLOR.darker());
                        repaint();
                        return;
                    }
                }

                 {
                    super();
                }
            });
            jpanel.add(jbutton);
            jpanel.add(new JLabel("Range Color"));
            jbutton = new JButton("change");
            jbutton.setBackground(GenePane.ALTERNATIVERANGE_COLOR);
            jbutton.addActionListener(new ActionListener() {

                public void actionPerformed(ActionEvent actionevent) {
                    Color color = JColorChooser.showDialog(null, "Choose Range Color", GenePane.ALTERNATIVERANGE_COLOR);
                    if(color == null) {
                        return;
                    } else {
                        GenePane.ALTERNATIVERANGE_COLOR = color;
                        optionPane.setVisible(false);
                        repaint();
                        return;
                    }
                }

                 {
                    super();
                }
            });
            jpanel.add(jbutton);
            optionPane = new JFrame();
            optionPane.getContentPane().add(jpanel);
        }

        void showOptionPane() {
            optionPane.pack();
            optionPane.setVisible(true);
        }


        public Option() {
            super("Option", GenePane.optionIcon);
            optionPane = null;
            mnemonic = 'o';
            toolTipText = "Show Options";
            init();
        }
    }

    class ScrollRight extends Action {

        public void actionPerformed(ActionEvent actionevent) {
            hScroll(-50);
        }

        public ScrollRight() {
            super("Right", GenePane.rightIcon);
            mnemonic = 'f';
            toolTipText = "Scroll Right";
        }
    }

    class ScrollLeft extends Action {

        public void actionPerformed(ActionEvent actionevent) {
            hScroll(50);
        }

        public ScrollLeft() {
            super("Left", GenePane.leftIcon);
            mnemonic = 'b';
            toolTipText = "Scroll Left";
        }
    }

    class Overview extends Action {

        public void actionPerformed(ActionEvent actionevent) {
            ovFrame.setVisible(true);
        }

        public Overview() {
            super("Overview", GenePane.viewIcon);
            mnemonic = 'o';
            toolTipText = "Show Overview";
        }
    }

    abstract class Action extends AbstractAction {

        protected char mnemonic;
        protected String toolTipText;

        char getMnemonic() {
            return mnemonic;
        }

        String getToolTipText() {
            return toolTipText;
        }

        public Action(String s, ImageIcon imageicon) {
            super(s, imageicon);
        }
    }


    private static final int CENTER_PANE_PREFERRED_WIDTH = 600;
    private static final int CENTER_PANE_MAX_WIDTH = 1200;
    private static final int DEFAULT_MAGPOWER = 100;
    private static final int MIN_MAGPOWER = 50;
    private static final int MAX_MAGPOWER = 5000;
    private static final int MAGPOWER_SPINNER_INCREMENT = 50;
    private static final int INITIAL_LEFT_MARGIN = 5;
    private static final int MAX_VOFFSET = 20;
    protected static Color BGCOLOR = new Color(208, 208, 208);
    protected static Color ORFCOLOR1 = new Color(0, 240, 0);
    protected static Color ORFCOLOR2;
    protected static Color UTRCOLOR1 = new Color(0, 0, 240);
    protected static Color UTRCOLOR2;
    protected static Color INTRONCOLOR1;
    protected static Color INTRONCOLOR2;
    protected static Color HIGHLIGHTCOLOR1 = new Color(240, 0, 0);
    protected static Color HIGHLIGHTCOLOR2 = new Color(160, 0, 0);
    protected static Color ALTERNATIVERANGE_COLOR = new Color(255, 225, 255);
    private static ImageIcon zoominIcon = null;
    private static ImageIcon zoomoutIcon = null;
    private static ImageIcon leftIcon = null;
    private static ImageIcon rightIcon = null;
    private static ImageIcon searchIcon = null;
    private static ImageIcon viewIcon = null;
    private static ImageIcon optionIcon = null;
    private JTextField messagePanel;
    private GenomeViewPane regionPanel;
    private labelPane labelPanel;
    private rulerPane rulerPanelA;
    private rulerPane rulerPanelB;
    private overviewPane ovPanel;
    private JFrame ovFrame;
    private JScrollPane spatternScrollPane;
    private SplicingPatternPane alternativePatternPanel;
    private ZoomSpinner zoomSpinner;
    private AnnotationFrame annotationFrame;
    private int hOffset;
    private int vOffset;
    private int panelHeight;
    private int panelWidth;
    private int basePerPixel;
    private int maxZoomScale;
    private int magPower;
    private int seqFrom;
    private int seqTo;
    private Point mouseFrom;
    private HashMap gffRecordsMap;
    private ArrayList drawnGeneSegmentList;
    private ArrayList alternativeExonPixelRectList;
    private GffRecord focusedGff;
    private GffRecord pointedGff;
    private Gene focusedGene;
    private int focusedRow;
    private String species;
    private URL baseURL;
    private AppletContext appletContext;
    private Locus locus;
    private ArrayList alternaLineList;
    private ArrayList lineInAlternaRange;
    private TreeSet lineNoAlternaRange;
    private JTextArea text;
    private JScrollPane scroll;
    private boolean equalScaleView;
    private EqualViewCalculator equalViewCalculator;
    private int oneExonLenEqualView;

    public GenePane() {
        messagePanel = new JTextField();
        regionPanel = new GenomeViewPane();
        labelPanel = new labelPane();
        rulerPanelA = new rulerPane();
        rulerPanelB = new rulerPane();
        ovPanel = new overviewPane();
        ovFrame = new JFrame();
        spatternScrollPane = null;
        alternativePatternPanel = new SplicingPatternPane();
        zoomSpinner = new ZoomSpinner(new SpinnerNumberModel(100, 50, 5000, 50));
        annotationFrame = null;
        hOffset = 0;
        vOffset = 20;
        panelHeight = 0;
        panelWidth = 600;
        basePerPixel = 0;
        maxZoomScale = 0;
        magPower = 100;
        seqFrom = -1;
        seqTo = -1;
        mouseFrom = null;
        gffRecordsMap = new HashMap();
        drawnGeneSegmentList = new ArrayList();
        alternativeExonPixelRectList = new ArrayList();
        focusedGff = null;
        pointedGff = null;
        focusedGene = null;
        focusedRow = -1;
        species = null;
        baseURL = null;
        appletContext = null;
        locus = null;
        alternaLineList = null;
        lineInAlternaRange = null;
        lineNoAlternaRange = null;
        text = null;
        scroll = null;
        equalScaleView = false;
        equalViewCalculator = null;
        oneExonLenEqualView = -1;
    }

    public void init(String s, URL url, URL url1) {
        try {
            species = s;
            readLocusData(url1.toString());
            readGff(url);
            ArrayList arraylist = locus.getSplicingPatternList();
            ArrayList arraylist1 = locus.getGeneList();
            equalViewCalculator = new EqualViewCalculator(arraylist, arraylist1, getStrand());
            equalViewCalculator.calculate();
            setOneExonLenEqualView(equalViewCalculator.getLineList().size() - 1);
            messagePanel.setEditable(false);
            layoutComponents();
            ovPanel.init();
            ovFrame.setTitle("Overview");
            ovFrame.getContentPane().add(ovPanel);
            ovFrame.pack();
            ovFrame.setVisible(true);
            ovFrame.setResizable(false);
            createAnnotationFrame();
        }
        catch(Exception exception) {
            exception.printStackTrace();
        }
    }

    private void setHOffset(int i) {
        hOffset = i;
    }

    private void setOneExonLenEqualView(int i) {
        oneExonLenEqualView = locus.getGenomeLen() / i;
        maxZoomScale = 4 * basePerPixel;
    }

    private int getOneExonLenEqualView() {
        return oneExonLenEqualView;
    }

    private int translateGenomePos(int i) {
        int j = -1;
        if(getStrand() == 1)
            j = locus.getMinGenomePos() + i;
        else
        if(getStrand() == -1)
            j = locus.getMaxGenomePos() - i;
        return j;
    }

    private int getLeftLinePos(int i, int j) {
        Exon aexon[] = locus.getGene(i).getExons();
        return getLeftLinePos(i, aexon[j]);
        Exception exception;
        exception;
        return 0;
    }

    private int getRightLinePos(int i, int j) {
        Exon aexon[] = locus.getGene(i).getExons();
        return getRightLinePos(i, aexon[j]);
        Exception exception;
        exception;
        return 0;
    }

    private int getLeftLinePos(int i, Exon exon) {
        int j = exon.getEqualViewLeft();
        if(j == -1)
            j = exon.getGenomeLeft();
        return genomePos2equalViewColumnNo(j);
    }

    private int getRightLinePos(int i, Exon exon) {
        int j = exon.getEqualViewRight();
        if(j == -1)
            j = exon.getGenomeRight();
        return genomePos2equalViewColumnNo(j);
    }

    private int genomePos2equalViewColumnNo(int i, int j) {
        TreeMap treemap = equalViewCalculator.getLineList();
        int k = -1;
        Integer integer = (Integer)treemap.get(new Integer(i));
        if(integer != null)
            k = integer.intValue();
        else
        if(j != 0) {
            int l = 1;
            do {
                int i1 = i - l;
                if(i1 <= locus.getMinGenomePos()) {
                    k = 0;
                    break;
                }
                int j1 = i + l;
                if(j1 >= locus.getMaxGenomePos()) {
                    k = treemap.size();
                    break;
                }
                Integer integer1 = (Integer)treemap.get(new Integer(i1));
                if(integer1 != null) {
                    k = integer1.intValue();
                    break;
                }
                Integer integer2 = (Integer)treemap.get(new Integer(j1));
                if(integer2 == null)
                    continue;
                k = integer2.intValue();
                break;
            } while(++l != j);
        }
        if(getStrand() == -1)
            k = treemap.size() - 1 - k;
        return k;
    }

    private int genomePos2equalViewColumnNo(int i) {
        return genomePos2equalViewColumnNo(i, -1);
    }

    public void createImageIcons(ImageIconMaker imageiconmaker) {
        zoominIcon = imageiconmaker.getZoominImageIcon();
        zoomoutIcon = imageiconmaker.getZoomoutImageIcon();
        leftIcon = imageiconmaker.getLeftImageIcon();
        rightIcon = imageiconmaker.getRightImageIcon();
        searchIcon = imageiconmaker.getSearchImageIcon();
        viewIcon = imageiconmaker.getViewImageIcon();
        optionIcon = imageiconmaker.getOptionImageIcon();
    }

    private JComponent constructToolBar() {
        JToolBar jtoolbar = new JToolBar();
        jtoolbar.setLayout(new FlowLayout(0));
        jtoolbar.setRollover(true);
        Action aaction[] = {
            new Overview(), new ScrollLeft(), new ScrollRight(), new Option()
        };
        for(int i = 0; i < aaction.length; i++) {
            JButton jbutton1 = jtoolbar.add(aaction[i]);
            jbutton1.setMargin(new Insets(0, 0, 0, 0));
            jbutton1.setMnemonic(aaction[i].getMnemonic());
            jbutton1.setToolTipText(aaction[i].getToolTipText());
        }

        jtoolbar.addSeparator();
        jtoolbar.add(new JLabel("Zoom : "));
        jtoolbar.add(zoomSpinner);
        jtoolbar.add(new JLabel(" %"));
        JButton jbutton = new JButton("Show monospaced view");
        jbutton.setToolTipText("Toggle for proportional or monospaced view");
        jbutton.setMnemonic('v');
        jbutton.addActionListener(new ChangeView());
        jtoolbar.add(jbutton);
        if(locus.numHomoloGenes() > 1) {
            JButton jbutton2 = new JButton("HomoloGene");
            jbutton2.setToolTipText("HomoloGene");
            jbutton2.setMnemonic('h');
            jbutton2.addActionListener(new HomoloGene());
            jtoolbar.add(jbutton2);
        }
        return jtoolbar;
    }

    private JComponent constructCenterPane() {
        JPanel jpanel = new JPanel();
        jpanel.setLayout(new BorderLayout());
        jpanel.add(rulerPanelA, "North");
        Object obj = new JPanel();
        ((JPanel) (obj)).setLayout(new BoxLayout(((Container) (obj)), 0));
        ((JPanel) (obj)).add(labelPanel);
        ((JPanel) (obj)).add(regionPanel);
        jpanel.add(((Component) (obj)), "Center");
        jpanel.add(rulerPanelB, "South");
        alternativePatternPanel.layoutComponents();
        spatternScrollPane = new JScrollPane(alternativePatternPanel, 20, 30);
        obj = new JSplitPane(0);
        ((JSplitPane) (obj)).setDividerLocation(300);
        ((JSplitPane) (obj)).add(jpanel, "top");
        ((JSplitPane) (obj)).add(spatternScrollPane, "bottom");
        ((JSplitPane) (obj)).setPreferredSize(new Dimension(600, 400));
        ((JSplitPane) (obj)).setMaximumSize(new Dimension(1200, 600));
        ((JSplitPane) (obj)).setAlignmentX(0.0F);
        ((JSplitPane) (obj)).setOneTouchExpandable(true);
        return ((JComponent) (obj));
    }

    private void createAnnotationFrame() {
        annotationFrame = new AnnotationFrame(this);
        annotationFrame.setSize(780, 500);
    }

    private void layoutComponents() {
        setLayout(new BorderLayout());
        setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
        add(constructToolBar(), "North");
        add(messagePanel, "South");
        add(constructCenterPane(), "Center");
    }

    private void readGff(URL url) {
        try {
            BufferedReader bufferedreader;
label0:
            {
                bufferedreader = new BufferedReader(new InputStreamReader(url.openStream()));
                Object obj = null;
                do {
                    String s;
                    do {
                        if((s = bufferedreader.readLine()) == null)
                            break label0;
                        s = s.trim();
                    } while(s.startsWith("#"));
                    GffRecord gffrecord = GffRecord.parseGFF(s);
                    if(gffrecord == null)
                        break;
                    if(getStrand() == -1) {
                        int i = gffrecord.getStart();
                        int j = gffrecord.getEnd();
                        if(i > j) {
                            gffrecord.setStart(complementPos2directPos(i));
                            gffrecord.setEnd(complementPos2directPos(j));
                        } else {
                            gffrecord.setStart(complementPos2directPos(j));
                            gffrecord.setEnd(complementPos2directPos(i));
                        }
                    }
                    addGffData(gffrecord);
                } while(true);
                System.err.println("NO GFF DATA: abort.");
            }
            bufferedreader.close();
            break MISSING_BLOCK_LABEL_169;
        }
        catch(Exception exception) {
            exception.printStackTrace();
        }
    }

    private void addGffData(GffRecord gffrecord) {
        ArrayList arraylist = (ArrayList)gffRecordsMap.get(gffrecord.getLabel());
        if(arraylist == null) {
            arraylist = new ArrayList();
            gffRecordsMap.put(gffrecord.getLabel(), arraylist);
        }
        arraylist.add(gffrecord);
    }

    public int getStrand() {
        int i = 0;
        if(locus != null)
            i = locus.getStrand();
        return i;
    }

    public int complementPos2directPos(int i) {
        int j = i;
        if(locus != null)
            j = locus.complementPos2directPos(i);
        return j;
    }

    public int directPos2complementPos(int i) {
        return complementPos2directPos(i);
    }

    private void readLocusData(String s) throws ParserConfigurationException, SAXException, IOException {
        SAXParserFactory saxparserfactory = SAXParserFactory.newInstance();
        SAXParser saxparser = saxparserfactory.newSAXParser();
        LocusSaxHandler locussaxhandler = new LocusSaxHandler();
        saxparser.parse(s, locussaxhandler);
        locus = locussaxhandler.getLocus();
    }

    private int searchPosition(String s, String s1, int i) {
        if((i = s.indexOf(s1, ++i)) >= 0)
            return i;
        if(seqFrom == -1 && (i = s.indexOf(s1)) >= 0)
            return i;
        else
            return -1;
    }

    private void hScroll(int i) {
        int j = getBasePerPixel();
        if(j * (panelWidth - 40) > locus.getMaxGenomePos() - locus.getMinGenomePos()) {
            hOffset = (locus.getMaxGenomePos() + locus.getMinGenomePos()) / 2 - (j * panelWidth) / 2;
            return;
        }
        hOffset -= i * j;
        if(hOffset < locus.getMinGenomePos() - 20 * j) {
            hOffset = locus.getMinGenomePos() - 20 * j;
            messagePanel.setText("Left scroll limit!");
        } else
        if(hOffset > locus.getMaxGenomePos() - (panelWidth - 20) * j) {
            hOffset = locus.getMaxGenomePos() - (panelWidth - 20) * j;
            messagePanel.setText("Right scroll limit!");
        }
        repaintPanel();
    }

    private void vScroll(int i) {
        vOffset += i;
        if(vOffset > 20) {
            vOffset = 20;
            messagePanel.setText("Top scroll limit!");
        } else {
            int j = panelHeight - regionPanel.getHeight();
            if(vOffset + j < 0) {
                vOffset = -j;
                messagePanel.setText("Bottom scroll limit!");
            }
        }
        repaintPanel();
    }

    private void initZoomScale() {
        int i = panelWidth - 10;
        basePerPixel = (int)Math.rint((double)locus.getGenomeLen() / (double)i);
        if(basePerPixel > 1)
            for(; locus.getGenomeLen() / (basePerPixel - 1) > i; basePerPixel++);
        maxZoomScale = 2 * basePerPixel;
    }

    private int getBasePerPixel() {
        int i = (int)Math.ceil((double)basePerPixel / ((double)magPower / 100D));
        return i > 1 ? i - 1 : 1;
    }

    private void zoom(int i, int j) {
        int k = getBasePerPixel();
        if(k == 1 && j > 0)
            return;
        int l = magPower + j;
        if(l > 5000)
            l = 5000;
        else
        if(l < 50)
            l = 50;
        zoomSpinner.setValue(new Integer(l));
        if(l <= 100) {
            hOffset = locus.getMinGenomePos();
        } else {
            int i1 = getBasePerPixel();
            if(i == -1)
                i = focusedGff == null ? (locus.getMaxGenomePos() + locus.getMinGenomePos()) / 2 : (focusedGff.getEnd() + focusedGff.getStart()) / 2;
            int k1 = (i - hOffset) / k;
            hOffset = i - i1 * k1;
        }
        int j1 = (l - magPower) / 5;
        Dimension dimension = getSize();
        for(int l1 = 0; l1 < 5; l1++) {
            magPower += j1;
            repaintPanel();
            try {
                Thread.sleep(10L);
            }
            catch(Exception exception) { }
            paintImmediately(0, 0, dimension.width, dimension.height);
        }

    }

    private void repaintPanel() {
        regionPanel.repaint();
        regionPanel.revalidate();
    }

    public void doKeyPressed(KeyEvent keyevent) {
        if(!keyevent.isControlDown() && !keyevent.isAltDown())
            return;
        switch(keyevent.getKeyCode()) {
        default:
            break;

        case 66: // 'B'
            if(keyevent.isControlDown()) {
                hScroll(-10);
                break;
            }
            if(keyevent.isAltDown())
                hScroll(-10);
            break;

        case 70: // 'F'
            if(keyevent.isControlDown()) {
                hScroll(10);
                break;
            }
            if(keyevent.isAltDown())
                hScroll(10);
            break;

        case 78: // 'N'
            if(keyevent.isControlDown())
                zoom(-1, -10);
            break;

        case 80: // 'P'
            if(keyevent.isControlDown())
                zoom(-1, 10);
            break;

        case 81: // 'Q'
            if(!keyevent.isControlDown());
            break;
        }
    }

    public void setBaseURL(URL url) {
        baseURL = url;
    }

    public URL getBaseURL() {
        return baseURL;
    }

    public void setAppletContext(AppletContext appletcontext) {
        appletContext = appletcontext;
    }

    public AppletContext getAppletContext() {
        return appletContext;
    }

    public String getSpecies() {
        return species;
    }

    public Locus getLocus() {
        return locus;
    }

    public Gene getGeneByCdnaId(String s) {
        for(int i = 0; i < locus.getGeneList().size(); i++) {
            Gene gene = (Gene)locus.getGeneList().get(i);
            if(gene.getCdnaId().equals(s))
                return gene;
        }

        return null;
    }

    static  {
        ORFCOLOR2 = Color.darkGray;
        UTRCOLOR2 = Color.darkGray;
        INTRONCOLOR1 = Color.lightGray;
        INTRONCOLOR2 = Color.BLACK;
    }




















































}