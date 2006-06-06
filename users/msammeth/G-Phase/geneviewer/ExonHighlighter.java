// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:04 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   ExonHighlighter.java

import javax.swing.JScrollPane;

public class ExonHighlighter extends AbstractSequenceHighlighter {

    private static final int INCREMENT_CARET_LINES = 8;
    private Exon _alternativeExon;
    private boolean _bHighlightAlternativeExon;

    public ExonHighlighter(Locus locus, GeneShape geneshape, SequenceTextPane sequencetextpane, JScrollPane jscrollpane, Exon exon, boolean flag) {
        _alternativeExon = null;
        _bHighlightAlternativeExon = true;
        _locus = locus;
        _highlightGeneShape = geneshape;
        _sequenceTextPane = sequencetextpane;
        _scrollPane = jscrollpane;
        _alternativeExon = exon;
        _bHighlightAlternativeExon = flag;
    }

    public int highlight() {
        Exon exon = _highlightGeneShape.getExon();
        String s = String.valueOf(exon.getGenomeLeft()) + ".." + String.valueOf(exon.getGenomeRight()) + " exon";
        int i = 0;
        String s1 = _sequenceTextPane.getText();
        int j = s1.indexOf(s);
        if(j >= 0) {
            j += s.length() + 1;
            Range arange[] = getHighlightRanges();
            for(int k = 0; k < arange.length; k++)
                highlight(arange[k], j);

            i = calcCaretPosition(j, getIncrementCaretPosition());
        }
        _sequenceTextPane.setPreviousCaretPosition(i);
        _sequenceTextPane.setCaretPosition(i);
        return j;
    }

    public int getHighlightLength(Range range) {
        Exon exon = _highlightGeneShape.getExon();
        int i = range.length();
        int j = Math.abs(range.getLeft() - exon.getGenomeLeft()) % 50;
        int k = 0;
        if(range.getLeft() < range.getRight()) {
            for(int l = range.getLeft() + (50 - j); l <= range.getRight();) {
                l += 50;
                k++;
            }

        } else {
            for(int i1 = range.getLeft() - (50 - j); i1 >= range.getRight();) {
                i1 -= 50;
                k++;
            }

        }
        return i + k;
    }

    public Range[] getHighlightRanges() {
        Exon exon = _highlightGeneShape.getExon();
        Range arange[] = null;
        int i = exon.getGenomeLeft();
        int j = exon.getGenomeRight();
        if(_alternativeExon != null) {
            int k = _alternativeExon.getGenomeLeft();
            int l = _alternativeExon.getGenomeRight();
            if(_bHighlightAlternativeExon) {
                arange = new Range[1];
                arange[0] = new Range();
                arange[0].setLeft(k);
                arange[0].setRight(l);
            } else
            if(i != k && j != l) {
                arange = new Range[2];
                arange[0] = new Range();
                arange[1] = new Range();
                arange[0].setLeft(i);
                arange[1].setRight(j);
                if(i < j) {
                    arange[0].setRight(k - 1);
                    arange[1].setLeft(l + 1);
                } else {
                    arange[0].setRight(k + 1);
                    arange[1].setLeft(l - 1);
                }
            } else {
                arange = new Range[1];
                arange[0] = new Range();
                if(i != k) {
                    arange[0].setLeft(i);
                    if(i < j)
                        arange[0].setRight(k - 1);
                    else
                        arange[0].setRight(k + 1);
                } else
                if(j != l) {
                    if(i < j)
                        arange[0].setLeft(l + 1);
                    else
                        arange[0].setLeft(l - 1);
                    arange[0].setRight(j);
                }
            }
        } else {
            arange = new Range[1];
            arange[0] = new Range();
            arange[0].setLeft(i);
            arange[0].setRight(j);
        }
        return arange;
    }

    private int highlight(Range range, int i) {
        Exon exon = _highlightGeneShape.getExon();
        int j = Math.abs(range.getLeft() - exon.getGenomeLeft());
        int k = j + j / 50;
        int l = i + k;
        int i1 = getHighlightLength(range);
        _sequenceTextPane.highlight(l, i1);
        return l;
    }

    private int getIncrementCaretPosition() {
        int i = _highlightGeneShape.getGffRecord().getLength();
        char c = '\u0198';
        return c >= i ? i + 50 : c;
    }
}