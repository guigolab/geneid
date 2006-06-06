// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:07 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   IntronHighlighter.java

import javax.swing.JScrollPane;

public class IntronHighlighter extends AbstractSequenceHighlighter {

    private static final int INCREMENT_CARET_LINES = 18;

    public IntronHighlighter(Locus locus, GeneShape geneshape, SequenceTextPane sequencetextpane, JScrollPane jscrollpane) {
        _locus = locus;
        _highlightGeneShape = geneshape;
        _sequenceTextPane = sequencetextpane;
        _scrollPane = jscrollpane;
    }

    public int highlight() {
        int i = 0;
        Range arange[] = getHighlightRanges();
        String s = String.valueOf(arange[0].getLeft()) + "," + String.valueOf(arange[0].getRight()) + " intron";
        String s1 = _sequenceTextPane.getText();
        int j = s1.indexOf(s);
        if(j >= 0) {
            j += s.length() + 1;
            int k = arange[0].length() + (arange[0].length() - 1) / 50;
            _sequenceTextPane.highlight(j, k);
            i = calcCaretPosition(j, getIncrementCaretPosition(k));
        }
        _sequenceTextPane.setPreviousCaretPosition(i);
        _sequenceTextPane.setCaretPosition(i);
        return j;
    }

    public Range[] getHighlightRanges() {
        Range arange[] = new Range[1];
        arange[0] = new Range();
        switch(_locus.getStrand()) {
        case -1: 
            arange[0].setLeft(_locus.directPos2complementPos(_highlightGeneShape.getGffRecord().getStart()));
            arange[0].setRight(_locus.directPos2complementPos(_highlightGeneShape.getGffRecord().getEnd()));
            break;

        case 1: // '\001'
        default:
            arange[0].setLeft(_highlightGeneShape.getGffRecord().getStart());
            arange[0].setRight(_highlightGeneShape.getGffRecord().getEnd());
            break;
        }
        return arange;
    }

    private int getIncrementCaretPosition(int i) {
        char c = '\u0396';
        return c >= i ? i + 2 : c;
    }
}