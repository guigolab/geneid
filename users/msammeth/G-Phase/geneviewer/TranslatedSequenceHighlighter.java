// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:07 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   TranslatedSequenceHighlighter.java

import java.util.StringTokenizer;
import javax.swing.JScrollPane;

public class TranslatedSequenceHighlighter extends ExonHighlighter {

    public TranslatedSequenceHighlighter(Locus locus, GeneShape geneshape, SequenceTextPane sequencetextpane, JScrollPane jscrollpane, Exon exon, boolean flag) {
        super(locus, geneshape, sequencetextpane, jscrollpane, exon, flag);
    }

    public int highlight() {
        int i = 0;
        boolean flag = false;
        Range arange[] = getHighlightCdnaRanges();
        if(arange.length > 0) {
            int j = highlight(arange[0].getLeft(), arange[0].getRight());
            for(int k = 1; k < arange.length; k++)
                highlight(arange[k].getLeft(), arange[k].getRight());

            i = calcCaretPosition(j, 0);
        }
        _sequenceTextPane.setPreviousCaretPosition(i);
        _sequenceTextPane.setCaretPosition(i);
        return i;
    }

    public Range[] getHighlightCdnaRanges() {
        Exon exon = _highlightGeneShape.getExon();
        int i = exon.getCdnaRange().getLeft();
        int j = exon.getCdnaRange().getRight();
        Range arange[] = getHighlightRanges();
        Range arange1[] = new Range[arange.length];
        for(int k = 0; k < arange.length; k++) {
            arange1[k] = new Range();
            arange1[k].setLeft(i + Math.abs(exon.getGenomeLeft() - arange[k].getLeft()));
            arange1[k].setRight(j - Math.abs(exon.getGenomeRight() - arange[k].getRight()));
        }

        return arange1;
    }

    private int highlight(int i, int j) {
        String s = _sequenceTextPane.getText();
        int k = 0;
        int l = 0;
        int i1 = 0;
        int j1 = -1;
        boolean flag = false;
        for(StringTokenizer stringtokenizer = new StringTokenizer(s, "\n", true); stringtokenizer.hasMoreTokens();) {
            String s1 = stringtokenizer.nextToken();
            if(k == 3) {
                if(s1.charAt(0) == ' ') {
                    i++;
                    j++;
                }
                if(s1.charAt(1) == ' ') {
                    i++;
                    j++;
                }
                if(i <= 60) {
                    int k1 = i + l;
                    int i2 = j > 60 ? 60 : j;
                    i2 += l;
                    _sequenceTextPane.highlight(k1 - 1, (i2 - k1) + 1);
                    j1 = k1 - 1;
                }
                i1 += s1.length();
            } else
            if(k > 3 && isNucleotideLine(s1)) {
                if(i <= i1 + 60 && i1 + 1 <= j) {
                    int l1 = i > i1 ? l + (i - i1) : l;
                    int j2 = j > i1 + 60 ? l + 60 : l + (j - i1);
                    _sequenceTextPane.highlight(l1 - 1, (j2 - l1) + 1);
                    if(j1 == -1)
                        j1 = l1 - 1;
                }
                i1 += s1.length();
            }
            l += s1.length();
            k++;
        }

        return j1;
    }

    private boolean isNucleotideLine(String s) {
        boolean flag = false;
        char c = s.charAt(0);
        if(c == 'A' || c == 'T' || c == 'G' || c == 'C' || c == 'a' || c == 't' || c == 'g' || c == 'c')
            flag = true;
        return flag;
    }
}