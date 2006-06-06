// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:04 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   Exon.java

import java.awt.Color;

public class Exon extends GeneSegment {

    public static final int EXONTYPE_UNKNOWN = 0;
    public static final int EXONTYPE_UTRONLY = 1;
    public static final int EXONTYPE_UTR_ORF = 2;
    public static final int EXONTYPE_ORFONLY = 3;
    public static final int EXONTYPE_ORF_UTR = 4;
    protected static Color ORFCOLOR1;
    protected static Color ORFCOLOR2;
    protected static Color UTRCOLOR1;
    protected static Color UTRCOLOR2;
    private static final int EXON_SHAPE_HEIGHT = 15;
    private Range _cdnaRange;
    protected Range _orfRange;
    protected Range _utrRange;
    private int _exonType;
    private int _equalViewLeft;
    private int _equalViewRight;

    public Exon() {
        _cdnaRange = null;
        _orfRange = null;
        _utrRange = null;
        _equalViewLeft = -1;
        _equalViewRight = -1;
        _exonType = 0;
        _segmentType = 1;
    }

    public Exon(int i, int j) {
        this();
        _genomeLeft = i;
        _genomeRight = j;
    }

    public ColoredRange[] getColoredRanges() {
        return null;
    }

    public int getHeight() {
        return 15;
    }

    public void setExonType(Range range) {
        int i = _cdnaRange.getLeft();
        int j = _cdnaRange.getRight();
        if(range.getLeft() <= i && j <= range.getRight())
            _exonType = 3;
        else
        if(j < range.getLeft())
            _exonType = 1;
        else
        if(range.getRight() < i)
            _exonType = 1;
        else
        if(i < range.getLeft())
            _exonType = 2;
        else
        if(range.getRight() < j)
            _exonType = 4;
    }

    public Range getCdnaRange() {
        return _cdnaRange;
    }

    public void setCdnaRange(Range range) {
        _cdnaRange = new Range(range.getLeft(), range.getRight());
    }

    public void setCdnaRange(int i, int j) {
        _cdnaRange = new Range(i, j);
    }

    public void setUtrRange(int i, int j) {
        _utrRange = new Range(i, j);
    }

    public void setOrfRange(int i, int j) {
        _orfRange = new Range(i, j);
    }

    public int getExonType() {
        return _exonType;
    }

    public void setExonType(int i) {
        _exonType = i;
    }

    public void setEqualViewLeft(int i) {
        _equalViewLeft = i;
    }

    public int getEqualViewLeft() {
        return _equalViewLeft;
    }

    public void setEqualViewRight(int i) {
        _equalViewRight = i;
    }

    public int getEqualViewRight() {
        return _equalViewRight;
    }

    static  {
        ORFCOLOR1 = GenePane.ORFCOLOR1;
        ORFCOLOR2 = GenePane.ORFCOLOR2;
        UTRCOLOR1 = GenePane.UTRCOLOR1;
        UTRCOLOR2 = GenePane.UTRCOLOR2;
    }
}