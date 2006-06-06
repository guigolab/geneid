// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:06 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   GeneSegment.java


public abstract class GeneSegment
    implements Cloneable {

    public static final int SEGMENT_TYPE_UNKNOWN = 0;
    public static final int SEGMENT_TYPE_EXON = 1;
    public static final int SEGMENT_TYPE_INTRON = 2;
    protected int _genomeLeft;
    protected int _genomeRight;
    protected int _segmentType;

    public GeneSegment() {
        _genomeLeft = -1;
        _genomeRight = -1;
        _segmentType = 0;
    }

    public GeneSegment(int i, int j, int k) {
        _genomeLeft = i;
        _genomeRight = j;
        _segmentType = k;
    }

    public Object clone() throws CloneNotSupportedException {
        return super.clone();
    }

    public abstract ColoredRange[] getColoredRanges();

    public abstract int getHeight();

    public void acceptDrawer(GeneSegmentDrawer genesegmentdrawer) {
        genesegmentdrawer.draw(this);
    }

    public int getGenomeLeft() {
        return _genomeLeft;
    }

    public void setGenomeLeft(int i) {
        _genomeLeft = i;
    }

    public int getGenomeRight() {
        return _genomeRight;
    }

    public void setGenomeRight(int i) {
        _genomeRight = i;
    }

    public int getSegmentType() {
        return _segmentType;
    }

    public void setSegmentType(int i) {
        _segmentType = i;
    }
}