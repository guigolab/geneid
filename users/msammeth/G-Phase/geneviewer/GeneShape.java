// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:06 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   GeneShape.java

import java.awt.Rectangle;

public class GeneShape {

    public static final int GENESEGMENT_UNKNOWN = 0;
    public static final int GENESEGMENT_EXON = 1;
    public static final int GENESEGMENT_INTRON = 2;
    protected Rectangle _segmentRect;
    protected int _segmentType;
    GffRecord _gffRecord;
    Exon _exon;
    ColoredRange _coloredRectangle;
    Rectangle highlightedArea;

    public GeneShape() {
        _segmentRect = null;
        _segmentType = 0;
        _gffRecord = null;
        _exon = null;
        _coloredRectangle = null;
        highlightedArea = null;
    }

    public boolean isHighlightedArea(int i, int j) {
        boolean flag = false;
        if(highlightedArea != null && highlightedArea.contains(i, j))
            flag = true;
        return flag;
    }

    public boolean isHighlighted() {
        return highlightedArea != null;
    }

    public String toString() {
        StringBuffer stringbuffer = new StringBuffer();
        switch(_segmentType) {
        case 0: // '\0'
            stringbuffer.append("Unknown");
            break;

        case 1: // '\001'
            stringbuffer.append("EXON");
            break;

        case 2: // '\002'
            stringbuffer.append("INTRON");
            break;
        }
        stringbuffer.append(':');
        stringbuffer.append(_gffRecord.toString());
        return stringbuffer.toString();
    }

    public Rectangle getSegmentRect() {
        return _segmentRect;
    }

    public void setSegmentRect(Rectangle rectangle) {
        _segmentRect = rectangle;
    }

    public int getSegmentType() {
        return _segmentType;
    }

    public void setSegmentType(int i) {
        _segmentType = i;
    }

    public GffRecord getGffRecord() {
        return _gffRecord;
    }

    public void setGffRecord(GffRecord gffrecord) {
        _gffRecord = gffrecord;
    }

    public Rectangle getHighlightedArea() {
        return highlightedArea;
    }

    public void setHighlightedArea(Rectangle rectangle) {
        highlightedArea = rectangle;
    }

    public Exon getExon() {
        return _exon;
    }

    public void setExon(Exon exon) {
        _exon = exon;
    }
}