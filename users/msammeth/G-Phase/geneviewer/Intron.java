// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:07 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   Intron.java

import java.awt.Color;

public class Intron extends GeneSegment {

    private static final int INTRON_SHAPE_HEIGHT = 3;
    private static final Color INTRONCOLOR = new Color(0, 0, 0);

    public Intron() {
        _segmentType = 2;
    }

    public ColoredRange[] getColoredRanges() {
        ColoredRange acoloredrange[] = new ColoredRange[1];
        acoloredrange[0] = new ColoredRange();
        acoloredrange[0].setColor1(INTRONCOLOR);
        acoloredrange[0].setColor2(INTRONCOLOR);
        acoloredrange[0].setRange(new Range(_genomeLeft, _genomeRight));
        return acoloredrange;
    }

    public int getHeight() {
        return 3;
    }

}