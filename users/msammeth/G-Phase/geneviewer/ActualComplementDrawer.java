package geneviewer;

// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:03 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   ActualComplementDrawer.java

import java.awt.Rectangle;

public class ActualComplementDrawer extends GeneSegmentDrawerImpl {

    public ActualComplementDrawer() {
    }

    protected Rectangle calcDrawingRect(Range range, int i) {
        int j = genomePosition2Pixel(complementPos2directPos(range.getLeft()));
        int k = genomePosition2Pixel(complementPos2directPos(range.getRight()));
        return new Rectangle(j, _vOffset, (k - j) + 1, i);
    }
}