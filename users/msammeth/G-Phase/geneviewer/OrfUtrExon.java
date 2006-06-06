// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:07 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   OrfUtrExon.java


public class OrfUtrExon extends Exon {

    public OrfUtrExon() {
    }

    public ColoredRange[] getColoredRanges() {
        ColoredRange acoloredrange[] = new ColoredRange[2];
        acoloredrange[0] = new ColoredRange();
        acoloredrange[0].setColor1(ORFCOLOR1);
        acoloredrange[0].setColor2(ORFCOLOR2);
        acoloredrange[0].setRange(_orfRange);
        acoloredrange[1] = new ColoredRange();
        acoloredrange[1].setColor1(UTRCOLOR1);
        acoloredrange[1].setColor2(UTRCOLOR2);
        acoloredrange[1].setRange(_utrRange);
        return acoloredrange;
    }
}