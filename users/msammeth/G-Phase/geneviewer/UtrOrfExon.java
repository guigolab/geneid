// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:03 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   UtrOrfExon.java


public class UtrOrfExon extends Exon {

    public UtrOrfExon() {
    }

    public ColoredRange[] getColoredRanges() {
        ColoredRange acoloredrange[] = new ColoredRange[2];
        acoloredrange[0] = new ColoredRange();
        acoloredrange[0].setColor1(UTRCOLOR1);
        acoloredrange[0].setColor2(UTRCOLOR2);
        acoloredrange[0].setRange(_utrRange);
        acoloredrange[1] = new ColoredRange();
        acoloredrange[1].setColor1(ORFCOLOR1);
        acoloredrange[1].setColor2(ORFCOLOR2);
        acoloredrange[1].setRange(_orfRange);
        return acoloredrange;
    }
}