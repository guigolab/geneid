// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:04 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   ColoredRange.java

import java.awt.Color;

public class ColoredRange {

    private Color _color1;
    private Color _color2;
    private Range _range;

    public ColoredRange() {
    }

    public Color getColor1() {
        return _color1;
    }

    public void setColor1(Color color) {
        _color1 = color;
    }

    public Color getColor2() {
        return _color2;
    }

    public void setColor2(Color color) {
        _color2 = color;
    }

    public Range getRange() {
        return _range;
    }

    public void setRange(Range range) {
        _range = range;
    }
}