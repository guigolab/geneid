// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:06 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   GeneSegmentDrawerImpl.java

import java.awt.*;

public abstract class GeneSegmentDrawerImpl
    implements GeneSegmentDrawer {

    protected Range _locusRange;
    protected Graphics2D _graphics2D;
    protected int _hOffset;
    protected int _vOffset;
    protected float _basePerPixel;
    protected int _drawingAreaWidth;
    protected int _drawaingAreaHeight;
    protected float _magPower;

    public GeneSegmentDrawerImpl() {
    }

    public void draw(GeneSegment genesegment) {
        ColoredRange acoloredrange[] = genesegment.getColoredRanges();
        for(int i = 0; i < acoloredrange.length; i++) {
            Rectangle rectangle = calcDrawingRect(acoloredrange[i].getRange(), genesegment.getHeight());
            GradientPaint gradientpaint = new GradientPaint(rectangle.x, rectangle.y, acoloredrange[i].getColor1(), rectangle.x, (rectangle.y + rectangle.height) - 1, acoloredrange[i].getColor2());
            _graphics2D.setPaint(gradientpaint);
            _graphics2D.fill(rectangle);
        }

    }

    protected abstract Rectangle calcDrawingRect(Range range, int i);

    protected int genomePosition2Pixel(int i) {
        int j = _hOffset;
        int k = _hOffset + (int)((float)_drawingAreaWidth * _basePerPixel);
        if(i < j)
            i = j;
        else
        if(i > k)
            i = k;
        return (int)((float)(i - _hOffset) / _basePerPixel);
    }

    protected int complementPos2directPos(int i) {
        return SequenceUtils.complementPos2directPos(_locusRange.getRight(), _locusRange.getLeft(), i);
    }

    public void setGraphics2D(Graphics2D graphics2d) {
        _graphics2D = graphics2d;
    }

    public void setVOffset(int i) {
        _vOffset = i;
    }
}