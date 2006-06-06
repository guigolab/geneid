// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:07 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   NullSequenceHighlighter.java


public class NullSequenceHighlighter extends AbstractSequenceHighlighter {

    private static NullSequenceHighlighter _uniqueInstance = null;

    private NullSequenceHighlighter() {
    }

    public static NullSequenceHighlighter getInstance() {
        if(_uniqueInstance == null)
            _uniqueInstance = new NullSequenceHighlighter();
        return _uniqueInstance;
    }

    public int highlight() {
        _sequenceTextPane.setPreviousCaretPosition(0);
        _sequenceTextPane.setCaretPosition(0);
        return 0;
    }

    public Range[] getHighlightRanges() {
        return new Range[0];
    }

}