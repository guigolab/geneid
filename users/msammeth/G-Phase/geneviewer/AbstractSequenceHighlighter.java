package geneviewer;

// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:03 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   AbstractSequenceHighlighter.java

import javax.swing.JScrollPane;

public abstract class AbstractSequenceHighlighter
    implements SequenceHighlighter {

    static final int TRANSLATED_SEQ_COLUMNS = 60;
    static final int SPLICING_SEQ_COLUMNS = 50;
    Locus _locus;
    GeneShape _highlightGeneShape;
    SequenceTextPane _sequenceTextPane;
    JScrollPane _scrollPane;

    public AbstractSequenceHighlighter() {
        _locus = null;
        _highlightGeneShape = null;
        _sequenceTextPane = null;
        _scrollPane = null;
    }

    int calcCaretPosition(int i, int j) {
        int k;
        if(i < _sequenceTextPane.getPreviousCaretPosition())
            k = i - 1;
        else
            k = i + j;
        if(k < 0)
            k = 0;
        else
        if(k > _sequenceTextPane.getText().length())
            k = _sequenceTextPane.getText().length();
        return k;
    }
}