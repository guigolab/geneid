// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:07 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   SequenceTextPane.java

import java.awt.*;
import javax.swing.JTextPane;
import javax.swing.plaf.ComponentUI;
import javax.swing.text.*;

public class SequenceTextPane extends JTextPane {

    public static final String MONOSPACE_FONT_NAME = "monospaced";
    public static final int FONT_SIZE = 14;
    public static final String HIGHLIGHT_STYLE_NAME = "highlightSequence";
    static final Color HIGHLIGHT_COLOR = new Color(200, 0, 0);
    private int _previousCaretPosition;
    private int _currentScrollPosition;
    private Style _highlightStyle;

    public SequenceTextPane() {
        _previousCaretPosition = 0;
        _currentScrollPosition = 0;
        _highlightStyle = null;
        setFont(new Font("monospaced", 0, 14));
        setEditable(false);
        _highlightStyle = getStyledDocument().addStyle("highlightSequence", null);
        StyleConstants.setForeground(_highlightStyle, HIGHLIGHT_COLOR);
    }

    public void highlight(int i, int j) {
        getStyledDocument().setCharacterAttributes(i, j, _highlightStyle, true);
    }

    public boolean getScrollableTracksViewportWidth() {
        java.awt.Container container = getParent();
        javax.swing.plaf.TextUI textui = getUI();
        return textui.getPreferredSize(this).width <= container.getSize().width;
    }

    public void setPreviousCaretPosition(int i) {
        _previousCaretPosition = i;
    }

    public int getPreviousCaretPosition() {
        return _previousCaretPosition;
    }

}