// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:04 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   AnnotationFrame.java

import javax.swing.border.CompoundBorder;
import javax.swing.border.TitledBorder;

class AnnotationFrame$BottomSpacingTitledBorder extends CompoundBorder {

    AnnotationFrame$BottomSpacingTitledBorder(String s) {
        super(new AnnotationFrame$BottomSpacingEmptyBorder(AnnotationFrame.this), new TitledBorder(s));
    }
}