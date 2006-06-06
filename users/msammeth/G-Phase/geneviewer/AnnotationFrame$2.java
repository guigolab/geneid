// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:03 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   AnnotationFrame.java

import java.applet.AppletContext;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.PrintStream;
import java.net.URL;

class AnnotationFrame$2
    implements ActionListener {

    public void actionPerformed(ActionEvent actionevent) {
        try {
            AppletContext appletcontext = AnnotationFrame.access$000(AnnotationFrame.this).getAppletContext();
            appletcontext.showDocument(new URL(AnnotationFrame.access$200(AnnotationFrame.this)), "_blank");
        }
        catch(Exception exception) {
            System.err.println(exception);
        }
    }

    AnnotationFrame$2() {
        super();
    }
}