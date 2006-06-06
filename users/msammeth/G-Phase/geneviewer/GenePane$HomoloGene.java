// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:05 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   GenePane.java

import java.applet.AppletContext;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.*;

class GenePane$HomoloGene
    implements ActionListener {

    public void actionPerformed(ActionEvent actionevent) {
        showHomoloGene();
    }

    private void showHomoloGene() {
        try {
            HashMap hashmap = GenePane.access$1300(GenePane.this).getHomoloGeneMap();
            Iterator iterator = hashmap.keySet().iterator();
            do {
                if(!iterator.hasNext())
                    break;
                String s = (String)iterator.next();
                Integer integer = (Integer)hashmap.get(s);
                if(!s.equals(getSpecies())) {
                    String s1 = "viewer.php?sp=" + s + "&lociId=" + integer.toString();
                    getAppletContext().showDocument(new URL(getBaseURL(), s1), "_blank");
                }
            } while(true);
        }
        catch(MalformedURLException malformedurlexception) { }
    }

    GenePane$HomoloGene() {
        super();
    }
}