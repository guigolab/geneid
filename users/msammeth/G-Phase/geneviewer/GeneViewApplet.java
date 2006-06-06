// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:06 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   GeneViewApplet.java

import java.awt.Container;
import java.io.PrintStream;
import java.net.MalformedURLException;
import java.net.URL;
import javax.swing.JApplet;

public class GeneViewApplet extends JApplet {

    protected GenePane genePane;
    private static final String FETCH_LOCUS_FILE_NAME = "fetchlocus.php";

    public GeneViewApplet() {
    }

    public void init() {
        try {
            String s = getParameter("species");
            String s1 = getParameter("lociId");
            genePane = new GenePane();
            URL url = getDocumentBase();
            if(url.getProtocol().equals("file"))
                url = new URL("http://alterna.cbrc.jp/");
            if(url != null)
                genePane.createImageIcons(new AppletImageIconMaker(url));
            genePane.setBaseURL(url);
            genePane.setAppletContext(getAppletContext());
            genePane.init(s, getFetchLocusURL(url, "gff", s1, s), getFetchLocusURL(url, "xml", s1, s));
            getContentPane().add(genePane);
        }
        catch(Exception exception) {
            System.err.println(exception);
        }
    }

    private URL getFetchLocusURL(URL url, String s, String s1, String s2) throws MalformedURLException {
        return new URL(url, "fetchlocus.php?type=" + s + "&lociId=" + s1 + "&sp=" + s2);
    }
}