// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:04 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   AppletImageIconMaker.java

import java.net.MalformedURLException;
import java.net.URL;
import javax.swing.ImageIcon;

public class AppletImageIconMaker
    implements ImageIconMaker {

    private URL baseUrl;

    public AppletImageIconMaker(URL url) {
        baseUrl = url;
    }

    public ImageIcon getLeftImageIcon() {
        return createImageIcon("left.gif", "Left");
    }

    public ImageIcon getRightImageIcon() {
        return createImageIcon("right.gif", "Right");
    }

    public ImageIcon getSearchImageIcon() {
        return createImageIcon("search.gif", "Search Pattern");
    }

    public ImageIcon getViewImageIcon() {
        return createImageIcon("view.gif", "Genes Overview");
    }

    public ImageIcon getZoominImageIcon() {
        return createImageIcon("zoomin.gif", "Zoom in");
    }

    public ImageIcon getZoomoutImageIcon() {
        return createImageIcon("zoomout.gif", "Zoom out");
    }

    public ImageIcon getOptionImageIcon() {
        return createImageIcon("option.gif", "Options");
    }

    private ImageIcon createImageIcon(String s, String s1) {
        ImageIcon imageicon = null;
        try {
            imageicon = new ImageIcon(new URL(baseUrl, s), s1);
        }
        catch(MalformedURLException malformedurlexception) { }
        return imageicon;
    }
}