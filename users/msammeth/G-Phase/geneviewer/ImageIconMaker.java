// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:07 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   ImageIconMaker.java

import javax.swing.ImageIcon;

public interface ImageIconMaker {

    public static final String ZOOMIN_IMG = "zoomin.gif";
    public static final String ZOOMIN_DESC = "Zoom in";
    public static final String ZOOMOUT_IMG = "zoomout.gif";
    public static final String ZOOMOUT_DESC = "Zoom out";
    public static final String LEFT_IMG = "left.gif";
    public static final String LEFT_DESC = "Left";
    public static final String RIGHT_IMG = "right.gif";
    public static final String RIGHT_DESC = "Right";
    public static final String SEARCH_IMG = "search.gif";
    public static final String SEARCH_DESC = "Search Pattern";
    public static final String VIEW_IMG = "view.gif";
    public static final String VIEW_DESC = "Genes Overview";
    public static final String OPTION_IMG = "option.gif";
    public static final String OPTION_DESC = "Options";

    public abstract ImageIcon getZoominImageIcon();

    public abstract ImageIcon getZoomoutImageIcon();

    public abstract ImageIcon getLeftImageIcon();

    public abstract ImageIcon getRightImageIcon();

    public abstract ImageIcon getSearchImageIcon();

    public abstract ImageIcon getViewImageIcon();

    public abstract ImageIcon getOptionImageIcon();
}