// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:07 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   GffRecord.java

import java.io.PrintStream;
import java.util.StringTokenizer;
import java.util.TreeMap;

public class GffRecord
    implements Cloneable {

    public static final String FEATURE_ORF = "ORF";
    public static final String FEATURE_UTR = "UTR";
    public static final char STRAND_DIRECT = 43;
    public static final char STRAND_COMPLEMENT = 45;
    private String seqname;
    private String source;
    private String feature;
    private int start;
    private int end;
    private double score;
    private int strand;
    private int frame;
    private String label;
    private TreeMap attribute;

    public GffRecord() {
        seqname = null;
        source = null;
        feature = null;
        start = -1;
        end = -1;
        score = -1D;
        strand = -1;
        frame = -1;
        label = null;
        attribute = null;
    }

    public GffRecord(String s, String s1, String s2, int i, int j, double d, 
            int k, int l, String s3, TreeMap treemap) {
        seqname = null;
        source = null;
        feature = null;
        start = -1;
        end = -1;
        score = -1D;
        strand = -1;
        frame = -1;
        label = null;
        attribute = null;
        seqname = s;
        source = s1;
        feature = s2;
        start = i;
        end = j;
        score = d;
        strand = k;
        frame = l;
        label = s3;
        attribute = treemap;
    }

    public static GffRecord parseGFF(String s) {
        StringTokenizer stringtokenizer = new StringTokenizer(s, "\t");
        if(stringtokenizer.countTokens() < 8) {
            System.err.println("GFF: less than 8 fields in line:" + s);
            return null;
        }
        String s1 = stringtokenizer.nextToken();
        String s2 = stringtokenizer.nextToken();
        String s3 = stringtokenizer.nextToken();
        int i = -1;
        int j = -1;
        byte byte0 = -1;
        int k = -1;
        double d = -1D;
        try {
            i = Integer.parseInt(stringtokenizer.nextToken());
        }
        catch(NumberFormatException numberformatexception) { }
        try {
            j = Integer.parseInt(stringtokenizer.nextToken());
        }
        catch(NumberFormatException numberformatexception1) { }
        try {
            d = Double.parseDouble(stringtokenizer.nextToken());
        }
        catch(NumberFormatException numberformatexception2) { }
        String s4 = stringtokenizer.nextToken();
        if(s4.equals("+"))
            byte0 = 0;
        else
        if(s4.equals("-"))
            byte0 = 1;
        else
        if(s4.equals("."))
            byte0 = 2;
        else
            byte0 = 3;
        try {
            k = Integer.parseInt(stringtokenizer.nextToken());
        }
        catch(NumberFormatException numberformatexception3) {
            k = 3;
        }
        String s5 = null;
        if(stringtokenizer.hasMoreTokens())
            s5 = stringtokenizer.nextToken();
        TreeMap treemap = null;
        if(stringtokenizer.hasMoreTokens()) {
            String s6 = stringtokenizer.nextToken();
            treemap = new TreeMap();
            stringtokenizer = new StringTokenizer(s6, "; ");
            do {
                if(!stringtokenizer.hasMoreTokens())
                    break;
                String s7 = stringtokenizer.nextToken();
                if(stringtokenizer.countTokens() == 0)
                    break;
                String s8 = stringtokenizer.nextToken();
                treemap.put(s7, s8);
            } while(true);
        }
        return new GffRecord(s1, s2, s3, i, j, d, byte0, k, s5, treemap);
    }

    public String toString() {
        StringBuffer stringbuffer = new StringBuffer();
        stringbuffer.append("'" + getSource() + "' ");
        stringbuffer.append(getFeature() + " ");
        stringbuffer.append(getStart() + "-" + getEnd() + " (length " + getLength() + ") ");
        return stringbuffer.toString();
    }

    public Object clone() {
        return super.clone();
        CloneNotSupportedException clonenotsupportedexception;
        clonenotsupportedexception;
        return new GffRecord();
    }

    public int getLength() {
        return Math.abs(getEnd() - getStart()) + 1;
    }

    public TreeMap getAttribute() {
        return attribute;
    }

    public int getEnd() {
        return end;
    }

    public String getFeature() {
        return feature;
    }

    public int getFrame() {
        return frame;
    }

    public String getLabel() {
        return label;
    }

    public double getScore() {
        return score;
    }

    public String getSeqname() {
        return seqname;
    }

    public String getSource() {
        return source;
    }

    public int getStart() {
        return start;
    }

    public int getStrand() {
        return strand;
    }

    public void setAttribute(TreeMap treemap) {
        attribute = treemap;
    }

    public void setEnd(int i) {
        end = i;
    }

    public void setFeature(String s) {
        feature = s;
    }

    public void setFrame(int i) {
        frame = i;
    }

    public void setLabel(String s) {
        label = s;
    }

    public void setScore(double d) {
        score = d;
    }

    public void setSeqname(String s) {
        seqname = s;
    }

    public void setSource(String s) {
        source = s;
    }

    public void setStart(int i) {
        start = i;
    }

    public void setStrand(int i) {
        strand = i;
    }
}