package geneviewer;

// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:03 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   Alternative.java

import java.util.ArrayList;
import java.util.StringTokenizer;

public class Alternative {

    public static final String PATTERNIMG_FILE = "patternimg.php";
    public static final String RANGE_SEPARATOR = "..";
    public static final String EXON_SEPARATOR = " / ";
    public static final String ALTERNATYPE_A = "Alt. Donor site";
    public static final String ALTERNATYPE_B = "Alt. Acceptor site";
    public static final String ALTERNATYPE_C = "Cassette";
    public static final String ALTERNATYPE_D = "Mutually Exclusive";
    public static final String ALTERNATYPE_E1 = "Alt. poly(A) 2";
    public static final String ALTERNATYPE_E2 = "Alt. poly(A) 1";
    public static final String ALTERNATYPE_F = "Retained Intron";
    public static final String ALTERNATYPE_OTHERS = "Others";
    public static final String ALTERNATYPE_5END = "5'end";
    private int _id;
    private String _type;
    private int _genomeLeft;
    private int _genomeRight;
    private String _pattern1;
    private String _pattern2;
    private int _geneIds1[];
    private int _geneIds2[];
    private boolean _nagnag;
    private ArrayList exonList;
    private int startExonID1;
    private int startExonID2;
    private int numExon1;
    private int numExon2;

    public Alternative() {
        _geneIds1 = null;
        _geneIds2 = null;
        _pattern1 = null;
        _pattern2 = null;
        exonList = new ArrayList();
        startExonID1 = -1;
        startExonID2 = -1;
        numExon1 = -1;
        numExon2 = -1;
    }

    public boolean isAlternativeGene(int i) {
        boolean flag = false;
        int j = 0;
        do {
            if(j >= _geneIds1.length)
                break;
            if(_geneIds1[j] == i) {
                flag = true;
                break;
            }
            j++;
        } while(true);
        if(!flag) {
            for(int k = 0; k < _geneIds2.length; k++)
                if(_geneIds2[k] == i)
                    flag = true;

        }
        return flag;
    }

    public int getNumAlternativeGenes() {
        return _geneIds1.length + _geneIds2.length;
    }

    public void addExon(Exon exon) {
        exonList.add(exon);
    }

    public String getPatternImageFileName() {
        String s;
        String s1;
        if(_pattern1.compareTo(_pattern2) <= 0) {
            s = _pattern1;
            s1 = _pattern2;
        } else {
            s = _pattern2;
            s1 = _pattern1;
        }
        StringBuffer stringbuffer = new StringBuffer("patternimg.php");
        stringbuffer.append('?');
        stringbuffer.append("p1=");
        stringbuffer.append(s);
        stringbuffer.append("&p2=");
        stringbuffer.append(s1);
        return stringbuffer.toString();
    }

    public String getCaption() {
        String s = "";
        if(_type.equals("c"))
            s = "Cassette";
        else
        if(_type.equals("e2"))
            s = "Alt. poly(A) 1";
        else
        if(_type.equals("f"))
            s = "Retained Intron";
        else
        if(_type.equals("b"))
            s = "Alt. Acceptor site";
        else
        if(_type.equals("e1"))
            s = "Alt. poly(A) 2";
        else
        if(_type.equals("a"))
            s = "Alt. Donor site";
        else
        if(_type.equals("d"))
            s = "Mutually Exclusive";
        else
        if(_type.equals("5end"))
            s = "5'end";
        else
        if(_type.equals("others"))
            s = "Others";
        return s;
    }

    public String getExonsString() {
        StringBuffer stringbuffer = new StringBuffer();
        for(int i = 0; i < exonList.size(); i++) {
            stringbuffer.append(getExonRangeString(i));
            if(i != exonList.size() - 1)
                stringbuffer.append(" / ");
        }

        return stringbuffer.toString();
    }

    public String getExonRangeString(int i) {
        Exon exon = (Exon)exonList.get(i);
        return String.valueOf(exon.getGenomeLeft()) + ".." + String.valueOf(exon.getGenomeRight());
    }

    public String getLengthSubtypeString() {
        StringBuffer stringbuffer = new StringBuffer();
        for(int i = 0; i < exonList.size(); i++) {
            stringbuffer.append(getLengthSubtypeString(i));
            if(i != exonList.size() - 1)
                stringbuffer.append(" / ");
        }

        return stringbuffer.toString();
    }

    public String getLengthSubtypeString(int i) {
        StringBuffer stringbuffer = new StringBuffer();
        Exon exon = (Exon)exonList.get(i);
        int j = Math.abs(exon.getGenomeLeft() - exon.getGenomeRight()) + 1;
        int k = j % 3;
        stringbuffer.append("len=");
        stringbuffer.append(j);
        stringbuffer.append(" subtype=");
        stringbuffer.append(k);
        return stringbuffer.toString();
    }

    public int getAlternaLength() {
        return Math.abs(_genomeLeft - _genomeRight) + 1;
    }

    public ArrayList getExonList() {
        return exonList;
    }

    public int getId() {
        return _id;
    }

    public String getType() {
        return _type;
    }

    public void setExonList(ArrayList arraylist) {
        exonList = arraylist;
    }

    public void setId(int i) {
        _id = i;
    }

    public void setType(String s) {
        _type = s;
    }

    public String getPattern1() {
        return _pattern1;
    }

    public String getPattern2() {
        return _pattern2;
    }

    public void setGene1(String s) {
        ArrayList arraylist = new ArrayList();
        for(StringTokenizer stringtokenizer = new StringTokenizer(s, ":"); stringtokenizer.hasMoreTokens(); arraylist.add(stringtokenizer.nextToken()));
        _geneIds1 = new int[arraylist.size()];
        for(int i = 0; i < arraylist.size(); i++)
            _geneIds1[i] = Integer.parseInt((String)arraylist.get(i));

    }

    public void setGene2(String s) {
        ArrayList arraylist = new ArrayList();
        for(StringTokenizer stringtokenizer = new StringTokenizer(s, ":"); stringtokenizer.hasMoreTokens(); arraylist.add(stringtokenizer.nextToken()));
        _geneIds2 = new int[arraylist.size()];
        for(int i = 0; i < arraylist.size(); i++)
            _geneIds2[i] = Integer.parseInt((String)arraylist.get(i));

    }

    public void setPattern1(String s) {
        _pattern1 = s;
    }

    public void setPattern2(String s) {
        _pattern2 = s;
    }

    public int getGenomeLeft() {
        return _genomeLeft;
    }

    public int getGenomeRight() {
        return _genomeRight;
    }

    public void setGenomeLeft(int i) {
        _genomeLeft = i;
    }

    public void setGenomeRight(int i) {
        _genomeRight = i;
    }

    public void setStartExonID1(int i) {
        startExonID1 = i;
    }

    public void setStartExonID2(int i) {
        startExonID2 = i;
    }

    public int getStartExonID1() {
        return startExonID1;
    }

    public int getStartExonID2() {
        return startExonID2;
    }

    public void setNumExon1(int i) {
        i = calcNumExon(1);
        numExon1 = i;
    }

    public void setNumExon2(int i) {
        i = calcNumExon(2);
        numExon2 = i;
    }

    public int getNumExon1() {
        return numExon1;
    }

    public int getNumExon2() {
        return numExon2;
    }

    private int calcNumExon(int i) {
        String s = "";
        if(i == 1)
            s = getPattern1();
        else
            s = getPattern2();
        int j = 0;
        for(int k = 0; k < s.length(); k++)
            if(s.charAt(k) == '1' && (k == 0 || k > 0 && s.charAt(k - 1) == '0'))
                j++;

        return j;
    }

    public int[] getGeneIds1() {
        return _geneIds1;
    }

    public int[] getGeneIds2() {
        return _geneIds2;
    }

    public void setNagnag(boolean flag) {
        _nagnag = flag;
    }

    public boolean isNagnag() {
        return _nagnag;
    }
}