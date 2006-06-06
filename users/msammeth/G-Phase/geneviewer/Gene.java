// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:04 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   Gene.java

import java.util.*;

public class Gene {

    private Locus _locus;
    private int _id;
    private String _cdnaId;
    private String _geneName;
    private int _cdsLeft;
    private int _cdsRight;
    private String _genBankAccession;
    private String _giNumber;
    private String _uniGeneId;
    private int _length;
    private boolean _nmd;
    private ArrayList _exonList;
    private Exon _equalViewExons[];
    private ArrayList _dataOrganizeAlterna;

    public Gene(Locus locus) {
        _exonList = new ArrayList();
        _equalViewExons = null;
        _dataOrganizeAlterna = new ArrayList();
        _locus = locus;
    }

    public boolean hasExon(int i, int j, boolean flag) {
        boolean flag1 = false;
        Iterator iterator = _exonList.iterator();
        do {
            if(!iterator.hasNext())
                break;
            Exon exon = (Exon)iterator.next();
            if(flag) {
                if(exon.getGenomeLeft() != i || exon.getGenomeRight() != j)
                    continue;
                flag1 = true;
                break;
            }
            Range range = new Range(exon.getGenomeLeft(), exon.getGenomeRight());
            if(!range.contains(new Range(i, j)))
                continue;
            flag1 = true;
            break;
        } while(true);
        return flag1;
    }

    public void setExonOrganizeAlternative(int i, int j, int k, int l) {
        TreeMap treemap = new TreeMap();
        treemap.put("AlternativeID", new Integer(i));
        treemap.put("StartExonID", new Integer(j));
        treemap.put("EndExonID", new Integer(k));
        treemap.put("WhichPattern", new Integer(l));
        _dataOrganizeAlterna.add(treemap);
    }

    public ArrayList getExonOrganizeAlternative() {
        return _dataOrganizeAlterna;
    }

    public boolean isAlterna(int i, int j, int k) {
        ArrayList arraylist = getExonOrganizeAlternative();
        for(int l = 0; l < arraylist.size(); l++) {
            TreeMap treemap = (TreeMap)arraylist.get(l);
            int i1 = ((Integer)treemap.get("AlternativeID")).intValue();
            int j1 = ((Integer)treemap.get("WhichPattern")).intValue();
            int k1 = ((Integer)treemap.get("StartExonID")).intValue();
            if(i == i1 && k == j1 && j == k1)
                return true;
        }

        return false;
    }

    public int getIndexEqualLeftRight(int i, int j) {
        int k = -1;
        for(int l = 0; l < _exonList.size(); l++) {
            int i1 = ((Exon)_exonList.get(l)).getGenomeLeft();
            int j1 = ((Exon)_exonList.get(l)).getGenomeRight();
            if(i1 == i && j1 == j)
                k = l;
        }

        return k;
    }

    public int getAlternaExon(int i, int j, int k) {
        ArrayList arraylist = getExonOrganizeAlternative();
        int l = -1;
        for(int i1 = 0; i1 < arraylist.size(); i1++) {
            TreeMap treemap = (TreeMap)arraylist.get(i1);
            int j1 = ((Integer)treemap.get("AlternativeID")).intValue();
            int k1 = ((Integer)treemap.get("WhichPattern")).intValue();
            int l1 = ((Integer)treemap.get("StartExonID")).intValue();
            if(j1 == i && k1 == j)
                l = l1 + k;
        }

        return l;
    }

    public int getExonIdIncludingRegion(int i, int j) {
        int k = -1;
        for(int l = 0; l < _exonList.size(); l++) {
            Exon exon = (Exon)_exonList.get(l);
            int i1 = exon.getGenomeLeft();
            int j1 = exon.getGenomeRight();
            if(i1 < j1) {
                if((i1 > i || j > j1) && (i > i1 || j1 > j))
                    continue;
                k = l;
                break;
            }
            if((i1 < i || j < j1) && (i < i1 || j1 < j))
                continue;
            k = l;
            break;
        }

        return k;
    }

    public Exon getExon(int i, int j) {
        Exon exon = null;
        Iterator iterator = _exonList.iterator();
        do {
            if(!iterator.hasNext())
                break;
            Exon exon1 = (Exon)iterator.next();
            if(exon1.getGenomeLeft() != i || exon1.getGenomeRight() != j)
                continue;
            exon = exon1;
            break;
        } while(true);
        return exon;
    }

    public void addExon(Exon exon) {
        _exonList.add(exon);
    }

    public String getCdnaId() {
        return _cdnaId;
    }

    public int getCdsLeft() {
        return _cdsLeft;
    }

    public int getCdsRight() {
        return _cdsRight;
    }

    public String getGeneName() {
        return _geneName;
    }

    public int getId() {
        return _id;
    }

    public void setCdnaId(String s) {
        _cdnaId = s;
    }

    public void setCdsLeft(int i) {
        _cdsLeft = i;
    }

    public void setCdsRight(int i) {
        _cdsRight = i;
    }

    public void setGeneName(String s) {
        _geneName = s;
    }

    public void setId(int i) {
        _id = i;
    }

    public String getGenBankAccession() {
        return _genBankAccession;
    }

    public String getGiNumber() {
        return _giNumber;
    }

    public int getLength() {
        return _length;
    }

    public String getUniGeneId() {
        return _uniGeneId;
    }

    public void setGenBankAccession(String s) {
        _genBankAccession = s;
    }

    public void setGiNumber(String s) {
        _giNumber = s;
    }

    public void setLength(int i) {
        _length = i;
    }

    public void setUniGeneId(String s) {
        _uniGeneId = s;
    }

    public void setNmd(boolean flag) {
        _nmd = flag;
    }

    public boolean isNmd() {
        return _nmd;
    }

    public Range getCdsRange() {
        return new Range(_cdsLeft, _cdsRight);
    }

    public Exon[] getExons() {
        return (Exon[])(Exon[])_exonList.toArray(new Exon[0]);
    }

    public Exon[] getEqualViewExons() {
        return _equalViewExons;
    }

    public void setEqualViewExons(Exon aexon[]) {
        _equalViewExons = aexon;
    }
}