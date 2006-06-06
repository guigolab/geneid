// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:07 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   Locus.java

import java.util.ArrayList;
import java.util.HashMap;

public class Locus {

    public static final int STRAND_UNKNOWN = 0;
    public static final int STRAND_DIRECT = 1;
    public static final int STRAND_COMPLEMENT = -1;
    private int _lociId;
    private String _chNo;
    private int _strand;
    private final ArrayList _defaultGeneList = new ArrayList();
    private ArrayList _geneList;
    private ArrayList _splicingPatternList;
    private int _minGenomePos;
    private int _maxGenomePos;
    private HashMap _homoloGeneMap;

    public Locus() {
        _geneList = new ArrayList();
        _splicingPatternList = new ArrayList();
        _homoloGeneMap = new HashMap();
    }

    public int directPos2complementPos(int i) {
        return SequenceUtils.directPos2complementPos(_minGenomePos, _maxGenomePos, i);
    }

    public int complementPos2directPos(int i) {
        return SequenceUtils.complementPos2directPos(_minGenomePos, _maxGenomePos, i);
    }

    public void sortGeneListBySplicingPattern(int i) {
        boolean aflag[] = new boolean[_defaultGeneList.size()];
        for(int j = 0; j < aflag.length; j++)
            aflag[j] = false;

        _geneList = new ArrayList();
        Alternative alternative = (Alternative)_splicingPatternList.get(i);
        int ai[] = alternative.getGeneIds1();
        for(int k = 0; k < ai.length; k++) {
            _geneList.add(_defaultGeneList.get(ai[k]));
            aflag[ai[k]] = true;
        }

        int ai1[] = alternative.getGeneIds2();
        for(int l = 0; l < ai1.length; l++) {
            _geneList.add(_defaultGeneList.get(ai1[l]));
            aflag[ai1[l]] = true;
        }

        for(int i1 = 0; i1 < aflag.length; i1++)
            if(!aflag[i1])
                _geneList.add(_defaultGeneList.get(i1));

    }

    public String getChnoAndStrandString() {
        StringBuffer stringbuffer = new StringBuffer("Chr.");
        stringbuffer.append(_chNo);
        stringbuffer.append(" (");
        if(_strand == 1)
            stringbuffer.append('+');
        else
        if(_strand == -1)
            stringbuffer.append('-');
        stringbuffer.append(')');
        return stringbuffer.toString();
    }

    public void addHomoloGene(String s, int i) {
        _homoloGeneMap.put(s, new Integer(i));
    }

    public int numHomoloGenes() {
        return _homoloGeneMap.size();
    }

    public Gene getGene(int i) {
        return (Gene)_defaultGeneList.get(i);
    }

    public void addGene(Gene gene) {
        _defaultGeneList.add(gene);
        _geneList.add(gene);
    }

    public void addSplicingPattern(Alternative alternative) {
        _splicingPatternList.add(alternative);
    }

    public int getGenomeLen() {
        return (_maxGenomePos - _minGenomePos) + 1;
    }

    public ArrayList getGeneList() {
        return _geneList;
    }

    public ArrayList getSplicingPatternList() {
        return _splicingPatternList;
    }

    public String getChNo() {
        return _chNo;
    }

    public int getLociId() {
        return _lociId;
    }

    public int getStrand() {
        return _strand;
    }

    public void setChNo(String s) {
        _chNo = s;
    }

    public void setLociId(int i) {
        _lociId = i;
    }

    public void setStrand(String s) {
        if(s.equals("+"))
            _strand = 1;
        else
        if(s.equals("-"))
            _strand = -1;
    }

    public int getMaxGenomePos() {
        return _maxGenomePos;
    }

    public int getMinGenomePos() {
        return _minGenomePos;
    }

    public void setMaxGenomePos(int i) {
        _maxGenomePos = i;
    }

    public void setMinGenomePos(int i) {
        _minGenomePos = i;
    }

    public HashMap getHomoloGeneMap() {
        return _homoloGeneMap;
    }
}