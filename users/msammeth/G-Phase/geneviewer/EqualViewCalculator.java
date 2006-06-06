// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:04 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   EqualViewCalculator.java

import java.util.*;

public class EqualViewCalculator {

    private ArrayList _alternativeList;
    private ArrayList _geneList;
    private int _strand;
    private ArrayList _alternaLineList;
    private ArrayList _lineInAlternaRange;
    private TreeSet _lineNoAlternaRange;
    private TreeMap _genomepos2linenoMap;
    private ArrayList _exonNoAlternaRange;
    private String _leftORright;
    private int _changeExonId;
    private int _changeExonPos;

    EqualViewCalculator(ArrayList arraylist, ArrayList arraylist1, int i) {
        _alternativeList = null;
        _geneList = null;
        _strand = 0;
        _alternaLineList = null;
        _lineInAlternaRange = null;
        _lineNoAlternaRange = null;
        _genomepos2linenoMap = null;
        _exonNoAlternaRange = null;
        _leftORright = null;
        _changeExonId = -1;
        _changeExonPos = -1;
        _alternativeList = arraylist;
        _geneList = arraylist1;
        _strand = i;
        _alternaLineList = new ArrayList();
        _lineInAlternaRange = new ArrayList();
        _lineNoAlternaRange = new TreeSet();
        _exonNoAlternaRange = new ArrayList();
    }

    public void calculate() {
        try {
            for(int i = 0; i < _alternativeList.size(); i++)
                setAlternaLine(i);

            for(int j = 0; j < _geneList.size(); j++) {
                Gene gene = (Gene)_geneList.get(j);
                Exon aexon[] = gene.getExons();
                Exon aexon1[] = new Exon[aexon.length];
                for(int i1 = 0; i1 < aexon.length; i1++) {
                    Exon exon = calculateExon(aexon[i1], j, i1);
                    aexon1[i1] = exon;
                }

                gene.setEqualViewExons(aexon1);
            }

            for(int k = 0; k < _alternativeList.size(); k++)
                cutEdge(k);

            for(int l = 0; l < _geneList.size(); l++) {
                setLineInAlternaRange(l, "left");
                setLineInAlternaRange(l, "right");
            }

            setLineNoAlternaRange();
            setExonNoAlternaRange();
            setLineList();
        }
        catch(Exception exception) {
            exception.printStackTrace();
        }
    }

    private Exon calculateExon(Exon exon, int i, int j) throws CloneNotSupportedException {
        Exon exon1 = (Exon)exon.clone();
        if(isAlternaRange(i, j)) {
            for(int k = 0; k < _alternativeList.size(); k++) {
                if(isIntersectAlternaFirstExon(k, i, j, 1) && isAlterna(k, i, j, 1)) {
                    setDataOrganizeAlterna(k, i, j, 1);
                    continue;
                }
                if(isIntersectAlternaFirstExon(k, i, j, 2) && isAlterna(k, i, j, 2))
                    setDataOrganizeAlterna(k, i, j, 2);
            }

            int l;
            for(l = 0; l < _alternativeList.size() && !isAlternaExon(l, i, j); l++);
            if(l == _alternativeList.size()) {
                addLineInAlternaRange(exon.getGenomeLeft());
                addLineInAlternaRange(exon.getGenomeRight());
                exon.setEqualViewLeft(exon.getGenomeLeft());
                exon.setEqualViewRight(exon.getGenomeRight());
                exon1.setGenomeLeft(exon.getGenomeLeft());
                exon1.setGenomeRight(exon.getGenomeRight());
            }
        } else {
            addExonNoAlterna(exon.getGenomeLeft(), exon.getGenomeRight());
        }
        return exon1;
    }

    private void setExonNoAlternaRange() {
        int i = 0;
        do {
            if(i >= _geneList.size())
                break;
            Exon aexon[] = ((Gene)_geneList.get(i)).getExons();
            for(int j = 0; j < aexon.length; j++) {
                if(isAlternaRange(i, j))
                    continue;
                Exon exon = aexon[j];
                int k = exon.getGenomeLeft();
                int l = exon.getGenomeRight();
                if(getStrand() == 1) {
                    Iterator iterator = _lineNoAlternaRange.iterator();
                    int i1 = -1;
                    int k1 = -1;
                    do {
                        if(!iterator.hasNext())
                            break;
                        k1 = ((Integer)iterator.next()).intValue();
                        if(k1 > k)
                            break;
                        i1 = k1;
                    } while(true);
                    if(i1 == -1)
                        i1 = k1;
                    exon.setEqualViewLeft(i1);
                    iterator = _lineNoAlternaRange.iterator();
                    i1 = -1;
                    do {
                        if(!iterator.hasNext())
                            break;
                        int j2 = ((Integer)iterator.next()).intValue();
                        if(j2 < l)
                            continue;
                        i1 = j2;
                        break;
                    } while(true);
                    exon.setEqualViewRight(i1);
                    continue;
                }
                if(getStrand() != -1)
                    continue;
                Iterator iterator1 = _lineNoAlternaRange.iterator();
                int j1 = -1;
                do {
                    if(!iterator1.hasNext())
                        break;
                    int l1 = ((Integer)iterator1.next()).intValue();
                    if(k > l1)
                        continue;
                    j1 = l1;
                    break;
                } while(true);
                exon.setEqualViewLeft(j1);
                iterator1 = _lineNoAlternaRange.iterator();
                j1 = -1;
                int i2 = -1;
                do {
                    if(!iterator1.hasNext())
                        break;
                    i2 = ((Integer)iterator1.next()).intValue();
                    if(l < i2)
                        break;
                    j1 = i2;
                } while(true);
                if(j1 == -1)
                    j1 = i2;
                exon.setEqualViewRight(j1);
            }

            i++;
        } while(true);
    }

    private void setLineInAlternaRange(int i, String s) throws Exception {
        Gene gene = (Gene)_geneList.get(i);
        ArrayList arraylist = gene.getExonOrganizeAlternative();
        for(int j = 0; j < arraylist.size(); j++) {
            TreeMap treemap = (TreeMap)arraylist.get(j);
            int k = ((Integer)treemap.get("AlternativeID")).intValue();
            int l = ((Integer)treemap.get("StartExonID")).intValue();
            int i1 = ((Integer)treemap.get("EndExonID")).intValue();
            int j1 = ((Integer)treemap.get("WhichPattern")).intValue();
            int k1 = -1;
            int l1 = -1;
            String s1 = "";
            if(j1 == 1) {
                k1 = ((Alternative)_alternativeList.get(k)).getGeneIds1()[0];
                l1 = ((Alternative)_alternativeList.get(k)).getStartExonID1();
                s1 = ((Alternative)_alternativeList.get(k)).getPattern1();
            } else {
                k1 = ((Alternative)_alternativeList.get(k)).getGeneIds2()[0];
                l1 = ((Alternative)_alternativeList.get(k)).getStartExonID2();
                s1 = ((Alternative)_alternativeList.get(k)).getPattern2();
            }
            Exon aexon[] = ((Gene)_geneList.get(i)).getExons();
            Exon aexon1[] = ((Gene)_geneList.get(k1)).getExons();
            int i2 = l;
            for(int j2 = l1; i2 <= i1; j2++) {
                Exon exon = aexon[i2];
                Exon exon1 = aexon1[j2];
                int k2 = -1;
                if(s1.charAt(0) == '1' && i2 == l)
                    try {
                        k2 = getAlternaRangeEdge(k, "left");
                    }
                    catch(Exception exception) {
                        throw new Exception(exception);
                    }
                else
                    k2 = exon1.getGenomeLeft();
                if(s.equals("left") && exon.getGenomeLeft() != k2) {
                    setChangeExonInfo("left", i2, k2);
                    if(!isEqualAlternaTerm(i)) {
                        setChangeExonInfo(null, -1, -1);
                        addLineInAlternaRange(exon.getGenomeLeft());
                        exon.setEqualViewLeft(exon.getGenomeLeft());
                    } else {
                        exon.setEqualViewLeft(k2);
                    }
                } else
                if(s.equals("left") && exon.getGenomeLeft() == k2)
                    exon.setEqualViewLeft(exon.getGenomeLeft());
                int l2 = -1;
                if(s1.charAt(s1.length() - 1) == '1' && i2 == i1)
                    try {
                        l2 = getAlternaRangeEdge(k, "right");
                    }
                    catch(Exception exception1) {
                        throw new Exception(exception1);
                    }
                else
                    l2 = exon1.getGenomeRight();
                if(s.equals("right") && exon.getGenomeRight() != l2) {
                    setChangeExonInfo("right", i2, l2);
                    if(!isEqualAlternaTerm(i)) {
                        setChangeExonInfo(null, -1, -1);
                        addLineInAlternaRange(exon.getGenomeRight());
                        exon.setEqualViewRight(exon.getGenomeRight());
                    } else {
                        exon.setEqualViewRight(l2);
                    }
                } else
                if(s.equals("right") && exon.getGenomeRight() == l2)
                    exon.setEqualViewRight(exon.getGenomeRight());
                i2++;
            }

        }

    }

    private void cutEdge(int i) {
        Alternative alternative = (Alternative)_alternativeList.get(i);
        int ai[] = alternative.getGeneIds1();
        int ai1[] = alternative.getGeneIds2();
        String s = alternative.getPattern1();
        String s1 = alternative.getPattern2();
        int j = alternative.getStartExonID1();
        int k = alternative.getStartExonID2();
        int l = alternative.getNumExon1();
        int i1 = alternative.getNumExon2();
        Exon aexon[] = ((Gene)_geneList.get(ai[0])).getExons();
        Exon aexon1[] = ((Gene)_geneList.get(ai1[0])).getExons();
        if(j == -1 || k == -1 || aexon.length == j + 1 || aexon1.length == k + 1)
            return;
        if(s.charAt(0) == '1' && s1.charAt(0) == '1') {
            int j1 = aexon[j].getGenomeLeft();
            int l1 = aexon1[k].getGenomeLeft();
            if(j1 != l1) {
                int j2 = -1;
                if(getStrand() == 1)
                    j2 = Math.max(j1, l1);
                else
                    j2 = Math.min(j1, l1);
                int l2 = -1;
                int j3 = -1;
                int l3 = -1;
                Exon aexon2[];
                if(j1 != j2) {
                    l2 = ai[0];
                    j3 = j1;
                    l3 = j;
                    aexon2 = aexon;
                } else {
                    l2 = ai1[0];
                    j3 = l1;
                    l3 = k;
                    aexon2 = aexon1;
                }
                setChangeExonInfo("left", l3, j2);
                if(!isEqualAlternaTerm(l2)) {
                    addLineInAlternaRange(j3);
                    aexon2[l3].setEqualViewLeft(j3);
                } else {
                    aexon2[l3].setEqualViewLeft(j2);
                }
                setChangeExonInfo(null, -1, -1);
            }
        }
        if(s.charAt(s.length() - 1) == '1' && s1.charAt(s1.length() - 1) == '1') {
            int k1 = aexon[(j + l) - 1].getGenomeRight();
            int i2 = aexon1[(k + i1) - 1].getGenomeRight();
            if(k1 != i2) {
                int k2 = -1;
                if(getStrand() == 1)
                    k2 = Math.min(k1, i2);
                else
                    k2 = Math.max(k1, i2);
                int i3 = -1;
                int k3 = -1;
                int i4 = -1;
                int j4 = -1;
                Exon aexon3[];
                if(k1 != k2) {
                    i3 = ai[0];
                    k3 = k1;
                    i4 = j;
                    aexon3 = aexon;
                    j4 = l;
                } else {
                    i3 = ai1[0];
                    k3 = i2;
                    i4 = k;
                    aexon3 = aexon1;
                    j4 = i1;
                }
                setChangeExonInfo("right", (i4 + j4) - 1, k2);
                if(!isEqualAlternaTerm(i3)) {
                    addLineInAlternaRange(k2);
                    aexon3[(i4 + j4) - 1].setEqualViewRight(k3);
                } else {
                    aexon3[(i4 + j4) - 1].setEqualViewRight(k2);
                }
                setChangeExonInfo(null, -1, -1);
            }
        }
    }

    private boolean isEqualAlternaTerm(int i) {
        int j = 0;
        do {
            if(j >= _alternativeList.size())
                break;
            int k;
            for(k = 0; k < ((Gene)_geneList.get(i)).getExons().length && (!isIntersectAlternaFirstExon(j, i, k, 1) || isAlterna(j, i, k, 1) == ((Gene)_geneList.get(i)).isAlterna(j, k, 1) && isAlterna(j, i, k, 2) == ((Gene)_geneList.get(i)).isAlterna(j, k, 2)) && (!isIntersectAlternaFirstExon(j, i, k, 2) || isAlterna(j, i, k, 1) == ((Gene)_geneList.get(i)).isAlterna(j, k, 1) && isAlterna(j, i, k, 2) == ((Gene)_geneList.get(i)).isAlterna(j, k, 2)); k++);
            if(k != ((Gene)_geneList.get(i)).getExons().length)
                break;
            j++;
        } while(true);
        return j == _alternativeList.size();
    }

    private void setChangeExonInfo(String s, int i, int j) {
        _leftORright = s;
        _changeExonId = i;
        _changeExonPos = j;
    }

    private int getAlternaRangeEdge(int i, String s) {
        if(s.equals("left"))
            return ((Integer)((ArrayList)_alternaLineList.get(i)).get(0)).intValue();
        if(s.equals("right"))
            return ((Integer)((ArrayList)_alternaLineList.get(i)).get(((ArrayList)_alternaLineList.get(i)).size() - 1)).intValue();
        break MISSING_BLOCK_LABEL_81;
        Exception exception;
        exception;
        return 0;
    }

    protected boolean isAlternaRange(int i, int j) {
        Exon aexon[] = ((Gene)_geneList.get(i)).getExons();
        Exon exon = aexon[j];
        int k = exon.getGenomeLeft();
        int l = exon.getGenomeRight();
        for(int i1 = 0; i1 < _alternativeList.size(); i1++) {
            ArrayList arraylist = (ArrayList)_alternaLineList.get(i1);
            if(arraylist.size() == 0)
                continue;
            int j1 = ((Integer)((ArrayList)_alternaLineList.get(i1)).get(0)).intValue();
            int k1 = ((Integer)((ArrayList)_alternaLineList.get(i1)).get(((ArrayList)_alternaLineList.get(i1)).size() - 1)).intValue();
            if(getStrand() == 1) {
                if(l >= j1 && k1 >= k)
                    return true;
                continue;
            }
            if(getStrand() == -1 && l <= j1 && k1 <= k)
                return true;
        }

        return false;
    }

    private void addExonNoAlterna(int i, int j) {
        int k = 0;
        int l = 0;
        do {
            if(l >= _exonNoAlternaRange.size())
                break;
            TreeMap treemap1 = (TreeMap)_exonNoAlternaRange.get(l);
            int i1 = ((Integer)treemap1.get("left")).intValue();
            int j1 = ((Integer)treemap1.get("right")).intValue();
            if(i == i1 && j == j1) {
                k = -1;
                break;
            }
            if(i < i1) {
                k = l;
                break;
            }
            k = l + 1;
            l++;
        } while(true);
        if(k != -1) {
            TreeMap treemap = new TreeMap();
            treemap.put("left", new Integer(i));
            treemap.put("right", new Integer(j));
            _exonNoAlternaRange.add(k, treemap);
        }
    }

    private void setLineNoAlternaRange() {
        if(getStrand() == 1) {
            int k2;
            for(int i = 0; i < _exonNoAlternaRange.size(); i = k2) {
                TreeMap treemap = (TreeMap)_exonNoAlternaRange.get(i);
                int k = ((Integer)treemap.get("left")).intValue();
                if(i == 0)
                    _lineNoAlternaRange.add(new Integer(k));
                int i1 = ((Integer)treemap.get("right")).intValue();
                int k1 = -1;
                int i2 = i1;
                k2 = i + 1;
                int i3 = -1;
                TreeSet treeset = new TreeSet();
                treeset.add(new Integer(i1));
                do {
                    if(k2 >= _exonNoAlternaRange.size())
                        break;
                    TreeMap treemap2 = (TreeMap)_exonNoAlternaRange.get(k2);
                    i3 = ((Integer)treemap2.get("left")).intValue();
                    int k3 = ((Integer)treemap2.get("right")).intValue();
                    if(i1 < i3) {
                        k1 = i3;
                        break;
                    }
                    treeset.add(new Integer(k3));
                    k2++;
                } while(true);
                Iterator iterator = treeset.iterator();
                do {
                    if(!iterator.hasNext())
                        break;
                    int l3 = ((Integer)iterator.next()).intValue();
                    if(k1 == -1) {
                        i2 = ((Integer)treeset.last()).intValue();
                        break;
                    }
                    if(i3 <= l3)
                        break;
                    i2 = l3;
                } while(true);
                if(k1 != -1 && !_lineNoAlternaRange.contains(new Integer(k1)))
                    _lineNoAlternaRange.add(new Integer(k1));
                if(!_lineNoAlternaRange.contains(new Integer(i2)))
                    _lineNoAlternaRange.add(new Integer(i2));
            }

        } else
        if(getStrand() == -1) {
            int l2;
            for(int j = _exonNoAlternaRange.size() - 1; j >= 0; j = l2) {
                TreeMap treemap1 = (TreeMap)_exonNoAlternaRange.get(j);
                int l = ((Integer)treemap1.get("left")).intValue();
                if(j == _exonNoAlternaRange.size() - 1)
                    _lineNoAlternaRange.add(new Integer(l));
                int j1 = ((Integer)treemap1.get("right")).intValue();
                int l1 = -1;
                int j2 = j1;
                l2 = j - 1;
                int j3 = -1;
                TreeSet treeset1 = new TreeSet();
                treeset1.add(new Integer(j1));
                do {
                    if(l2 < 0)
                        break;
                    TreeMap treemap3 = (TreeMap)_exonNoAlternaRange.get(l2);
                    j3 = ((Integer)treemap3.get("left")).intValue();
                    int i4 = ((Integer)treemap3.get("right")).intValue();
                    if(j1 > j3) {
                        l1 = j3;
                        break;
                    }
                    treeset1.add(new Integer(i4));
                    l2--;
                } while(true);
                Iterator iterator1 = treeset1.iterator();
                do {
                    if(!iterator1.hasNext())
                        break;
                    int j4 = ((Integer)iterator1.next()).intValue();
                    if(l1 == -1) {
                        j2 = ((Integer)treeset1.first()).intValue();
                        break;
                    }
                    if(j3 >= j4)
                        break;
                    j2 = j4;
                } while(true);
                if(l1 != -1)
                    _lineNoAlternaRange.add(new Integer(l1));
                _lineNoAlternaRange.add(new Integer(j2));
            }

        }
    }

    private void addAlternaLine(int i, int j, int k) {
        if(!((ArrayList)_alternaLineList.get(i)).contains(new Integer(j)))
            ((ArrayList)_alternaLineList.get(i)).add(new Integer(j));
        if(!((ArrayList)_alternaLineList.get(i)).contains(new Integer(k)))
            ((ArrayList)_alternaLineList.get(i)).add(new Integer(k));
    }

    private void addLineInAlternaRange(int i) {
        if(!_lineInAlternaRange.contains(new Integer(i)))
            _lineInAlternaRange.add(new Integer(i));
    }

    private boolean isIntersectAlternaFirstExon(int i, int j, int k, int l) {
        Alternative alternative = (Alternative)_alternativeList.get(i);
        int i1 = 0;
        Object obj = null;
        int j1 = 0;
        if(l == 1) {
            i1 = alternative.getGeneIds1()[0];
            String s = alternative.getPattern1();
            j1 = alternative.getStartExonID1();
        } else
        if(l == 2) {
            i1 = alternative.getGeneIds2()[0];
            String s1 = alternative.getPattern2();
            j1 = alternative.getStartExonID2();
        }
        if(j1 == -1)
            return false;
        Exon aexon[] = ((Gene)_geneList.get(i1)).getExons();
        int k1 = aexon[j1].getGenomeLeft();
        int l1 = aexon[j1].getGenomeRight();
        Exon aexon1[] = ((Gene)_geneList.get(j)).getExons();
        int i2 = aexon1[k].getGenomeLeft();
        int j2 = aexon1[k].getGenomeRight();
        if(_leftORright != null && _changeExonId == k)
            if(_leftORright.equals("left"))
                i2 = _changeExonPos;
            else
                j2 = _changeExonPos;
        if(getStrand() == 1)
            return j2 >= k1 && l1 >= i2;
        if(getStrand() == -1)
            return j2 <= k1 && l1 <= i2;
        else
            return false;
    }

    private void setDataOrganizeAlterna(int i, int j, int k, int l) {
        int i1 = 0;
        if(l == 1)
            i1 = ((Alternative)_alternativeList.get(i)).getNumExon1();
        else
        if(l == 2)
            i1 = ((Alternative)_alternativeList.get(i)).getNumExon2();
        ((Gene)_geneList.get(j)).setExonOrganizeAlternative(i, k, (k + i1) - 1, l);
    }

    private boolean isAlternaExon(int i, int j, int k) {
        Gene gene = (Gene)_geneList.get(j);
        ArrayList arraylist = gene.getExonOrganizeAlternative();
        for(int l = 0; l < arraylist.size(); l++) {
            TreeMap treemap = (TreeMap)arraylist.get(l);
            int i1 = ((Integer)treemap.get("AlternativeID")).intValue();
            int j1 = ((Integer)treemap.get("StartExonID")).intValue();
            int k1 = ((Integer)treemap.get("EndExonID")).intValue();
            if(i1 == i && j1 <= k && k <= k1)
                return true;
        }

        return false;
    }

    private int getStrand() {
        return _strand;
    }

    private boolean isAlterna(int i, int j, int k, int l) {
        int i1 = ((Alternative)_alternativeList.get(i)).getGeneIds1()[0];
        int j1 = ((Alternative)_alternativeList.get(i)).getGeneIds2()[0];
        String s = ((Alternative)_alternativeList.get(i)).getPattern1();
        String s1 = ((Alternative)_alternativeList.get(i)).getPattern2();
        int k1 = ((Alternative)_alternativeList.get(i)).getStartExonID1();
        int l1 = ((Alternative)_alternativeList.get(i)).getStartExonID2();
        if(k1 == -1 || l1 == -1)
            return false;
        Exon aexon[] = ((Gene)_geneList.get(i1)).getExons();
        Exon aexon1[] = ((Gene)_geneList.get(j1)).getExons();
        Exon aexon2[] = ((Gene)_geneList.get(j)).getExons();
        String s2 = s;
        String s3 = s1;
        int i2 = k1;
        int j2 = l1;
        int k2 = k;
        int l2 = 0;
        if(l == 1) {
            String s4 = s2;
            s2 = s3;
            s3 = s4;
            Exon aexon3[] = aexon;
            aexon = aexon1;
            aexon1 = aexon3;
            int k3 = i2;
            i2 = j2;
            j2 = k3;
        }
        for(l2 = 0; l2 < s2.length() - 1; l2++) {
            if(k2 > aexon2.length - 1) {
                String s5 = null;
                if(l == 1)
                    s5 = s;
                else
                    s5 = s1;
                for(; l2 < s5.length() && s5.charAt(l2) == '0'; l2++);
                return l2 == s5.length();
            }
            int i3 = aexon2[k2].getGenomeRight();
            int j3 = aexon2[k2].getGenomeLeft();
            if(_leftORright != null && _changeExonId == k2)
                if(_leftORright.equals("left"))
                    j3 = _changeExonPos;
                else
                    i3 = _changeExonPos;
            if(s2.charAt(l2) == '1' && s2.charAt(l2 + 1) == '0' && s3.charAt(l2) == '1' && s3.charAt(l2 + 1) == '1') {
                int l3 = aexon[i2].getGenomeRight();
                if(getStrand() == 1) {
                    if(j3 >= l3 || l3 >= i3)
                        break;
                    i2++;
                    continue;
                }
                if(getStrand() != -1)
                    continue;
                if(j3 <= l3 || l3 <= i3)
                    break;
                i2++;
                continue;
            }
            if(s2.charAt(l2) == '1' && s2.charAt(l2 + 1) == '1' && s3.charAt(l2) == '1' && s3.charAt(l2 + 1) == '0') {
                int i4 = aexon1[j2].getGenomeRight();
                int l6 = aexon[i2].getGenomeLeft();
                if(getStrand() == 1) {
                    if(l6 >= i3 || i3 > i4)
                        break;
                    k2++;
                    j2++;
                    continue;
                }
                if(getStrand() != -1)
                    continue;
                if(l6 <= i3 || i3 < i4)
                    break;
                k2++;
                j2++;
                continue;
            }
            if(s2.charAt(l2) == '1' && s2.charAt(l2 + 1) == '0' && s3.charAt(l2) == '1' && s3.charAt(l2 + 1) == '0') {
                int j4 = aexon[i2].getGenomeRight();
                if(j4 != i3)
                    break;
                j2++;
                k2++;
                i2++;
                continue;
            }
            if(s2.charAt(l2) == '1' && s2.charAt(l2 + 1) == '1' && s3.charAt(l2) == '0' && s3.charAt(l2 + 1) == '1') {
                int k4 = aexon1[j2].getGenomeLeft();
                int i7 = aexon[i2].getGenomeRight();
                if(getStrand() != 1 ? getStrand() == -1 && (k4 < j3 || j3 <= i7) : k4 > j3 || j3 >= i7)
                    break;
                continue;
            }
            if(s2.charAt(l2) == '1' && s2.charAt(l2 + 1) == '0' && s3.charAt(l2) == '0' && s3.charAt(l2 + 1) == '1') {
                int l4 = aexon[i2].getGenomeRight();
                if(l4 != j3)
                    break;
                i2++;
                continue;
            }
            if(s2.charAt(l2) == '1' && s2.charAt(l2 + 1) == '0' && s3.charAt(l2) == '0' && s3.charAt(l2 + 1) == '0') {
                int i5 = aexon[i2].getGenomeRight();
                int j7 = aexon[i2].getGenomeLeft();
                int i8 = -1;
                if(k2 != 0)
                    i8 = aexon2[k2 - 1].getGenomeRight();
                if(getStrand() == 1) {
                    if(i5 >= j3 || (l2 == 0 || s2.charAt(l2 - 1) != '1' || s3.charAt(l2 - 1) != '1') && k2 != 0 && i8 > j7)
                        break;
                    i2++;
                    continue;
                }
                if(getStrand() != -1)
                    continue;
                if(i5 <= j3 || (l2 == 0 || s2.charAt(l2 - 1) != '1' || s3.charAt(l2 - 1) != '1') && k2 != 0 && i8 < j7)
                    break;
                i2++;
                continue;
            }
            if(s2.charAt(l2) == '0' && s2.charAt(l2 + 1) == '1' && s3.charAt(l2) == '1' && s3.charAt(l2 + 1) == '1') {
                int j5 = aexon[i2].getGenomeLeft();
                if(getStrand() != 1 ? getStrand() == -1 && (j3 <= j5 || j5 <= i3) : j3 >= j5 || j5 >= i3)
                    break;
                continue;
            }
            if(s2.charAt(l2) == '0' && s2.charAt(l2 + 1) == '1' && s3.charAt(l2) == '1' && s3.charAt(l2 + 1) == '0') {
                int k5 = aexon[i2].getGenomeLeft();
                if(k5 != i3)
                    break;
                j2++;
                k2++;
                continue;
            }
            if(s2.charAt(l2) == '0' && s2.charAt(l2 + 1) == '0' && s3.charAt(l2) == '1' && s3.charAt(l2 + 1) == '0') {
                int l5 = -1;
                if(i2 != 0)
                    l5 = aexon[i2 - 1].getGenomeRight();
                int k7 = aexon[i2].getGenomeLeft();
                int j8 = aexon1[j2].getGenomeLeft();
                if(getStrand() != 1 ? getStrand() == -1 && (i3 <= k7 || (l2 == 0 || s2.charAt(l2 - 1) != '1' || s3.charAt(l2 - 1) != '1') && i2 != 0 && l5 < j3) : i3 >= k7 || (l2 == 0 || s2.charAt(l2 - 1) != '1' || s3.charAt(l2 - 1) != '1') && i2 != 0 && l5 > j3)
                    break;
                j2++;
                k2++;
                continue;
            }
            if(s2.charAt(l2) == '0' && s2.charAt(l2 + 1) == '1' && s3.charAt(l2) == '0' && s3.charAt(l2 + 1) == '1') {
                int i6 = aexon[i2].getGenomeLeft();
                if(i6 != j3)
                    break;
                continue;
            }
            if(s2.charAt(l2) == '0' && s2.charAt(l2 + 1) == '1' && s3.charAt(l2) == '0' && s3.charAt(l2 + 1) == '0') {
                int j6 = aexon[i2].getGenomeLeft();
                if(getStrand() != 1 ? getStrand() == -1 && j6 <= j3 : j6 >= j3)
                    break;
                continue;
            }
            if(s2.charAt(l2) != '0' || s2.charAt(l2 + 1) != '0' || s3.charAt(l2) != '0' || s3.charAt(l2 + 1) != '1')
                continue;
            int k6 = aexon1[j2].getGenomeLeft();
            int l7 = aexon1[j2].getGenomeRight();
            if(getStrand() != 1 ? getStrand() == -1 && (i3 > k6 || l7 > j3) : i3 < k6 || l7 < j3)
                break;
            if(i2 == aexon.length)
                continue;
            int k8 = aexon[i2].getGenomeLeft();
            if(getStrand() != 1 ? getStrand() == -1 && j3 <= k8 : j3 >= k8)
                break;
        }

        return l2 == s2.length() - 1;
    }

    private void setAlternaLine(int i) throws Exception {
        Alternative alternative = (Alternative)_alternativeList.get(i);
        int j = alternative.getGeneIds1()[0];
        int k = alternative.getGeneIds2()[0];
        String s = alternative.getPattern1();
        String s1 = alternative.getPattern2();
        int l = alternative.getStartExonID1();
        int i1 = alternative.getStartExonID2();
        int j1 = alternative.getNumExon1();
        int k1 = alternative.getNumExon2();
        Exon aexon[] = ((Gene)_geneList.get(j)).getExons();
        Exon aexon1[] = ((Gene)_geneList.get(k)).getExons();
        if(l == -1 || l + j1 > aexon.length)
            l = 0;
        if(i1 == -1 || i1 + k1 > aexon1.length)
            i1 = 0;
        int l1 = l;
        int i2 = i1;
        int j2 = 0;
        int k2 = 0;
        ArrayList arraylist = new ArrayList();
        if(s.charAt(j2) == '1' && s1.charAt(j2) == '1') {
            if(getStrand() == 1)
                k2 = Math.max(aexon[l1].getGenomeLeft(), aexon1[i2].getGenomeLeft());
            else
                k2 = Math.min(aexon[l1].getGenomeLeft(), aexon1[i2].getGenomeLeft());
        } else
        if(s.charAt(j2) == '1' && s1.charAt(j2) == '0')
            k2 = aexon[l1].getGenomeLeft();
        else
        if(s.charAt(j2) == '0' && s1.charAt(j2) == '1')
            k2 = aexon1[i2].getGenomeLeft();
        else
            throw new Exception("p1(0)=0 and p2(0)=0");
        arraylist.add(new Integer(k2));
        for(j2 = 0; j2 < s.length() - 1; j2++) {
            if(s.charAt(j2) == '1' && s.charAt(j2 + 1) == '0' && s1.charAt(j2) == '1' && s1.charAt(j2 + 1) == '1') {
                k2 = aexon[l1].getGenomeRight();
                l1++;
            } else
            if(s.charAt(j2) == '1' && s.charAt(j2 + 1) == '1' && s1.charAt(j2) == '1' && s1.charAt(j2 + 1) == '0') {
                k2 = aexon1[i2].getGenomeRight();
                i2++;
            } else
            if(s.charAt(j2) == '1' && s.charAt(j2 + 1) == '0' && s1.charAt(j2) == '1' && s1.charAt(j2 + 1) == '0') {
                k2 = aexon[l1].getGenomeRight();
                l1++;
                i2++;
            } else
            if(s.charAt(j2) == '1' && s.charAt(j2 + 1) == '1' && s1.charAt(j2) == '0' && s1.charAt(j2 + 1) == '1')
                k2 = aexon1[i2].getGenomeLeft();
            else
            if(s.charAt(j2) == '1' && s.charAt(j2 + 1) == '0' && s1.charAt(j2) == '0' && s1.charAt(j2 + 1) == '1') {
                k2 = aexon[l1].getGenomeRight();
                l1++;
            } else
            if(s.charAt(j2) == '1' && s.charAt(j2 + 1) == '0' && s1.charAt(j2) == '0' && s1.charAt(j2 + 1) == '0') {
                k2 = aexon[l1].getGenomeRight();
                l1++;
            } else
            if(s.charAt(j2) == '0' && s.charAt(j2 + 1) == '1' && s1.charAt(j2) == '1' && s1.charAt(j2 + 1) == '1')
                k2 = aexon[l1].getGenomeLeft();
            else
            if(s.charAt(j2) == '0' && s.charAt(j2 + 1) == '1' && s1.charAt(j2) == '1' && s1.charAt(j2 + 1) == '0') {
                k2 = aexon[l1].getGenomeLeft();
                i2++;
            } else
            if(s.charAt(j2) == '0' && s.charAt(j2 + 1) == '0' && s1.charAt(j2) == '1' && s1.charAt(j2 + 1) == '0') {
                k2 = aexon1[i2].getGenomeRight();
                i2++;
            } else
            if(s.charAt(j2) == '0' && s.charAt(j2 + 1) == '1' && s1.charAt(j2) == '0' && s1.charAt(j2 + 1) == '1')
                k2 = aexon[l1].getGenomeLeft();
            else
            if(s.charAt(j2) == '0' && s.charAt(j2 + 1) == '1' && s1.charAt(j2) == '0' && s1.charAt(j2 + 1) == '0')
                k2 = aexon[l1].getGenomeLeft();
            else
            if(s.charAt(j2) == '0' && s.charAt(j2 + 1) == '0' && s1.charAt(j2) == '0' && s1.charAt(j2 + 1) == '1')
                k2 = aexon1[i2].getGenomeLeft();
            else
                throw new Exception("p1(0)=0 and p2(0)=0");
            arraylist.add(new Integer(k2));
        }

        if(s.charAt(j2) == '1' && s1.charAt(j2) == '1') {
            if(getStrand() == 1)
                k2 = Math.min(aexon[l1].getGenomeRight(), aexon1[i2].getGenomeRight());
            else
                k2 = Math.max(aexon[l1].getGenomeRight(), aexon1[i2].getGenomeRight());
        } else
        if(s.charAt(j2) == '1' && s1.charAt(j2) == '0')
            k2 = aexon[l1].getGenomeRight();
        else
        if(s.charAt(j2) == '0' && s1.charAt(j2) == '1')
            k2 = aexon1[i2].getGenomeRight();
        else
            throw new Exception("p1(0)=0 and p2(0)=0");
        arraylist.add(new Integer(k2));
        _alternaLineList.add(arraylist);
    }

    protected ArrayList getAlternativeLineList() {
        ArrayList arraylist = new ArrayList();
        for(int i = 0; i < _alternativeList.size(); i++) {
            ArrayList arraylist1 = (ArrayList)_alternaLineList.get(i);
            for(int j = 0; j < arraylist1.size(); j++)
                if(!arraylist.contains(arraylist1.get(j)))
                    arraylist.add(arraylist1.get(j));

        }

        return arraylist;
    }

    protected ArrayList getAlternalineList(int i) {
        ArrayList arraylist = new ArrayList();
        arraylist = (ArrayList)_alternaLineList.get(i);
        return arraylist;
    }

    protected ArrayList getLineInAlternaRange() {
        return _lineInAlternaRange;
    }

    protected TreeSet getLineNoAlternaRange() {
        return _lineNoAlternaRange;
    }

    private void setLineList() {
        TreeMap treemap = new TreeMap();
        for(int i = 0; i < getAlternativeLineList().size(); i++) {
            int l = ((Integer)getAlternativeLineList().get(i)).intValue();
            if(!treemap.containsKey(new Integer(l)))
                treemap.put(new Integer(l), new Integer(0));
        }

        for(int j = 0; j < getLineInAlternaRange().size(); j++) {
            int i1 = ((Integer)getLineInAlternaRange().get(j)).intValue();
            if(!treemap.containsKey(new Integer(i1)))
                treemap.put(new Integer(i1), new Integer(0));
        }

label0:
        for(int k = 0; k < getLineNoAlternaRange().size(); k++) {
            Iterator iterator = getLineNoAlternaRange().iterator();
            do {
                if(!iterator.hasNext())
                    continue label0;
                int j1 = ((Integer)iterator.next()).intValue();
                if(!treemap.containsKey(new Integer(j1)))
                    treemap.put(new Integer(j1), new Integer(0));
            } while(true);
        }

        Set set = treemap.entrySet();
        Iterator iterator1 = set.iterator();
        int k1 = 0;
        java.util.Map.Entry entry;
        for(; iterator1.hasNext(); treemap.put(entry.getKey(), new Integer(k1++)))
            entry = (java.util.Map.Entry)iterator1.next();

        _genomepos2linenoMap = treemap;
    }

    protected TreeMap getLineList() {
        return _genomepos2linenoMap;
    }
}