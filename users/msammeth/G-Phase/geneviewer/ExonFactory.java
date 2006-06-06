// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:04 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   ExonFactory.java


public class ExonFactory {

    public ExonFactory() {
    }

    public static Exon createExon(Range range, Range range1) {
        Object obj = null;
        int i = range.getLeft();
        int j = range.getRight();
        int k = range.getLeft();
        int l = range.getRight();
        if(k <= i && j <= l)
            obj = new OrfOnlyExon();
        else
        if(j < k)
            obj = new UtrOnlyExon();
        else
        if(l < i)
            obj = new UtrOnlyExon();
        else
        if(i < k)
            obj = new UtrOrfExon();
        else
        if(l < j)
            obj = new OrfUtrExon();
        if(obj != null)
            ((Exon) (obj)).setCdnaRange(range);
        return ((Exon) (obj));
    }
}