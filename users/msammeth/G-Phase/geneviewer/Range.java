// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:07 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   Range.java


public class Range {

    private int _left;
    private int _right;

    Range() {
        _left = 0;
        _right = 0;
    }

    Range(int i, int j) {
        _left = 0;
        _right = 0;
        _left = i;
        _right = j;
    }

    public boolean contains(Range range) {
        boolean flag = false;
        int i;
        int j;
        if((_left - _right) * (range.getLeft() - range.getRight()) < 0) {
            i = range.getRight();
            j = range.getLeft();
        } else {
            i = range.getLeft();
            j = range.getRight();
        }
        if(_left - _right < 0) {
            if(_left <= i && j <= _right)
                flag = true;
        } else
        if(_left >= i && j >= _right)
            flag = true;
        return flag;
    }

    public int length() {
        return Math.abs(_left - _right) + 1;
    }

    public void setRange(int i, int j) {
        setLeft(i);
        setRight(j);
    }

    public int getLeft() {
        return _left;
    }

    public void setLeft(int i) {
        _left = i;
    }

    public int getRight() {
        return _right;
    }

    public void setRight(int i) {
        _right = i;
    }
}