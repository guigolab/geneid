// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:07 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   SequenceUtils.java


public class SequenceUtils {

    public SequenceUtils() {
    }

    static char translateCodon(String s, int i, boolean flag) {
        if(s == null)
            return ' ';
        String s1 = s.substring(i, i + 3).toUpperCase();
        if(flag) {
            char ac[] = new char[3];
            ac[0] = s1.charAt(2);
            ac[1] = s1.charAt(1);
            ac[2] = s1.charAt(0);
            s1 = new String(ac);
        }
        if(s1.equals("ATT") || s1.equals("ATC") || s1.equals("ATA"))
            return 'I';
        if(s1.startsWith("CT"))
            return 'L';
        if(s1.equals("TTA") || s1.equals("TTG"))
            return 'L';
        if(s1.startsWith("GT"))
            return 'V';
        if(s1.equals("TTT") || s1.equals("TTC"))
            return 'F';
        if(s1.equals("ATG"))
            return 'M';
        if(s1.equals("TGT") || s1.equals("TGC"))
            return 'C';
        if(s1.startsWith("GC"))
            return 'A';
        if(s1.startsWith("GG"))
            return 'G';
        if(s1.startsWith("CC"))
            return 'P';
        if(s1.startsWith("AC"))
            return 'T';
        if(s1.startsWith("TC"))
            return 'S';
        if(s1.equals("AGT") || s1.equals("AGC"))
            return 'S';
        if(s1.equals("TAT") || s1.equals("TAC"))
            return 'Y';
        if(s1.equals("TGG"))
            return 'W';
        if(s1.equals("CAA") || s1.equals("CAG"))
            return 'Q';
        if(s1.equals("AAT") || s1.equals("AAC"))
            return 'N';
        if(s1.equals("CAT") || s1.equals("CAC"))
            return 'H';
        if(s1.equals("GAA") || s1.equals("GAG"))
            return 'E';
        if(s1.equals("GAT") || s1.equals("GAC"))
            return 'D';
        if(s1.equals("AAA") || s1.equals("AAG"))
            return 'K';
        if(s1.startsWith("CG"))
            return 'R';
        return !s1.equals("AGA") && !s1.equals("AGG") ? '$' : 'R';
    }

    static String getAminoAcidSequence(String s) {
        StringBuffer stringbuffer = new StringBuffer();
        for(int i = 0; i + 3 <= s.length(); i += 3)
            stringbuffer.append(translateCodon(s, i, false));

        return stringbuffer.toString();
    }

    static String makeComplement(String s) {
        if(s == null)
            return null;
        StringBuffer stringbuffer = new StringBuffer();
        for(int i = 0; i < s.length(); i++)
            switch(s.charAt(i)) {
            case 97: // 'a'
                stringbuffer.append('t');
                break;

            case 65: // 'A'
                stringbuffer.append('T');
                break;

            case 99: // 'c'
                stringbuffer.append('g');
                break;

            case 67: // 'C'
                stringbuffer.append('G');
                break;

            case 103: // 'g'
                stringbuffer.append('c');
                break;

            case 71: // 'G'
                stringbuffer.append('C');
                break;

            case 116: // 't'
                stringbuffer.append('a');
                break;

            case 84: // 'T'
                stringbuffer.append('A');
                break;
            }

        return stringbuffer.toString();
    }

    static String makeReverse(String s) {
        if(s == null)
            return null;
        StringBuffer stringbuffer = new StringBuffer();
        int i = s.length();
        for(int j = 0; j < i; j++)
            stringbuffer.append(s.charAt(i - j - 1));

        return stringbuffer.toString();
    }

    public static int complementPos2directPos(int i, int j, int k) {
        return (j - k) + i;
    }

    public static int directPos2complementPos(int i, int j, int k) {
        return complementPos2directPos(i, j, k);
    }
}