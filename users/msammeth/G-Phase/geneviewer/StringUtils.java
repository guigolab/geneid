// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:07 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   StringUtils.java


public class StringUtils {

    public static final String NEW_LINE = "\n";
    public static final String SPACE = " ";

    public StringUtils() {
    }

    public static String getLinedString(String s, int i) {
        StringBuffer stringbuffer = new StringBuffer();
        int j = s.length();
        for(int k = 0; k < j; k += i) {
            if(k + i <= j)
                stringbuffer.append(s.substring(k, k + i));
            else
                stringbuffer.append(s.substring(k));
            stringbuffer.append("\n");
        }

        return stringbuffer.toString();
    }

    public static boolean isLowerCase(String s) {
        for(int i = 0; i < s.length(); i++)
            if(!Character.isLowerCase(s.charAt(i)))
                return false;

        return true;
    }

    public static boolean isUpperCase(String s) {
        for(int i = 0; i < s.length(); i++)
            if(s.charAt(i) != '$' && !Character.isUpperCase(s.charAt(i)))
                return false;

        return true;
    }
}