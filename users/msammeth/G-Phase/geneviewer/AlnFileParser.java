package geneviewer;

// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:03 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   AlnFileParser.java

import java.io.*;
import java.util.StringTokenizer;

class AlnFileParser {

    AlnFileParser() {
    }

    static void parseHeaderLine(String s) {
        StringTokenizer stringtokenizer = new StringTokenizer(s);
        stringtokenizer.nextToken();
        String s1 = stringtokenizer.nextToken();
        String s2 = stringtokenizer.nextToken();
        stringtokenizer.nextToken();
        stringtokenizer.nextToken();
        String s3 = stringtokenizer.nextToken();
        String s4 = stringtokenizer.nextToken();
        int i;
        int j;
        for(StringTokenizer stringtokenizer1 = new StringTokenizer(s4, ","); stringtokenizer1.hasMoreTokens(); System.out.println(s3 + "\t" + "ALN" + "\t" + "exon" + "\t" + i + "\t" + j + "\t" + "0" + "\t" + "+" + "\t" + i % 3 + "\t" + s1)) {
            String s5 = stringtokenizer1.nextToken();
            i = Integer.parseInt(s5.substring(0, s5.indexOf("..")));
            j = Integer.parseInt(s5.substring(s5.indexOf("..") + 2));
        }

    }

    public static void main(String args[]) {
        if(args.length < 1) {
            System.out.println("Usage: program Goto-file");
            System.exit(0);
        }
        try {
            FileInputStream fileinputstream = new FileInputStream(args[0]);
            BufferedReader bufferedreader = new BufferedReader(new InputStreamReader(fileinputstream));
            Object obj = null;
            do {
                String s;
                if((s = bufferedreader.readLine()) == null)
                    break;
                s = s.trim();
                if(s.startsWith("Q:"))
                    parseHeaderLine(s);
            } while(true);
        }
        catch(IOException ioexception) {
            System.err.println(ioexception);
        }
    }
}