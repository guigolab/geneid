// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:07 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   NormalSequenceGetter.java

import java.io.*;
import java.net.*;
import java.util.HashMap;

public class NormalSequenceGetter extends HttpSequenceGetter {

    private static final String FETCH_FILE_NAME = "fetchseq.php";

    public NormalSequenceGetter(URL url) {
        super(url);
    }

    public Sequence getSequence(Object obj) {
        HashMap hashmap = (HashMap)obj;
        String s = (String)hashmap.get("sp");
        String s1 = (String)hashmap.get("cdnaId");
        String s2 = "";
        String s3 = "";
        try {
            s3 = fetchData(new URL(_baseURL, getQueryUrl(s, s1)));
            if(s3.length() > 0 && s3.charAt(0) == '>') {
                int i = s3.indexOf("\n");
                s2 = s3.substring(0, i);
                s3 = s3.substring(i + 1);
            }
        }
        catch(MalformedURLException malformedurlexception) {
            System.err.println(malformedurlexception);
        }
        catch(IOException ioexception) {
            System.err.println(ioexception);
        }
        catch(Exception exception) {
            System.err.println(exception);
        }
        return new Sequence(s2, s3);
    }

    private String getQueryUrl(String s, String s1) throws UnsupportedEncodingException {
        return "fetchseq.php?cdnaId=" + URLEncoder.encode(s1, "UTF-8") + "&sp=" + s;
    }
}