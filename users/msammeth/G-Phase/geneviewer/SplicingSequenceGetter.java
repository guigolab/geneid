// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:07 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   SplicingSequenceGetter.java

import java.io.IOException;
import java.io.PrintStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.HashMap;

public class SplicingSequenceGetter extends HttpSequenceGetter {

    private static final String FETCH_FILE_NAME = "fetchdia.php";

    public SplicingSequenceGetter(URL url) {
        super(url);
    }

    public Sequence getSequence(Object obj) {
        String s = "";
        try {
            URL url = new URL(_baseURL, getQueryUrl((HashMap)obj));
            s = fetchData(url);
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
        return new Sequence(s);
    }

    private String getQueryUrl(HashMap hashmap) {
        StringBuffer stringbuffer = new StringBuffer("fetchdia.php");
        stringbuffer.append('?');
        stringbuffer.append("sp");
        stringbuffer.append('=');
        stringbuffer.append((String)hashmap.get("sp"));
        stringbuffer.append('&');
        stringbuffer.append("chNo");
        stringbuffer.append('=');
        stringbuffer.append((String)hashmap.get("chNo"));
        stringbuffer.append('&');
        stringbuffer.append("lociId");
        stringbuffer.append('=');
        stringbuffer.append((String)hashmap.get("lociId"));
        stringbuffer.append('&');
        stringbuffer.append("geneId");
        stringbuffer.append('=');
        stringbuffer.append((String)hashmap.get("geneId"));
        return stringbuffer.toString();
    }
}