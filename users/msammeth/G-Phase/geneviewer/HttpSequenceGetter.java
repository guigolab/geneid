// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:07 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   HttpSequenceGetter.java

import java.io.*;
import java.net.URL;

public abstract class HttpSequenceGetter
    implements SequenceGetter {

    public static final String HTTP_KEY_BASEURL = "baseUrl";
    public static final String QUERY_KEY_SPECIES = "sp";
    public static final String QUERY_KEY_CHNO = "chNo";
    public static final String QUERY_KEY_LOCIID = "lociId";
    public static final String QUERY_KEY_GENEID = "geneId";
    public static final String QUERY_KEY_CDNAID = "cdnaId";
    private static final int BUF_SIZE = 4096;
    protected URL _baseURL;

    public abstract Sequence getSequence(Object obj);

    public HttpSequenceGetter(URL url) {
        _baseURL = null;
        _baseURL = url;
    }

    protected String fetchData(URL url) throws IOException {
        char ac[] = new char[4096];
        StringBuffer stringbuffer = new StringBuffer();
        java.io.InputStream inputstream = url.openStream();
        BufferedReader bufferedreader = new BufferedReader(new InputStreamReader(inputstream));
        int i;
        while((i = bufferedreader.read(ac, 0, 4096)) != -1) 
            stringbuffer.append(ac, 0, i);
        bufferedreader.close();
        return stringbuffer.toString();
    }
}