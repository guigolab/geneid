// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:04 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   AppletParameter.java


public class AppletParameter {

    private String gffUrl;
    private String seqUrl;
    private String splicingPatternUrl;
    private String seqenceName;

    public AppletParameter() {
    }

    public String getGffUrl() {
        return gffUrl;
    }

    public String getSeqenceName() {
        return seqenceName;
    }

    public String getSeqUrl() {
        return seqUrl;
    }

    public String getSplicingPatternUrl() {
        return splicingPatternUrl;
    }

    public void setGffUrl(String s) {
        gffUrl = s;
    }

    public void setSeqenceName(String s) {
        seqenceName = s;
    }

    public void setSeqUrl(String s) {
        seqUrl = s;
    }

    public void setSplicingPatternUrl(String s) {
        splicingPatternUrl = s;
    }
}