// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:07 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   LocusSaxHandler.java

import java.util.LinkedList;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

public class LocusSaxHandler extends DefaultHandler {

    private static final String LOCUS_TAG = "locus";
    private static final String GENES_TAG = "genes";
    private static final String GENE_TAG = "gene";
    private static final String EXON_TAG = "exon";
    private static final String ALTERNATIVES_TAG = "alternatives";
    private static final String ALTERNATIVE_TAG = "alternative";
    private static final String HOMOLOGENES_TAG = "homologenes";
    private static final String HOMOLOGENE_TAG = "homologene";
    private static final String ID_ATTR = "id";
    private static final String CH_NO_ATTR = "ch_no";
    private static final String STRAND_ATTR = "strand";
    private static final String CDNA_ID_ATTR = "cdna_id";
    private static final String GENE_NAME_ATTR = "gene_name";
    private static final String CDS_LEFT_ATTR = "cds_left";
    private static final String CDS_RIGHT_ATTR = "cds_right";
    private static final String GENBANK_ACC_ATTR = "gb_acc";
    private static final String GI_NO_ATTR = "gi_no";
    private static final String UNIGENE_ID_ATTR = "unigene_id";
    private static final String LENGTH_ATTR = "length";
    private static final String NMD_ATTR = "nmd";
    private static final String CDNA_LEFT_ATTR = "cdna_left";
    private static final String CDNA_RIGHT_ATTR = "cdna_right";
    private static final String GENOME_LEFT_ATTR = "genome_left";
    private static final String GENOME_RIGHT_ATTR = "genome_right";
    private static final String UTR_LEFT_ATTR = "utr_left";
    private static final String UTR_RIGHT_ATTR = "utr_right";
    private static final String ORF_LEFT_ATTR = "orf_left";
    private static final String ORF_RIGHT_ATTR = "orf_right";
    private static final String TYPE_ATTR = "type";
    private static final String NAGNAG_ATTR = "nagnag";
    private static final String START_EXON_ID_ATTR = "start_exon_id";
    private static final String NUM_EXON_ATTR = "num_exon";
    private static final String PATTERN_ATTR = "pattern";
    private static final String SPECIES_ATTR = "species";
    private static final String LOCIID_ATTR = "loci_id";
    private Locus _locus;
    private Gene _gene;
    private Alternative _alternative;
    private Exon _exon;
    private LinkedList _elemNameStack;
    private int _minGenomePos;
    private int _maxGenomePos;

    public LocusSaxHandler() {
        _locus = new Locus();
        _gene = null;
        _alternative = null;
        _exon = null;
        _elemNameStack = new LinkedList();
        _minGenomePos = 0x7fffffff;
        _maxGenomePos = 0x80000000;
    }

    public void startElement(String s, String s1, String s2, Attributes attributes) throws SAXException {
        if(s2.equals("locus")) {
            _locus.setLociId(intValue(attributes.getValue("id")));
            _locus.setChNo(attributes.getValue("ch_no"));
            _locus.setStrand(attributes.getValue("strand"));
        } else
        if(!s2.equals("genes"))
            if(s2.equals("gene")) {
                if(getParentElemName().equals("genes")) {
                    _gene = new Gene(_locus);
                    _gene.setId(intValue(attributes.getValue("id")));
                    _gene.setCdnaId(attributes.getValue("cdna_id"));
                    _gene.setGeneName(attributes.getValue("gene_name"));
                    _gene.setCdsLeft(intValue(attributes.getValue("cds_left")));
                    _gene.setCdsRight(intValue(attributes.getValue("cds_right")));
                    _gene.setGenBankAccession(attributes.getValue("gb_acc"));
                    _gene.setGiNumber(attributes.getValue("gi_no"));
                    _gene.setUniGeneId(attributes.getValue("unigene_id"));
                    _gene.setLength(intValue(attributes.getValue("length")));
                    _gene.setNmd(Boolean.valueOf(attributes.getValue("nmd")).booleanValue());
                } else
                if(getParentElemName().equals("alternative")) {
                    if(_alternative.getGeneIds1() == null)
                        _alternative.setGene1(attributes.getValue("id"));
                    else
                        _alternative.setGene2(attributes.getValue("id"));
                    if(_alternative.getPattern1() == null)
                        _alternative.setPattern1(attributes.getValue("pattern"));
                    else
                        _alternative.setPattern2(attributes.getValue("pattern"));
                    if(_alternative.getStartExonID1() == -1)
                        _alternative.setStartExonID1(intValue(attributes.getValue("start_exon_id")));
                    else
                        _alternative.setStartExonID2(intValue(attributes.getValue("start_exon_id")));
                    if(_alternative.getNumExon1() == -1)
                        _alternative.setNumExon1(intValue(attributes.getValue("num_exon")));
                    else
                        _alternative.setNumExon2(intValue(attributes.getValue("num_exon")));
                }
            } else
            if(s2.equals("exon")) {
                if(getParentElemName().equals("gene")) {
                    _exon = createExon(_gene, attributes);
                    int i = _exon.getGenomeLeft();
                    int j = _exon.getGenomeRight();
                    if(i < j) {
                        if(_minGenomePos > i)
                            _minGenomePos = i;
                        if(_maxGenomePos < j)
                            _maxGenomePos = j;
                    } else {
                        if(_minGenomePos > j)
                            _minGenomePos = j;
                        if(_maxGenomePos < i)
                            _maxGenomePos = i;
                    }
                } else
                if(getParentElemName().equals("alternative")) {
                    _exon = new Exon();
                    _exon.setGenomeLeft(intValue(attributes.getValue("genome_left")));
                    _exon.setGenomeRight(intValue(attributes.getValue("genome_right")));
                }
            } else
            if(!s2.equals("alternatives"))
                if(s2.equals("alternative")) {
                    _alternative = new Alternative();
                    _alternative.setId(intValue(attributes.getValue("id")));
                    _alternative.setType(attributes.getValue("type"));
                    _alternative.setGenomeLeft(intValue(attributes.getValue("genome_left")));
                    _alternative.setGenomeRight(intValue(attributes.getValue("genome_right")));
                    _alternative.setNagnag(Boolean.valueOf(attributes.getValue("nagnag")).booleanValue());
                } else
                if(!s2.equals("homologenes") && s2.equals("homologene")) {
                    String s3 = attributes.getValue("species");
                    int k = intValue(attributes.getValue("loci_id"));
                    _locus.addHomoloGene(s3, k);
                }
        _elemNameStack.addLast(s2);
    }

    public void endElement(String s, String s1, String s2) throws SAXException {
        _elemNameStack.removeLast();
        if(s2.equals("locus")) {
            _locus.setMinGenomePos(_minGenomePos);
            _locus.setMaxGenomePos(_maxGenomePos);
        } else
        if(s2.equals("gene")) {
            if(getParentElemName().equals("genes")) {
                _locus.addGene(_gene);
                _gene = null;
            }
        } else
        if(s2.equals("alternative")) {
            _locus.addSplicingPattern(_alternative);
            _alternative = null;
        } else
        if(s2.equals("exon"))
            if(getParentElemName().equals("gene")) {
                _gene.addExon(_exon);
                _exon = null;
            } else
            if(getParentElemName().equals("alternative")) {
                _alternative.addExon(_exon);
                _exon = null;
            }
    }

    public Locus getLocus() {
        return _locus;
    }

    private String getParentElemName() {
        return (String)_elemNameStack.getLast();
    }

    private int intValue(String s) {
        int i = 0x80000000;
        if(s != null)
            i = Integer.parseInt(s);
        return i;
    }

    private Exon createExon(Gene gene, Attributes attributes) {
        Exon exon = null;
        int i = intValue(attributes.getValue("cdna_left"));
        int j = intValue(attributes.getValue("cdna_right"));
        int k = intValue(attributes.getValue("genome_left"));
        int l = intValue(attributes.getValue("genome_right"));
        String s = attributes.getValue("utr_left");
        String s1 = attributes.getValue("utr_right");
        String s2 = attributes.getValue("orf_left");
        String s3 = attributes.getValue("orf_right");
        exon = ExonFactory.createExon(new Range(i, j), gene.getCdsRange());
        exon.setGenomeLeft(k);
        exon.setGenomeRight(l);
        if(s != null && s1 != null)
            exon.setUtrRange(intValue(s), intValue(s1));
        if(s2 != null && s3 != null)
            exon.setOrfRange(intValue(s2), intValue(s3));
        return exon;
    }
}