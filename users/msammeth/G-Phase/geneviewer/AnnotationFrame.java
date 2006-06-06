package geneviewer;

// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:04 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   AnnotationFrame.java

import java.applet.AppletContext;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.PrintStream;
import java.io.UnsupportedEncodingException;
import java.net.URL;
import java.net.URLEncoder;
import java.util.HashMap;
import javax.swing.*;
import javax.swing.border.*;

public class AnnotationFrame extends JFrame {
    class BottomSpacingTitledBorder extends CompoundBorder {

        BottomSpacingTitledBorder(String s) {
            super(new BottomSpacingEmptyBorder(), new TitledBorder(s));
        }
    }

    class BottomSpacingEmptyBorder extends EmptyBorder {

        BottomSpacingEmptyBorder() {
            super(0, 0, 5, 0);
        }

        BottomSpacingEmptyBorder(int i) {
            super(0, 0, i, 0);
        }
    }

    class BottomSpacingLabel extends JLabel {

        BottomSpacingLabel() {
            super();
            setBorder(new BottomSpacingEmptyBorder());
        }

        BottomSpacingLabel(String s) {
            this();
            setText(s);
        }
    }


    public static final int FRAME_WIDTH = 780;
    public static final int FRAME_HEIGHT = 500;
    private static final int LEFT_PANEL_WIDTH = 300;
    private static final int LEFT_PANEL_HEIGHT = 450;
    private static final String GENBANK_BASE_URL = "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=nucleotide&val=";
    private static final String GENBANK_PAGE_TARGET = "_blank";
    private static final String GO_PAGE_TARGET = "_blank";
    private static final String SEQBROWSER_PAGE_TARGET = "sequence";
    private static final String ENSEMBL_URL = "http://www.ensembl.org/Homo_sapiens/textview?species=Homo_sapiens&idx=All&q=";
    private static final String CDNAID_LABEL = "UniGene ID : ";
    private static final String UGCLUSTID_LABEL = "UniGene cluster ID : ";
    private static final String CDS_LABEL = "CDS : ";
    private static final String GINUMBER_LABEL = "Gene ID : ";
    private static final String LENGTH_LABEL = "Length : ";
    private static final String GENBANK_LABEL = "Accession number : ";
    private static final String GENBANK_BTN_LABEL = "GenBank";
    private static final String ENSEMBL_BTN_LABEL = "Ensembl";
    private static final String NUCLEOTIDE_SEQ_LABEL = "nucleotide sequence ";
    private static final String NUCLEOTIDE_BTN_LABEL = "nucleotide";
    private static final String AMINO_SEQ_LABEL = "amino acid sequence ";
    private static final String AMINOACID_BTN_LABEL = "amino acid";
    private static final String TRANSLATED_SEQ_LABEL = "translation ";
    private static final String TRANSLATED_BTN_LABEL = "translation";
    private static final String PUBDBBTN_BORDER_TITLE = "Public Database";
    private static final String SEQBTN_BORDER_TITLE = "Sequence";
    private static final String CLOSE_BTN_LABEL = "Close";
    private static final String NUCLEOTIDE_SEQ_TITLE = "nucleotide sequence";
    private static final String AMINOACID_SEQ_TITLE = "amino acid sequence";
    private static final String TRANSLATED_SEQ_TITLE = "translation\n\n";
    private static final String SPLICING_BTN_LABEL = "splicing";
    private static final String BROWSER_BTN_LABEL = "browser";
    private static final int SEQ_NUM_COLUMNS = 60;
    private static final Color LEFT_PANEL_COLOR = new Color(255, 255, 224);
    private static final int COMPONENT_SPACING = 5;
    private static final String NEW_LINE = "\n";
    private static final int SEQ_TYPE_NUCLEOTIDE = 1;
    private static final int SEQ_TYPE_AMINOACID = 2;
    private static final int SEQ_TYPE_TRANSLATION_NONE = 3;
    private static final int SEQ_TYPE_TRANSLATION_EXON = 4;
    private static final int SEQ_TYPE_SPLICING_NONE = 5;
    private static final int SEQ_TYPE_SPLICING_EXON = 6;
    private static final int SEQ_TYPE_SPLICING_INTRON = 7;
    private JPanel _leftPanel;
    private JPanel _southPanel;
    private JPanel _centerPanel;
    private JLabel _cdnaIdLabel;
    private JLabel _uniGeneClusterIdLabel;
    private JLabel _cdsRegionLabel;
    private JLabel _giNumberLabel;
    private JLabel _lengthLabel;
    private JLabel _genBankAccLabel;
    private JLabel _geneNameLabel;
    private JPanel _pubdbBtnPanel;
    private JButton _genBankButton;
    private JButton _ensemblButton;
    private JPanel _sequenceBtnPanel;
    private JButton _nucSequenceButton;
    private JButton _amiSequenceButton;
    private JButton _translationButton;
    private JButton _splicingButton;
    private JButton _browserButton;
    private JScrollPane _sequenceScrollPane;
    private SequenceTextPane _sequenceTextPane;
    private JButton _closeButton;
    private JPanel _contentPane;
    private GenePane _genePane;
    private Gene _noticedGene;
    private GeneShape _highlightGeneShape;
    private Exon _alternativeExon;
    private boolean _bHighlightAlternativeExon;
    private Sequence _genomeSequence;
    private Sequence _splicingSequence;
    private int _currentSequenceType;
    SequenceHighlighter _sequenceHighlighter;

    public AnnotationFrame(GenePane genepane) {
        _leftPanel = null;
        _southPanel = null;
        _centerPanel = null;
        _cdnaIdLabel = null;
        _uniGeneClusterIdLabel = null;
        _cdsRegionLabel = null;
        _giNumberLabel = null;
        _lengthLabel = null;
        _genBankAccLabel = null;
        _geneNameLabel = null;
        _pubdbBtnPanel = null;
        _genBankButton = null;
        _ensemblButton = null;
        _sequenceBtnPanel = null;
        _nucSequenceButton = null;
        _amiSequenceButton = null;
        _translationButton = null;
        _splicingButton = null;
        _browserButton = null;
        _sequenceScrollPane = null;
        _sequenceTextPane = null;
        _closeButton = null;
        _contentPane = null;
        _genePane = null;
        _noticedGene = null;
        _highlightGeneShape = null;
        _alternativeExon = null;
        _bHighlightAlternativeExon = false;
        _genomeSequence = null;
        _splicingSequence = null;
        _currentSequenceType = 1;
        _sequenceHighlighter = null;
        _genePane = genepane;
        _contentPane = (JPanel)getContentPane();
        createComponents();
        layoutComponents();
        setActionListeners();
    }

    private void createComponents() {
        _leftPanel = new JPanel();
        _leftPanel.setBackground(LEFT_PANEL_COLOR);
        _leftPanel.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
        _leftPanel.setMaximumSize(new Dimension(300, 450));
        _centerPanel = new JPanel();
        _southPanel = new JPanel();
        _cdnaIdLabel = new BottomSpacingLabel();
        _uniGeneClusterIdLabel = new BottomSpacingLabel();
        _cdsRegionLabel = new BottomSpacingLabel();
        _giNumberLabel = new BottomSpacingLabel();
        _lengthLabel = new BottomSpacingLabel();
        _genBankAccLabel = new BottomSpacingLabel();
        _geneNameLabel = new BottomSpacingLabel();
        _pubdbBtnPanel = new JPanel();
        _pubdbBtnPanel.setBackground(LEFT_PANEL_COLOR);
        _pubdbBtnPanel.setBorder(new BottomSpacingTitledBorder("Public Database"));
        _genBankButton = new JButton("GenBank");
        _ensemblButton = new JButton("Ensembl");
        _sequenceBtnPanel = new JPanel();
        _sequenceBtnPanel.setBackground(LEFT_PANEL_COLOR);
        _sequenceBtnPanel.setBorder(new BottomSpacingTitledBorder("Sequence"));
        _nucSequenceButton = new JButton("nucleotide");
        _amiSequenceButton = new JButton("amino acid");
        _translationButton = new JButton("translation");
        _splicingButton = new JButton("splicing");
        _browserButton = new JButton("browser");
        _sequenceTextPane = new SequenceTextPane();
        _sequenceScrollPane = new JScrollPane(_sequenceTextPane);
        _closeButton = new JButton("Close");
    }

    private void layoutLeftComponents() {
        _leftPanel.setLayout(new BoxLayout(_leftPanel, 1));
        _pubdbBtnPanel.setLayout(new BoxLayout(_pubdbBtnPanel, 0));
        _pubdbBtnPanel.setAlignmentX(0.0F);
        _pubdbBtnPanel.add(_genBankButton);
        _pubdbBtnPanel.add(_ensemblButton);
        _sequenceBtnPanel.setLayout(new GridLayout(2, 3));
        _sequenceBtnPanel.setAlignmentX(0.0F);
        _sequenceBtnPanel.setMaximumSize(new Dimension(400, 88));
        _sequenceBtnPanel.add(_nucSequenceButton);
        _sequenceBtnPanel.add(_amiSequenceButton);
        _sequenceBtnPanel.add(_translationButton);
        _sequenceBtnPanel.add(_splicingButton);
        _sequenceBtnPanel.add(_browserButton);
        _leftPanel.add(_cdnaIdLabel);
        _leftPanel.add(_uniGeneClusterIdLabel);
        _leftPanel.add(_cdsRegionLabel);
        _leftPanel.add(_giNumberLabel);
        _leftPanel.add(_lengthLabel);
        _leftPanel.add(_genBankAccLabel);
        _leftPanel.add(_pubdbBtnPanel);
        _leftPanel.add(_sequenceBtnPanel);
    }

    private void layoutCenterComponents() {
        _centerPanel.add(_sequenceScrollPane, "Center");
    }

    private void layoutSouthComponents() {
        _southPanel.add(_closeButton);
    }

    private void layoutComponents() {
        layoutLeftComponents();
        layoutCenterComponents();
        layoutSouthComponents();
        _contentPane.add(_southPanel, "South");
        _contentPane.add(_leftPanel, "West");
        _contentPane.add(_sequenceScrollPane, "Center");
    }

    private void setActionListeners() {
        _genBankButton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent actionevent) {
                try {
                    AppletContext appletcontext = _genePane.getAppletContext();
                    appletcontext.showDocument(new URL(getGenBankPageUrl()), "_blank");
                }
                catch(Exception exception) {
                    System.err.println(exception);
                }
            }

             {
                super();
            }
        });
        _ensemblButton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent actionevent) {
                try {
                    AppletContext appletcontext = _genePane.getAppletContext();
                    appletcontext.showDocument(new URL(getEnsemblPageUrl()), "_blank");
                }
                catch(Exception exception) {
                    System.err.println(exception);
                }
            }

             {
                super();
            }
        });
        _nucSequenceButton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent actionevent) {
                showNucleotideSequence();
            }

             {
                super();
            }
        });
        _amiSequenceButton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent actionevent) {
                showAminoAcidSequence();
            }

             {
                super();
            }
        });
        _translationButton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent actionevent) {
                showTranslatedSequence();
            }

             {
                super();
            }
        });
        _splicingButton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent actionevent) {
                showSplicingSequence();
            }

             {
                super();
            }
        });
        _browserButton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent actionevent) {
                try {
                    AppletContext appletcontext = _genePane.getAppletContext();
                    appletcontext.showDocument(new URL(_genePane.getBaseURL(), getHtmlSeqPath()), "sequence");
                }
                catch(Exception exception) {
                    System.err.println(exception);
                }
            }

             {
                super();
            }
        });
        _closeButton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent actionevent) {
                setVisible(false);
            }

             {
                super();
            }
        });
    }

    private void addGridBagComponent(JPanel jpanel, Component component, int i, int j, int k, int l) {
        GridBagConstraints gridbagconstraints = new GridBagConstraints();
        gridbagconstraints.fill = 0;
        gridbagconstraints.gridx = i;
        gridbagconstraints.gridy = j;
        gridbagconstraints.gridwidth = k;
        gridbagconstraints.gridheight = l;
        gridbagconstraints.anchor = 17;
        gridbagconstraints.insets = new Insets(1, 1, 1, 1);
        jpanel.add(component, gridbagconstraints);
    }

    public void init() {
        _highlightGeneShape = null;
        _alternativeExon = null;
        _bHighlightAlternativeExon = false;
        _sequenceHighlighter = NullSequenceHighlighter.getInstance();
    }

    public void setNoticedGene(Gene gene) {
        _noticedGene = gene;
        HashMap hashmap = getUrlParametersMap();
        Object obj = new NormalSequenceGetter(_genePane.getBaseURL());
        _genomeSequence = ((SequenceGetter) (obj)).getSequence(hashmap);
        obj = new SplicingSequenceGetter(_genePane.getBaseURL());
        _splicingSequence = ((SequenceGetter) (obj)).getSequence(hashmap);
        int i = gene.getLength();
        if(i == 0)
            i = _genomeSequence.getNucleotideSequence().length();
        _cdnaIdLabel.setText("UniGene ID : " + _noticedGene.getCdnaId());
        _uniGeneClusterIdLabel.setText("UniGene cluster ID : " + getShowingString(_noticedGene.getUniGeneId()));
        _cdsRegionLabel.setText("CDS : " + String.valueOf(_noticedGene.getCdsLeft()) + " .. " + String.valueOf(_noticedGene.getCdsRight()));
        _giNumberLabel.setText("Gene ID : " + getShowingString(_noticedGene.getGiNumber()));
        _lengthLabel.setText("Length : " + getShowingString(String.valueOf(i)));
        _genBankAccLabel.setText("Accession number : " + _noticedGene.getGenBankAccession());
        _geneNameLabel.setText(_noticedGene.getGeneName());
        setTitle(_noticedGene.getGeneName());
        showNucleotideSequence();
    }

    private String getShowingString(String s) {
        if(s.length() == 0 || s.equals("0"))
            return "(Unknown)";
        else
            return s;
    }

    public void showNucleotideSequence() {
        String s = ">" + _noticedGene.getCdnaId() + "\n" + "No sequence data.";
        if(_genomeSequence != null && _genomeSequence.hasSequence())
            s = _genomeSequence.getNucleotideSequence();
        StringBuffer stringbuffer = new StringBuffer();
        stringbuffer.append("nucleotide sequence");
        stringbuffer.append("\n");
        stringbuffer.append("\n");
        stringbuffer.append(_genomeSequence.getAnnotation());
        stringbuffer.append("\n");
        stringbuffer.append(StringUtils.getLinedString(s, 60));
        _sequenceTextPane.setText(stringbuffer.toString());
        _sequenceTextPane.setCaretPosition(0);
        _currentSequenceType = 1;
    }

    public void showAminoAcidSequence() {
        String s = "No sequence data.";
        if(_genomeSequence != null && _genomeSequence.hasSequence())
            s = _genomeSequence.getAminoAcidSequence(_noticedGene.getCdsLeft(), _noticedGene.getCdsRight());
        StringBuffer stringbuffer = new StringBuffer();
        stringbuffer.append("amino acid sequence");
        stringbuffer.append("\n");
        stringbuffer.append("\n");
        stringbuffer.append(">");
        stringbuffer.append(_noticedGene.getCdnaId());
        stringbuffer.append("\n");
        stringbuffer.append(StringUtils.getLinedString(s, 60));
        stringbuffer.append("\n");
        _sequenceTextPane.setText(stringbuffer.toString());
        _sequenceTextPane.setCaretPosition(0);
        _currentSequenceType = 2;
    }

    public void showTranslatedSequence() {
        String s = "No sequence data.";
        if(_genomeSequence != null && _genomeSequence.hasSequence())
            s = _genomeSequence.getTranslatedSequence(_noticedGene.getCdsLeft(), _noticedGene.getCdsRight(), 60);
        _sequenceTextPane.setText("translation\n\n" + s);
        if(_highlightGeneShape != null) {
            _currentSequenceType = 4;
            setSequenceHighlighter();
            _sequenceHighlighter.highlight();
        } else {
            _currentSequenceType = 3;
            _sequenceTextPane.setCaretPosition(0);
        }
    }

    public void showSplicingSequence() {
        setSequenceHighlighter();
        _sequenceTextPane.setText(_splicingSequence.getNucleotideSequence());
        if(_highlightGeneShape == null) {
            _currentSequenceType = 5;
            _sequenceTextPane.setCaretPosition(0);
        } else {
            if(_highlightGeneShape.getSegmentType() == 1)
                _currentSequenceType = 6;
            else
            if(_highlightGeneShape.getSegmentType() == 2)
                _currentSequenceType = 7;
            setSequenceHighlighter();
            _sequenceHighlighter.highlight();
        }
    }

    private void setSequenceHighlighter() {
        if(_highlightGeneShape != null) {
            if(_highlightGeneShape.getSegmentType() == 1) {
                if(_currentSequenceType == 6)
                    _sequenceHighlighter = new ExonHighlighter(_genePane.getLocus(), _highlightGeneShape, _sequenceTextPane, _sequenceScrollPane, _alternativeExon, _bHighlightAlternativeExon);
                else
                if(_currentSequenceType == 4)
                    _sequenceHighlighter = new TranslatedSequenceHighlighter(_genePane.getLocus(), _highlightGeneShape, _sequenceTextPane, _sequenceScrollPane, _alternativeExon, _bHighlightAlternativeExon);
                else
                    _sequenceHighlighter = NullSequenceHighlighter.getInstance();
            } else
            if(_highlightGeneShape.getSegmentType() == 2)
                _sequenceHighlighter = new IntronHighlighter(_genePane.getLocus(), _highlightGeneShape, _sequenceTextPane, _sequenceScrollPane);
            else
                _sequenceHighlighter = NullSequenceHighlighter.getInstance();
        } else {
            _sequenceHighlighter = NullSequenceHighlighter.getInstance();
        }
    }

    private int getSequenceTextPaneScrollPos() {
        return _sequenceScrollPane.getVerticalScrollBar().getValue();
    }

    private HashMap getUrlParametersMap() {
        Locus locus = _genePane.getLocus();
        HashMap hashmap = new HashMap();
        hashmap.put("sp", _genePane.getSpecies());
        hashmap.put("chNo", locus.getChNo());
        hashmap.put("lociId", String.valueOf(locus.getLociId()));
        hashmap.put("geneId", String.valueOf(_noticedGene.getId()));
        hashmap.put("cdnaId", _noticedGene.getCdnaId());
        return hashmap;
    }

    private String getHtmlSeqPath() throws UnsupportedEncodingException {
        Locus locus = _genePane.getLocus();
        StringBuffer stringbuffer = new StringBuffer("htmlseq.php?");
        stringbuffer.append("s=");
        stringbuffer.append(_genePane.getSpecies());
        stringbuffer.append("&c=");
        stringbuffer.append(locus.getChNo());
        stringbuffer.append("&l=");
        stringbuffer.append(locus.getLociId());
        stringbuffer.append("&g=");
        stringbuffer.append(_noticedGene.getId());
        stringbuffer.append("&d=");
        stringbuffer.append(URLEncoder.encode(_noticedGene.getCdnaId(), "UTF-8"));
        stringbuffer.append("&r=");
        stringbuffer.append(_noticedGene.getCdsLeft());
        stringbuffer.append(',');
        stringbuffer.append(_noticedGene.getCdsRight());
        stringbuffer.append("&t=");
        stringbuffer.append(_currentSequenceType);
        if(_highlightGeneShape != null) {
            Range arange[] = _sequenceHighlighter.getHighlightRanges();
            if(arange != null) {
                stringbuffer.append("&h=");
                for(int i = 0; i < arange.length; i++) {
                    if(i != 0)
                        stringbuffer.append(',');
                    stringbuffer.append(arange[i].getLeft());
                    stringbuffer.append(',');
                    stringbuffer.append(arange[i].getRight());
                }

            }
        }
        if(_highlightGeneShape != null && _highlightGeneShape.getExon() != null) {
            stringbuffer.append("&e=");
            stringbuffer.append(_highlightGeneShape.getExon().getGenomeLeft());
            stringbuffer.append(',');
            stringbuffer.append(_highlightGeneShape.getExon().getGenomeRight());
            stringbuffer.append("&n=");
            stringbuffer.append(_highlightGeneShape.getExon().getCdnaRange().getLeft());
            stringbuffer.append(',');
            stringbuffer.append(_highlightGeneShape.getExon().getCdnaRange().getRight());
        }
        stringbuffer.append("#h");
        return stringbuffer.toString();
    }

    public void setAlternativeExon(Exon exon) {
        _alternativeExon = exon;
    }

    public void setHighlightGeneShape(GeneShape geneshape) {
        _highlightGeneShape = geneshape;
    }

    public void setHighlightAlternativeExon(boolean flag) {
        _bHighlightAlternativeExon = flag;
    }

    private String getGenBankPageUrl() {
        return "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=nucleotide&val=" + _noticedGene.getGiNumber();
    }

    private String getEnsemblPageUrl() {
        StringBuffer stringbuffer = new StringBuffer();
        stringbuffer.append("http://www.ensembl.org/Homo_sapiens/textview?species=Homo_sapiens&idx=All&q=");
        stringbuffer.append(_noticedGene.getGenBankAccession());
        return stringbuffer.toString();
    }





}