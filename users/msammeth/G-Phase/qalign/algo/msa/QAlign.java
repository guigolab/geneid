package qalign.algo.msa;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;

import qalign.algo.CancelException;
import qalign.algo.CostTable;
import qalign.gui.GUIProxy;

/**
 * @author micha
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class QAlign extends MSA {



        /**
         * First edge of optimal path through hyperspace.
         */
        protected Edge rootEdge= null;

        /**
         * A flag indicating the need to calculate the epsilon values.<br>
         * <br>
         * As a  default (<code>eflag= 1</code>),  the  program  calculates
         * a progressive  multiple  alignment  and  uses it to set epsilons
         * for each pairwise  alignment.   Frequently the  "optimal multiple
         * alignment" (simultaneous alignment) will be found to have observed
         * epsilons exceeding those supplied  or calculated.<br>
         * When this is the case, it is advisable to  rerun  the  program
         * using  suitably  augmented epsilons. <br>
         * <br>
         * Alternatively user  specified  epsilons for each pairwise alignment
         * may be used (<code>eflag= 0</code>).<br>
         * @see #readEpsilons()
         */
        protected byte eflag= 1;


        /**
         * Switches between progressive and simultaneous multiple
         * alignment.<br>
         * The progressive alignment and other data always is produced
         * by the program before the it attempts to produce a simultaneous
         * multiple alignment.<br>
         * The default mode (<code>mflag= 1</code>) turns on the simultaneous
         * alignment.
         */
        protected byte mflag= 1;


        /**
         * Indicates, whether a progressive alignment layout has been found.
         */
        protected boolean foundProgressive= false;

        /**
         * Indicates, whether a simultaneous alignment layout has been found.
         */
        protected boolean foundSimultaneous= false;

        /**
         * Time the <code>run()</code> method needed.
         */
        protected long deltatime= 0l;

        /**
         * Reports whether a simultanous multiple alignment has been found.
         */
        public boolean isFoundSimultaneous() {

                return foundSimultaneous;
        }

        /**
         * Reports whether a progressive multiple alignment has been found.
         */
        public boolean isFoundProgressive() {

                return foundProgressive;
        }

        public static void main(String[] args) throws CancelException {

                /*
                                > APK   CAPK  BOVINE CARDIAC MUSCLE CYCLIC AMP-DEPENDENT (ALPHA)
                                > LCK       MLCK  RAT SKELETAL MUSCLE MYOSIN LIGHT CHAIN KINASE
                                > SKH   SKH  HELA MYSTERY PUTATIVE PROTEIN KINASE
                                > D28     CD28  S. CEREVISIAE CELL CYCLE CONTROL PROTEIN KINASE
                                > EE1   WEE1  S. POMBE MITOTIC INHIBITOR
                                > AF1   RAF1  HUMAN C-RAF-1 ONCOGENE
                                > MOS   CMOS  HUMAN C-MOS ONCOGENE
                                > SRC   CSRC  CHICKEN CELLULAR ONCOGENE c-src
                                > FES  THIS IS VFES TYROSINE KINASEzwz
                                > DGM  PDGF RECEPTOR, MOUSE KINASE REGION
                                > GFR   EGFR  HUMAN EPIDERMAL GROWTH FACTOR RECEPTOR
                                > SVK   HSVK  HERPES SIMPLEX VIRUS PUTATIVE PROTEIN KINASE
                */
/*String[] seqs =
                        {
                                "DQFERIKTLGTGSFGRVMLVKHMETGNHYAMKILDKQKVVKLKQIEHTLNEKRILQAVNFPFLVKLEFSFKDNSNLYMVMEYVPGGEMFSHLRRIGRFSEPHARFYAAQIVLTFEYLHSLDLIYRDLKPENLLIDQQGYIQVTDFGFAKRVKGRTWTLCGTPEYLAPEIILSKGYNKAVDWWALGVLIYEMAAGYPPFFADQPIQIYEKIVSGKVRFPSHFSSDLKDLLRNLLQVDLTKRFGNLKDGVNDIKNHK",
                                "FSMNSKEALGGGKFGAVCTCTEKSTGLKLAAKVIKKQTPKDKEMVMLEIEVMNQLNHRNLIQLYAAIETPHEIVLFMEYIEGGELFERIVDEDYHLTEVDTMVFVRQICDGILFMHKMRVLHLDLKPENILCVNTTGHLVKIIDFGLARRYNPNEKLKVNFGTPEFLSPEVVNYDQISDKTDMWSLGVITYMLLSGLSPFLGDDDTETLNNVLSGNWYFDEETFEAVSDEAKDFVSNLIVKEQGARMSAAQCLAHPWLNNL",
                                "AKYDIKALIGRGSFSRVVRVEHRATRQPYAIKMIETKYREGREVCESELRVLRRVRHANIIQLVEVFETQERVYMVMELATGGELFDRIIAKGSFTERDATRVLQMVLDGVRYLHALGITHRDLKPENLLYYHPGTDSKIIITDFGLASARKKGDDCLMKTTCGTPEYIAPEVLVRKPYTNSVDMWALGVIAYILLSGTMPFEDDNRTRLYRQILRGKYSYSGEPWPSVSNLAKDFIDRLLTVDPGARMTALQALRHPWVVSM",
                                "ANYKRLEKVGEGTYGVVYKALDLRPGQGQRVVALKKIRLESEDEGVPSTAIREISLLKELKDDNIVRLYDIVHSDAHKLYLVFEFLDLDLKRYMEGIPKDQPLGADIVKKFMMQLCKGIAYCHSHRILHRDLKPQNLLINKDGNLKLGDFGLARAFGVPLRAYTHEIVTLWYRAPEVLLGGKQYSTGVDTWSIGCIFAEMCNRKPIFSGDSEIDQIFKIFRVLGTPNEAIWPDIVYLPDFKPSFPQWRRKDLSQVVPSLDPRGIDLLDKLLAYDPINRISARRAAIHPYFQES",
                                "TRFRNVTLLGSGEFSEVFQVEDPVEKTLKYAVKKLKVKFSGPKERNRLLQEVSIQRALKGHDHIVELMDSWEHGGFLYMQVELCENGSLDRFLEEQGQLSRLDEFRVWKILVEVALGLQFIHHKNYVHLDLKPANVMITFEGTLKIGDFGMASVWPVPRGMEREGDCEYIAPEVLANHLYDKPADIFSLGITVFEAAANIVLPDNGQSWQKLRSGDLSDAPRLSSTDNGSSLTSSSRETPANSIIGQGGLDRVVEWMLSPEPRNRPTIDQILATDEVCWV",
                                "SEVMLSTRIGSGSFGTVYKGKWHGDVAVKILKVVDPTPEQFQAFRNEVAVLRKTRHVNILLFMGYMTKDNLAIVTQWCEGSSLYKHLHVQETKFQMFQLIDIARQTAQGMDYLHAKNIIHRDMKSNNIFLHEGLTVKIGDFGLATVKSRWSGSQQVEQPTGSVLWMAPEVIRMQDNNPFSFQSDVYSYGIVLYELMTGELPYSRDQIIFMVGRGYASPDLSKLYKNCPKAMKRLVADCVKKVKEERPLFPQILSSIELLQH",
                                "EQVCLLQRLGAGGFGSVYKATYRGVPVAIKQVNKCTKNRLASRRSFWAELNVARLRHDNIVRVVAASTRTPAGSNSLGTIIMEFGGNVTLHQVIYGAAGHPEGDAGEPHCRTGGQLSLGKCLKYSLDVVNGLLFLHSQSIVHLDLKPANILISEQDVCKISDFGCSEKLEDLLCFQTPSYPLGGTYTHRAPELLKGEGVTPKADIYSFAITLWQMTTKQAPYSGERQHILYAVVAYDLRPSLSAAVFEDSLPGQRLGDVIQRCWRPSAAQRPSARLLLVDLTSLKA",
                                "ESLRLEVKLGQGCFGEVWMGTWNGTTRVAIKTLKPGNMSPEAFLQEAQVMKKLRHEKLVQLYAVVSEEPIYIVTEYMSKGSLLDFLKGEMGKYLRLPQLVDMAAQIASGMAYVERMNYVHRDLRAANILVGENLVCKVADFGLARLIEDNEYTARQGAKFPIKWTAPEAALYGRFTIKSDVWSFGILLTELTTKGRVPYPGMVNREVLDQVERGYRMPCPPECPESLHDLMCQCWRRDPEERPTFEYLQAFLEDYFT",
                                "VLNRAVPKDKWVLNHEDLVLGEQIGRGNFGEVFSGRLRADNTLVAVKSCRETLPPDIKAKFLQEAKILKQYSHPNIVRLIGVCTQKQPIYIVMELVQGGDFLTFLRTEGARLRMKTLLQMVGDAAAGMEYLESKCCIHRDLAARNCLVTEKNVLKISDFGMSREAADGIYAASGGLRQVPVKWTAPEALNYGRYSSESDVWSFGILLWETFSLGASPYPNLSNQQTREFVEKGGRLPCPELCPDAVFRLMEQCWAYEPGQRPSFSAIYQEL",
                                "DQLVLGRTLGSGAFGQVVEATAHGLSHSQATMKVAVKMLKSTARSSEKQALMSELYGDLVDYLHRNKHTFLQRHSNKHCPPSAELYSNALPVGFSLPSHLNLTGESDGGYMDMSKDESIDYVPMLDMKGDIKYADIESPSYMAPYDNYVPSAPERTYRATLINDSPVLSYTDLVGFSYQVANGMDFLASKNCVHRDLAARNVLICEGKLVKICDFGLARDIMRDSNYISKGSTYLPLKWMAPESIFNSLYTTLSDVWSFGILLWEIFTLGGTPYPELPMNDQFYNAIKRGYRMAQPAHASDEIYEIMQKCWEEKFETRPPFSQLVLLLERLLGEGYKKKY",
                                "TEFKKIKVLGSGAFGTVYKGLWIPEGEKVKIPVAIKELREATSPKANKEILDEAYVMASVDNPHVCRLLGICLTSTVQLITQLMPFGCLLDYVREHKDNIGSQYLLNWCVQIAKGMNYLEDRRLVHRDLAARNVLVKTPQHVKITDFGLAKLLGAEEKEYHAEGGKVPIKWMALESILHRIYTHQSDVWSYGVTVWELMTFGSKPYDGIPASEISSILEKGERLPQPPICTIDVYMIMVKCWMIDADSRPKFRELIIEFSKMAR",
                                "MGFTIHGALTPGSEGCVFDSSHPDYPQRVIVKAGWYTSTSHEARLLRRLDHPAILPLLDLHVVSGVTCLVLPKYQADLYTYLSRRLNPLGRPQIAAVSRQLLSAVDYIHRQGIIHRDIKTENIFINTPEDICLGDFGAACFVQGSRSSPFPYGIAGTIDTNAPEVLAGDPYTTTVDIWSAGLVIFETAVHNASLFSAPRGPKRGPCDSQITRIIRQAQVHVDEFSPHPESRLTSRYRSRAAGNNRPPYTRPAWTRYYKMDIDV" };
*/
                        /*
                                >PIR:HAHU Hemoglobin alpha chain - Human
                                >PIR:HBHU Hemoglobin beta chain - Human
                                >PIR:MYHU Myoglobin - Human
                                >IGLOB Haemoglobin (Ctp HbVIIB-5) - Midge (Chironomus thummi)
                                >PIR:GPUGNI Nonlegume hemoglobin I - Swamp oak
                                >PIR:GGZLB Bacterial hemoglobin - Vitreoscilla sp.
                                >PIR:GPYL Leghemoglobin I - Yellow lupine
                        */
/*
        String[] seqs =
                        {
                                "VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
                                "VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH",
                                "GLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASEDLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLEFISECIIQVLQSKHPGDFGADAQGAMNKALELFRKDMASNYKELGFQG",
                                "MKFFAVLALCIVGAIASPLTADEASLVQSSWKAVSHNEVEILAAVFAAYPDIQNKFSQFAGKDLASIKDTGAFATHATRIVSFLSEVIALSGNTSNAAAVNSLVSKLGDDHKARGVSAAQFGEFRTALVAYLQANVSWGDNVAAAWNKALDNTFAIVVPRL",
                                "ALTEKQEALLKQSWEVLKQNIPAHSLRLFALIIEAAPESKYVFSFLKDSNEIPENNPKLKAHAAVIFKTICESATELRQKGHAVWDNNTLKRLGSIHLKNKITDPHFEVMKGALLGTIKEAIKENWSDEMGQAWTEAYNQLVATIKAEMKE",
                                "MLDQQTINIIKATVPVLKEHGVTITTTFYKNLFAKHPEVRPLFDMGRQESLEQPKALAMTVLAAAQNIENLPAILPAVKKIAVKHCQAGVAAAHYPIVGQELLGAIKEVLGDAATDDILDAWGKAYGVIADVFIQVEADLYAQAVE",
                                "GVLTDVQVALVKSSFEEFNANIPKNTHRFFTLVLEIAPGAKDLFSFLKGSSEVPQNNPDLQAHAGKVFKLTYEAAIQLEVNGAVASDATLKSLGSVHVSKGVVDAHFPVVKEAILKTIKEVVGDKWSEELNTAWTIAYDELAIIIKKEMKDAA" };
*/
        String[] seqs =
                        {
                                "VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
                                "VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH",
                                "GLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASEDLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLEFISECIIQVLQSKHPGDFGADAQGAMNKALELFRKDMASNYKELGFQG",
                                "MKFFAVLALCIVGAIASPLTADEASLVQSSWKAVSHNEVEILAAVFAAYPDIQNKFSQFAGKDLASIKDTGAFATHATRIVSFLSEVIALSGNTSNAAAVNSLVSKLGDDHKARGVSAAQFGEFRTALVAYLQANVSWGDNVAAAWNKALDNTFAIVVPRL",
                                "ALTEKQEALLKQSWEVLKQNIPAHSLRLFALIIEAAPESKYVFSFLKDSNEIPENNPKLKAHAAVIFKTICESATELRQKGHAVWDNNTLKRLGSIHLKNKITDPHFEVMKGALLGTIKEAIKENWSDEMGQAWTEAYNQLVATIKAEMKE",
                                "MLDQQTINIIKATVPVLKEHGVTITTTFYKNLFAKHPEVRPLFDMGRQESLEQPKALAMTVLAAAQNIENLPAILPAVKKIAVKHCQAGVAAAHYPIVGQELLGAIKEVLGDAATDDILDAWGKAYGVIADVFIQVEADLYAQAVE",
                                "GVLTDVQVALVKSSFEEFNANIPKNTHRFFTLVLEIAPGAKDLFSFLKGSSEVPQNNPDLQAHAGKVFKLTYEAAIQLEVNGAVASDATLKSLGSVHVSKGVVDAHFPVVKEAILKTIKEVVGDKWSEELNTAWTIAYDELAIIIKKEMKDAA",
                                "VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
                                "VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH",
                                "GLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASEDLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLEFISECIIQVLQSKHPGDFGADAQGAMNKALELFRKDMASNYKELGFQG",
                                "MKFFAVLALCIVGAIASPLTADEASLVQSSWKAVSHNEVEILAAVFAAYPDIQNKFSQFAGKDLASIKDTGAFATHATRIVSFLSEVIALSGNTSNAAAVNSLVSKLGDDHKARGVSAAQFGEFRTALVAYLQANVSWGDNVAAAWNKALDNTFAIVVPRL",
                                "ALTEKQEALLKQSWEVLKQNIPAHSLRLFALIIEAAPESKYVFSFLKDSNEIPENNPKLKAHAAVIFKTICESATELRQKGHAVWDNNTLKRLGSIHLKNKITDPHFEVMKGALLGTIKEAIKENWSDEMGQAWTEAYNQLVATIKAEMKE",
                                "MLDQQTINIIKATVPVLKEHGVTITTTFYKNLFAKHPEVRPLFDMGRQESLEQPKALAMTVLAAAQNIENLPAILPAVKKIAVKHCQAGVAAAHYPIVGQELLGAIKEVLGDAATDDILDAWGKAYGVIADVFIQVEADLYAQAVE",
                                "GVLTDVQVALVKSSFEEFNANIPKNTHRFFTLVLEIAPGAKDLFSFLKGSSEVPQNNPDLQAHAGKVFKLTYEAAIQLEVNGAVASDATLKSLGSVHVSKGVVDAHFPVVKEAILKTIKEVVGDKWSEELNTAWTIAYDELAIIIKKEMKDAA",
                                "VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
                                "VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH",
                                "GLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASEDLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLEFISECIIQVLQSKHPGDFGADAQGAMNKALELFRKDMASNYKELGFQG",
                                "MKFFAVLALCIVGAIASPLTADEASLVQSSWKAVSHNEVEILAAVFAAYPDIQNKFSQFAGKDLASIKDTGAFATHATRIVSFLSEVIALSGNTSNAAAVNSLVSKLGDDHKARGVSAAQFGEFRTALVAYLQANVSWGDNVAAAWNKALDNTFAIVVPRL",
                                "ALTEKQEALLKQSWEVLKQNIPAHSLRLFALIIEAAPESKYVFSFLKDSNEIPENNPKLKAHAAVIFKTICESATELRQKGHAVWDNNTLKRLGSIHLKNKITDPHFEVMKGALLGTIKEAIKENWSDEMGQAWTEAYNQLVATIKAEMKE",
                                "MLDQQTINIIKATVPVLKEHGVTITTTFYKNLFAKHPEVRPLFDMGRQESLEQPKALAMTVLAAAQNIENLPAILPAVKKIAVKHCQAGVAAAHYPIVGQELLGAIKEVLGDAATDDILDAWGKAYGVIADVFIQVEADLYAQAVE",
                                "GVLTDVQVALVKSSFEEFNANIPKNTHRFFTLVLEIAPGAKDLFSFLKGSSEVPQNNPDLQAHAGKVFKLTYEAAIQLEVNGAVASDATLKSLGSVHVSKGVVDAHFPVVKEAILKTIKEVVGDKWSEELNTAWTIAYDELAIIIKKEMKDAA",
                                "VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
                                "VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH",
                                "GLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASEDLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLEFISECIIQVLQSKHPGDFGADAQGAMNKALELFRKDMASNYKELGFQG",
                                "MKFFAVLALCIVGAIASPLTADEASLVQSSWKAVSHNEVEILAAVFAAYPDIQNKFSQFAGKDLASIKDTGAFATHATRIVSFLSEVIALSGNTSNAAAVNSLVSKLGDDHKARGVSAAQFGEFRTALVAYLQANVSWGDNVAAAWNKALDNTFAIVVPRL",
                                "ALTEKQEALLKQSWEVLKQNIPAHSLRLFALIIEAAPESKYVFSFLKDSNEIPENNPKLKAHAAVIFKTICESATELRQKGHAVWDNNTLKRLGSIHLKNKITDPHFEVMKGALLGTIKEAIKENWSDEMGQAWTEAYNQLVATIKAEMKE",
                                "MLDQQTINIIKATVPVLKEHGVTITTTFYKNLFAKHPEVRPLFDMGRQESLEQPKALAMTVLAAAQNIENLPAILPAVKKIAVKHCQAGVAAAHYPIVGQELLGAIKEVLGDAATDDILDAWGKAYGVIADVFIQVEADLYAQAVE",
                                "GVLTDVQVALVKSSFEEFNANIPKNTHRFFTLVLEIAPGAKDLFSFLKGSSEVPQNNPDLQAHAGKVFKLTYEAAIQLEVNGAVASDATLKSLGSVHVSKGVVDAHFPVVKEAILKTIKEVVGDKWSEELNTAWTIAYDELAIIIKKEMKDAA" };

                QAlign.setNumber(seqs.length);
                QAlign qalign = new QAlign();
                qalign.setAll(
                        0,	// weighting tree
                        1,	// output console
                        0,	// simultaneous alignment
                        seqs,
                        CostTable.PAM250,
                        5,	// minimal epsilon
                        null,
                        null,
                        null);

                qalign.run();

                if (qalign.foundProgressive) {

                        System.out.println();
                        System.out.println("PROGRESSIVE alignment found:");
                        String[] result= qalign.getProgressiveLayout();
                        for (int i= 0; i< result.length; ++i)
                                System.out.println(result[i]);
                }
                if (qalign.foundSimultaneous) {

                        System.out.println();
                        System.out.println("SIMULTANEOUS alignment found:");
                        String[] result= qalign.getSimultaneousLayout();
                        for (int i= 0; i< result.length; ++i)
                                System.out.println(result[i]);
                }
        }

        /**
         * Multiple setting of features.
         */
        public void setAll(
                int bf,
                int of,
                int mf,
                String[] newSeqs,
                int newCostType,
                int newMinE,
                int[][] fposSeq,
                int[][] fposStart,
                int[] fposLen) throws CancelException {


                        // process arguments
                setOflag((byte) of);
                setBflag((byte) bf);
                setGflag((byte) 1);	// penalize terminal gaps
                setFflag((byte) 1);	// provide fixed positions
                setMflag((byte) mf);
                setMinE(newMinE);

                        // set input
                setSequences(newSeqs);
                setCostTable(newCostType);
                setForcedPositions(fposSeq, fposStart, fposLen);

        }

        /**
         * Multiple setting of features.
         */
        public void setAll (
                int of,
                int mf,
                String[] newSeqs,
                int newCostType) throws CancelException {

                        // process arguments
                setOflag((byte) of);
                setMflag((byte) mf);

                        // set input
                setSequences(newSeqs);
                setCostTable(newCostType);
        }

        /**
         * Starts the algorithm with the given settings.
         */
        public void run() throws CancelException {

                int i, j, len, size;
                long starttime, endtime;

                starttime = System.currentTimeMillis();

                // *** ALLOCATE MEMORY ***

                for (len = i = 1; i <= K; i++) // find longest sequence
                        if (N[i] > len)
                                len = N[i];
                ++len; // take one more... (1- to 0-based)
                size = len * (1 + eflag);
                for (i = 0; i < LENGTH + 1; i++) { // init rows of 2D-Arrays
                        dd[i] = new int[size];
                        hh[i] = new int[size];
                        vv[i] = new int[size];
                }

                // *** COMPUTE MULTIPLE SEQUENCE ALIGNMENT ***

                if (K == 2) // pairwise alignment
                        bflag = 0; // no need of weights
                if (bflag == 0) // no sequence weights
                        for (i = 1; i < K; ++i)
                                for (j = i + 1; j <= K; ++j) // for all sequence combinations
                                        scale[i][j] = 1; // init scale factor with one (no scaling)

                if ((bflag != 0)
                        || (eflag != 0)) { // if weight factors or epsilons needed
                        if (oflag != 0)
                                System.out.println("Calculating pairwise alignments.");
                        if (proxy!= null)
                                proxy.setMessage("calculating pairwise alignments");
                        primer(); // init pairwise alignments
                }

                if (bflag != 0) { // weight factors needed
                        if (oflag != 0)
                                System.out.println("Calculating weights.");
                        if (proxy!= null)
                                proxy.setMessage("calculating weights");
                        bias(); // ~beeinflussen
                }

                if (eflag != 0) { // epsilons needed

                        if (oflag != 0)
                                System.out.println("Calculating epsilons.");

                        if (proxy!= null)
                                proxy.setMessage("calculating epsilons");

                        try {
                                ecalc(2 * len - 2); // calculate epsilons
                        } catch (Exception e) {
                                if (oflag!= 0)
                                        e.printStackTrace();
                                throw new CancelException(e.getMessage());
                        }

                        foundProgressive= true; // a progressive alignment was found

                        if (oflag != 0) {
                                System.out.println("Estimated epsilons:");
                                System.out.println("sequences\tepsilon");
                                for (i = 1; i < K; ++i)
                                        for (j = i + 1;
                                                j <= K;
                                                ++j) // all combinations of sequences
                                                System.out.println("   " + i + "   " + j + " \t" + epsi[i][j]);
                                System.out.println();
                        }
                        if (proxy!= null) {
                                String outStr= "Estimated epsilons:\nsequences\tepsilon";
                                for (i = 1; i < K; ++i)
                                        for (j = i + 1;
                                                j <= K;
                                                ++j) // all combinations of sequences
                                                outStr+= "   " + i + "   " + j + " \t" + epsi[i][j]+ "\n";
                                outStr+="\n";
                                proxy.setOutput(outStr);
                        }
                }

                if (mflag != 0) { // optimal multiple alignment required
                        if (oflag != 0)
                                System.out.println("Calculating pairwise projection costs.");
                        if (proxy!= null)
                                proxy.setMessage("calculating pairwise projection costs");
                        faces(); // projections of k-space = faces
                        // FREE vv, hh, dd = Garbage.collector?
                        if (delta < 0)
                                for (delta = 0, i = 1; i < K; ++i)
                                        for (j = i + 1;
                                                j <= K;
                                                ++j) // all combinations of sequences
                                                delta += (scale[i][j] * epsi[i][j]);
                        Upper = Lower + delta; // determine upper bound of costs
                        for (i = 0; i < K; ++i)
                                for (j = i + 1; j <= K; ++j) {
                                        proj[i][j] = 0; // init half of the proj[][] with 0
                                        scale[j][i] = scale[i][j]; // make symmetric for msa()
                                }
                        if (oflag != 0)
                                System.out.println("Calculating multiple alignment.");
                        if (proxy!= null)
                                proxy.setMessage("calculating multiple alignment");

                        rootEdge= msa();
                        if (rootEdge==null || rootEdge.dist> Upper) {	// heap ran empty || lastMan not good enough
                                fatal("Multiple alignment within bound does not exist.");
                                return;
                        }
                        foundSimultaneous= true; 	// a simultaneous alignment has been found
                        if (oflag!= 0)
                                display(rootEdge); // display result by edge computed by msa()
                }

                endtime = System.currentTimeMillis();
                deltatime = endtime - starttime;
                if (oflag!= 0)
                        System.out.println("Elapsed time = " + deltatime);
                getProgressiveLayout();
        }

        public void cancel() {

                cancel= true;
        }

        public int getCost() {
        	return alicost;
        }
        
        /**
         * Sets the cost table used for the alignment.<br>
         * Automatically sets terminal gap costs according to the <code>gflag</code>
         *
         * @param type the type of the cost table according to the static fields
         * @see #DNARNA
         * @see #PAM250
         * @see CostTable
         */
        public void setCostTable(int type) {

                        // get cost table
                byte[][] costs= CostTable.getCostTable(type);

                        // load costs into members of MSA
                G= costs[0][0];
//		System.out.println(G);
                for (int i= 1; i< costs.length; ++i)
                        for (int j= 1; j< costs[i].length; ++j) {
                                sub((char) costs[i][0], (char) costs[j][0], costs[i][j]);
//				System.out.println((char) costs[i][0]+" "+ (char) costs[j][0]+" "+ costs[i][j]);
                        }

                GG = (gflag != 0) ? 0 : G; // end gaps not penalized or penalized

                        // Altschul Gap Costs
                        // 1 means DASH
                        // 0 means any char
                T[0][0][0][0] = 0; // -- : --
                T[0][0][0][1] = G; // -- : -x
                T[0][0][1][0] = G; // -x : --
                T[0][0][1][1] = 0; // -x : -x
                T[0][1][0][0] = 0; // -- : x-
                T[0][1][0][1] = 0; // -- : xx
                T[0][1][1][0] = G; // -x : x-
                T[0][1][1][1] = 0; // -x : xx
                T[1][0][0][0] = 0; // x- : --
                T[1][0][0][1] = G; // x- : -x
                T[1][0][1][0] = 0; // xx : --	?
                T[1][0][1][1] = 0; // xx : -x	?
                T[1][1][0][0] = 0; // x- : x-	?
                T[1][1][0][1] = G; // x- : xx	?
                T[1][1][1][0] = G; // xx : x-
                T[1][1][1][1] = 0; // xx : xx
                T[2][0][2][0] = 0; // -- : --
                T[2][0][2][1] = 0; // -- : --
                T[2][1][2][0] = 0; // -- : --
                T[2][1][2][1] = 0; // -- : --
                T[0][2][0][2] = 0; // -- : --
                T[0][2][1][2] = 0; // -- : --
                T[1][2][0][2] = 0; // -- : --
                T[1][2][1][2] = 0; // -- : --
                T[2][2][2][2] = 0; // -- : --

        }

        /**
         * Sets the epsilon values.<br>
         * <br>
         * Epsilon values should be provided with  one
         * integer  for  each  pair  of sequences in the order
         * 1-2, 1-3, ... , 1-N, 2-3, ... , (N-1)-N.
         *
         * @param newEpsilons the epsilon values for all sequences
         * @deprecated not tested
         */
        public void setEpsilons(int[] newEpsilons) {

                        // nothing to do
                if (newEpsilons.length!= ((K*(K-1))/2)) {
                        System.err.println("Invalid number of epsilons in setEpsilons(int[])!");
                        return;
                }

                        // epsilons dont need to be calculated
                eflag= 0;

                        // load into members
                int count= 0;
                for (int i = 1; i < K; ++i)
                        for (int j = (i+1); j <= K; ++j) {

                                epsi[i][j] = newEpsilons[count++];
                                if (epsi[i][j] < 0)
                                        fatal("Epsilon must be positive.");
                        }

        }
        /**
         * Sets the input sequences.<br>
         * <br>
         * <i>(Adapted from data())</i>
         *
         * @param seqs the set of input sequences.
         */
        public void setSequences(String[] seqs) throws CancelException {

                int i, j;
                int n, nb;
                char a;
                char[][] buffer;;

                        // read sequences
                nb = seqs.length - 1;

                if ((nb + 1) > NUMBER) {
                        if (oflag!= 0)
                                System.out.println((nb + 1) + " is too many seqs");
                        throw new CancelException((nb + 1) + " is too many seqs");
                }

                        // load internal arrays
                buffer= new char[nb+2][];
                for (K = 1; K <= nb + 1; K++) {
                        buffer[K] = new char[seqs[K - 1].length() + 1];
                        for (n = 0, i = 0; i < seqs[K - 1].length(); i++) {
                                a = seqs[K - 1].charAt(i);
                                if (n > LENGTH + 1) {
                                        if (oflag!= 0)
                                                System.out.println("seq " + K + " too long (QAlign)");
                                        throw new CancelException("seq " + K + " too long ");
                                }
                                //if ((a >= 'A') && (a <= 'Z')) { // cut wrong characters
                                        buffer[K][i] = a; // convert to char[]
                                        n++;
                                //}
                        }
                        N[K] = n; // save in length array
                        buffer[K][i] = '\0'; // append '\0' to end of buffer
                        S[K] = new char[n + 2]; // initial DASH and terminal '\0' ??
                        S[K][0] = getDASH(); // append '-' to start of seq
                        for (i = 0; i < n + 1; i++)
                                S[K][i + 1] = buffer[K][i];
                }
                K--;

                        // test out
/*			for (i= 0; i< seqsStr.length; i++) {
                                for (j= 0; j< S[i].length; j++)
                                        System.out.print(S[i][j]);
                                System.out.print(" : "+N[i]);
                                System.out.println();
                        }
*/

        }

        /**
         * returns new String[] with pos 0..K-1 inited.
         * Erstellungsdatum: (20.02.2002 18:13:44)
         * @return java.lang.String[]
         */
        public String[] getProgressiveLayout() {

                        // see ecalc() for stopper init (from original code!)
                byte stopper= (byte) '*';
                int plength;
                for (plength= (AL[0].length-1); plength> 0; --plength)
                        if (AL[0][plength]== stopper)
                                break;

                        // convert to string[]
                String[] result= new String[K];
                for (int i= 0; i< result.length; ++i) {
                        result[i]= new String(AL[i+1],1,(plength-1));
//			System.err.println(result[i]);
                }
                return result;
        }

        public String[] getSimultaneousLayout() throws CancelException {

                        // nothing to do
                if (!foundSimultaneous)
                        return null;

                Edge e, f;
                int[] p= new int[NUMBER+1];
                int[] q= new int[NUMBER+1];
                int[] r= new int[NUMBER+1];
                e= rootEdge;

                        // recover shortest path to source tracing backward from sink
                for (e.next=null;e.tail!=presource;e=f) {
                f=e.backtrack;
                        coord(f.tail,p);
                        safe_coord(f.head,q);
                        safe_coord(e.head,r);
                        f.next = e;
                        project(p,q,r);
                }

                        // project path to alignment
                StringBuffer[] result= new StringBuffer[K];
                for (int i= 0; i< result.length; ++i)
                        result[i]= new StringBuffer();
                for (e=e.next; e!=null; e=e.next) {
                        coord(e.tail,p);
                        safe_coord(e.head,q);
                        for (int k=1;k<=K;k++)
                                result[k-1].append((char) S[k][q[k] *	(q[k]-p[k])]);
                }

                        // convert to string
                String[] res= new String[result.length];
                for (int i= 0; i< result.length; ++i)
                        res[i]= result[i].toString();

                return res;
        }

        /**
         * Sets the positions to be aligned fix in <b>progressive</b> alignment only!<br>
         * <br>
         * Example: to force<br>
         * <br>
         *   * **<br>
         * TAT-TT<br>
         * T-TATT<br>
         * <br>
         * one has to set:
         * seqNum= {{1,2}, {1,2}},
         * seqPos= {{3,2}, {4,4}} and
         * seqLen= {1,2}.<br>
         * <br>
         * <i>(Adapted from <code>fixread()</code>)<i>
         *
         * @param seqNum all arrays with the 1-based numbers of the sequences involved
         * in each fixed block
         * @param pos all arrays with the 1-based positions within each of the sequences
         * (in the same order as in parameter <code>seqNum</code>
         * @param length the lengths of the fixed blocks (in the correct order with respect
         * to <code>seqNum</code> and <code>pos</code>
         */
        public void setForcedPositions(
                int[][] seqNum,
                int[][] seqPos,
                int[] length) {

                int i, j, len, sqnum;
                int[] sqnames = new int[NUMBER];
                int[] sqstarts = new int[NUMBER];

                        // init fix array
                for (i = 1; i <= NUMBER; ++i)
                        for (j = 1; j <= NUMBER; ++j)
                                pFIX[i][j] = null;

                        // abort if not inited
                if ((seqNum == null) || (seqPos == null) || (length == null))
                        return;


                        // copy to members
                for (int x = 0; x < seqNum.length; x++) {
                        if (seqNum[x].length != seqPos[x].length) // error
                                fatal("Block starts do not equal sequence number.");

                        for (i = 0;
                                i < seqNum[x].length;
                                ++i) { // init program handled variables
                                sqnames[i] = seqNum[x][i];
                                sqstarts[i] = seqPos[x][i];
                        }
                        sqnum = sqnames.length;
                        len = length[x];

                        fill(sqnum, sqnames, sqstarts, len);
                }

        }

        /**
         * Switches between progressive and simultaneous multiple
         * alignment.<br>
         * The progressive alignment and other data always is produced
         * by the program before the it attempts to produce a simultaneous
         * multiple alignment.<br>
         * The default mode (<code>mflag= 1</code>) turns on the simultaneous
         * alignment.
         * @param mflag The mflag to set
         */
        public void setMflag(byte mflag) {
                this.mflag = mflag;
        }


        /**
         * Switches the simultaneous alignment (after aligning progressively)
         * on and off.
         * @param newSimultanousAlignment <code>true</code> if a simultaneous
         * alignment is desired, <code>false</code> else.
         * @see #mflag
         * @see #setMflag(byte)
         */
        public void setSimultaneousAlignment(boolean newSimultaneousAlignment) {

                if (newSimultaneousAlignment)
                        setMflag((byte) 1);
                else
                        setMflag((byte) 0);
        }

        /**
         * Sets the gflag.<br>
         * default 1: do not penalize terminal gaps.<br>
         * (-g)= 0= end gaps PENALIZED (charged the same as internal ones)<br>
         * @param gflag The gflag to set
         * @see #gflag
         */
        public void setGflag(byte gflag) {
                this.gflag = gflag;
        }

        /**
         * Sets the fflag.<br>
         * 0: (default) no fixed positions.
         * 1: file with forced positions provided.
         * @param fflag The fflag to set
         */
        public void setFflag(byte fflag) {
                this.fflag = fflag;
        }

        /**
         * Sets the oflag.<br>
         * default 1: output to console.<br>
         * set to 0: suppress output.
         * @param oflag The oflag to set
         */
        public void setOflag(byte oflag) {
                this.oflag = oflag;
        }

        /**
         * Sets the oflag.<br>
         * @param newOutput <code>true</code> if output to console should be enabled,
         * <code>false</code> else
         * @see #setOflag(byte)
         */
        public void setOutput(boolean newOutput) {

                if (newOutput)
                        setOflag((byte) 1);
                else
                        setOflag((byte) 0);
        }

        /**
         * It is strongly recommended not to use a simultaneous alignment
         * on a bigger set.<br>
         * Trying <code>java -Xmx512M</code> and hoping the best.
         */
        public static void setNumber(int newNumber) throws NumberFormatException {

                ++newNumber;	// 1-based

                        // accept only byte values
                if ( (newNumber< 0) )
                        throw new NumberFormatException("Too many sequences.");

                        // nothing to do:
                        // could conflict with multiple instances
                if (newNumber< NUMBER)
                        return;

                NUMBER= newNumber;
        }

        /**
         * It is strongly recommended not to use a simultaneous alignment
         * on a bigger set.<br>
         * Trying <code>java -Xmx512M</code> and hoping the best.
         */
        public static void setLength(int newLength) {

                newLength+= 2;	// 1-based and 0-terminated

                        // accept only byte values
                if ((newLength> Short.MAX_VALUE)
                        || (newLength< 0))
                        throw new NumberFormatException("Not a valid length: "+newLength);

                        // nothing to do:
                        // could conflict with multiple instances
                if (newLength< LENGTH)
                        return;

                LENGTH= (short) newLength;
        }

        /**
         * Sets the bflag.<br>
         * default 1: tree for weights is calculated.
         * set to 0: unweighted summing (DO NOT calculates evolutionary tree).
         * @param bflag The bflag to set
         * @deprecated Progressive alignment may crash (same as in original c-code).
         */
        public void setBflag(byte bflag) {
                this.bflag = bflag;
        }

        /**
         * Sets the bflag.<br>
         * default 1: tree for weights is calculated.
         * set to 0: unweighted summing (DO NOT calculates evolutionary tree).
         * @param newWeightingTree <code>true</code> if tree weights for the
         * alignment are desired, <code>false</code> else.
         * @deprecated Progressive alignment may crash (same as in original c-code).
         */
        public void setWeightingTree(boolean newWeightingTree) {

                if (newWeightingTree)
                        setBflag((byte) 1);
                else
                        setBflag((byte) 0);
        }

        /**
         * Returns the time (in [msec]) the <code>run()</code> method needed to solve the problem.
         */
        public long getTime() {

                return deltatime;
        }

        /**
         * Sets the heuristical minimum boundary.<br>
         * Empirically, <i>Gupta et al.</i> found 5 to be a good value,
         * which is used by default.<br>
         * @param newMinE the new minimal boundary of the Lippman heuristics
         * @deprecated WARNING: changing the default may be changing results!
         */
        public void setMinE(int newMinE) {

                this.MINE= newMinE;
        }

        /**
         * Set's the proxy for different GUI output.
         * @param the new GUI proxy
         */
        public void setProxy(GUIProxy newProxy) {

                this.proxy= newProxy;
        }
}
