/*
 * Created on Mar 6, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.db;

import gphase.Constants;

import gphase.algo.ASAnalyzer;
import gphase.io.gtf.GTFObject;
import gphase.model.ASMultiVariation;
import gphase.model.ASVariation;
import gphase.model.Exon;
import gphase.model.Gene;
import gphase.model.GeneHomology;
import gphase.model.Graph;
import gphase.model.EncodeRegion;
import gphase.model.GraphHandler;
import gphase.model.Phase;
import gphase.model.Species;
import gphase.model.Transcript;
import gphase.model.Translation;

import java.beans.Encoder;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.lang.reflect.Method;
import java.sql.Connection;
import java.sql.DatabaseMetaData;
import java.sql.SQLException;
import java.sql.Statement;
import java.sql.ResultSet;
import java.sql.ResultSetMetaData;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.Properties;
import java.util.StringTokenizer;
import java.util.Vector;

import org.ensembl.driver.ConfigurationException;
import org.ensembl.driver.Driver;
import org.ensembl.driver.DriverManager;
import org.ensembl.driver.EnsemblDriver;
import org.ensembl.driver.compara.HomologyAdaptor;
import org.ensembl.driver.plugin.standard.BaseDriver;
import org.ensembl.util.PropertiesUtil;

import com.sun.org.apache.bcel.internal.generic.GETSTATIC;

import qalign.tools.FASTAWrapper;
import qalign.tools.SequenceWrapper;

/**
 * 
 * 
 * @author micha
 */
// TODO: 54 threads open after reading in 3 species !!!
public class EnsemblDBAdaptor {

	public static boolean DEBUG= false;
	Driver driver;
	static int ensemblVersion= -1;
	
	Graph graph= null;
	HashMap tempTranscriptMap= null, tempExonMap= null;
	public static final String[] SPECIES_ALL_SHORT = { "agambiae", "amellifera", "celegans", "cfamiliaris", "dmelanogaster", "drerio", "frubripes", "hsapiens", "mmusculus", "ptroglodytes", "rnorvegicus", "scerevisiae", "tnigroviridis" };
	public static final String[] SPECIES_ENCODE = { Species.SP_NAMES_BINOMIAL[0], Species.SP_NAMES_BINOMIAL[1], Species.SP_NAMES_BINOMIAL[2], Species.SP_NAMES_BINOMIAL[3], Species.SP_NAMES_BINOMIAL[4], Species.SP_NAMES_BINOMIAL[5], Species.SP_NAMES_BINOMIAL[6], Species.SP_NAMES_BINOMIAL[7], Species.SP_NAMES_BINOMIAL[8], Species.SP_NAMES_BINOMIAL[9], Species.SP_NAMES_BINOMIAL[10], Species.SP_NAMES_BINOMIAL[11] };
	public static final String[] SPECIES_ENSEMBL = { Species.SP_NAMES_BINOMIAL[0], Species.SP_NAMES_BINOMIAL[1], Species.SP_NAMES_BINOMIAL[2], Species.SP_NAMES_BINOMIAL[3], Species.SP_NAMES_BINOMIAL[4], Species.SP_NAMES_BINOMIAL[5], Species.SP_NAMES_BINOMIAL[6], Species.SP_NAMES_BINOMIAL[7], Species.SP_NAMES_BINOMIAL[8], Species.SP_NAMES_BINOMIAL[9], Species.SP_NAMES_BINOMIAL[10], Species.SP_NAMES_BINOMIAL[11] };
	public static final String[] SPECIES_SMALL = { "mouse", "yeast" };
	public static final String[] SPECIES_ISMB = { "human", "mouse", "rat", "cow", "dog", "chicken", "frog", "zebrafish", "fruitfly", "mosquito" };
	public static final String[] SPECIES = SPECIES_ISMB;

	/**
	 * 
	 */
	public EnsemblDBAdaptor() {
		super();
		// TODO Auto-generated constructor stub
	}

	/**
	*Stellt die Verbindung zur Datenbank her.
	*@param url Die URL der Datenbank <i>ohne die Angabe von "jdbc:mysql://"</i>
	*@param user Der Benutzername
	*@param pass Das Passwort
	*@return Die hergestellte Connection con
	*/
	public Connection connect(String dbName) {
 
		Properties prop= new Properties();
		prop.put("host","ensembldb.ensembl.org");		// guaranteed to point @ latest db
		
		prop.put("user","anonymous");		// guaranteed to point @ latest db
		prop.put("database_prefix",dbName);		// guaranteed to point @ latest db
//		prop.put("database",dbName);		// guaranteed to point @ latest db
//		prop.put("ensembl_driver","org.ensembl.driver.plugin.compara.ComparaMySQLDriver");
//		prop.put("connection_pool_size","4");
//		prop.put("jdbc_driver","com.p6spy.engine.spy.P6SpyDriver");
//		prop.put("port","3333");
//		prop.put("jdbc_driver","com.p6spy.engine.spy.P6SpyDriver");
		
		Connection connection= null;

		driver= null;	// re-init doesnt work: ConcurrentModificationException
						// killed startCloserThread in ConnectionPoolDataSource 
						// to not get too many closing threads (>100)
		if (driver== null) 
			try {driver = DriverManager.loadDriver(prop);}
			catch (ConfigurationException e) {e.printStackTrace();}
		else {
			try {	
				driver.closeAllConnections();		// these guys throw some exceptions, dont worry
			} catch (Exception e) {; /*:-)*/}
			try {	
				driver.removeAllAdaptors();
			} catch (Exception e) {; /*:-)*/}
			try {	
				driver.clearAllCaches();
			} catch (Exception e) {; /*:-)*/}
			try {	
				driver.initialise(prop);
			} catch (Exception e) {; /*:-)*/}
		}
		try {
			connection = driver.getConnection();
//			System.out.println("Driver connection: " + driver.toString());
//			System.out.println("Driver configuration: " +driver.getConfiguration());
		} catch (Exception e) {
			System.err.println("connect() failed:");
			e.printStackTrace();
//			System.out.println("SQLState:     " + e.getSQLState());
//			System.out.println("VendorError:  " + e.getErrorCode());
		}

		return connection;
	}
	
	
	/**
	*@deprecated delme
	*/
	public Connection connect2(String dbName) {

		Connection connection= null;

		try {					
			Class.forName("org.gjt.mm.mysql.Driver").newInstance();	// loading the driver
				// org.ensembl.driver.plugin.standard.MySQLDriver
			
		} catch (Exception e) {
			e.printStackTrace();
		}

		try {					
			connection =						// logging into database
				java.sql.DriverManager.getConnection(
					"jdbc:mysql://" + "ensembldb.ensembl.org",
					"anonymous",
					"");
		} catch (Exception e) {
			System.err.println("Fehler bei Connection");
			System.err.println("SQLException: " + e.getMessage());
//			System.out.println("SQLState:     " + e.getSQLState());
//			System.out.println("VendorError:  " + e.getErrorCode());
		}

		return connection;
	}
	
	/**
	 * @deprecated delme ?!
	 */
	/*
	 * Modification in org.ensembl.plugin.driver.standard.MySQLAdaptor.loadAdaptors();
	 */
	public Connection connect3(String dbName) {

		Connection connection= null;

		try {					
			EnsemblDriver driver = (EnsemblDriver) Class.forName("org.ensembl.driver.plugin.standard.MySQLDriver").newInstance();
			Properties driverConfig= new Properties();
			driverConfig.put(
			  "ensembl_driver",
			  "org.ensembl.driver.plugin.standard.MySQLDriver");
			driverConfig.put("host","ensembldb.ensembl.org");
			driverConfig.put("user","anonymous");
			driverConfig.put("database_prefix",dbName);
			driver.initialise(driverConfig);
//			Class.forName("org.gjt.mm.mysql.Driver").newInstance();	// loading the driver
				// 
			
			connection =						// logging into database
				driver.getConnection();
		} catch (Exception e) {
			e.printStackTrace();
		}


		return connection;
	}	
	
	void testStatement(Connection con) {
		
		try {
			Statement stmt = con.createStatement(
						ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
											// only one ResultSet per Statement possible
			ResultSet rs = stmt.executeQuery(
				"select seq_region.name from seq_region,gene where gene.seq_region_start=46192 "
				+ "and seq_region.seq_region_id=gene.seq_region_id");
											// TYPE and CONCUR value given by Statement 
			ResultSetMetaData rsmd= rs.getMetaData();
			int cols= rsmd.getColumnCount();

			for (int i= 1; i<= cols; i++)	// get header of table
				System.out.print(rsmd.getColumnName(i)+ "\t");
			System.out.println("\n--");
		
			while (rs.next()) {				// get body of table
				for (int i= 1; i<= cols; i++)
					System.out.println(rs.getString(i)+ "\t");
			}		
		} catch (SQLException e) { 			// thrown by both methods
			System.err.println(
				"!PROBLEM: "
					+ e
					+ "\n\r"
					+ "\tA database error occured during the query.\n\r"
					+ "\t\tvendor-specific ErrorCode: "
					+ e.getErrorCode()
					+ "\n\r"
					+ "--------------------------------------------------");
		}
		
	}
	
	/**
	 * Checks version of the databases used, e.g.,
	 * homo_sapiens_core_31_35d, ensembl_compara_35, 
	 * ensembl_mart_35.
	 * @param x
	 */
	
	void checkVersion(int x) {
		
		if(ensemblVersion< 0) // not yet inited
			ensemblVersion= x;
		else if (ensemblVersion!= x) {
			System.err.println("Ensembl version error: "+ x+ " != "+ensemblVersion);
			System.exit(-1);
		}
			
	}
	
	Graph buildGraph(String[] specNames) {
		
		if (specNames== null || specNames.length< 1)	// or (< 2) ?
			return null;
		
		System.out.println("["+new Date(System.currentTimeMillis())+"] -- start building graph --");
		Species[] spec= new Species[specNames.length];
		for (int i = 0; i < spec.length; i++) 
			spec[i]= new Species(specNames[i]);
			
		graph= new Graph(spec);
		
		System.out.println("["+new Date(System.currentTimeMillis())+"] "+graph.getSpecies().length+ " species.");
		GeneHomology[] homols= retrieveHomologGenesAll_mart(graph);
		System.out.println("["+new Date(System.currentTimeMillis())+"] retrieving info for "+homols.length+" homologies.");
		homols= retrieveHomologyInfo_single(homols);			// cannot join two databases (mart, compara)
		
		System.out.println("["+new Date(System.currentTimeMillis())+"] retrieving info for "+graph.getGenes().length+" genes.");
		graph= retrieveGeneInfo(graph);
		System.out.println("["+new Date(System.currentTimeMillis())+"] # finished building graph #");
		
		return graph;
	}
	
	
	
	void retrieveExonsSingle(Connection con, Transcript[] transcr) {
		
		ResultSet rs= null;
		int pctCtr= 0;
		int pctBase= transcr.length/ 100;
		if (tempExonMap== null) 	// get or create transcript map
			tempExonMap= new HashMap(tempTranscriptMap.size()* 5, tempTranscriptMap.size()* 7);	// ?5,7
		for (int i= 0; i < transcr.length; i++) {
			Transcript trans= transcr[i];
			try {
				Statement stmt = con.createStatement(
							ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
												// only one ResultSet per Statement possible
				String query= "SELECT exon_stable_id.stable_id,exon.seq_region_start,exon.seq_region_end,phase,end_phase "+
								"FROM exon,exon_stable_id,exon_transcript,transcript,transcript_stable_id "+ 
								"WHERE exon_stable_id.exon_id=exon.exon_id"+
								" AND transcript.transcript_id=exon_transcript.transcript_id"+
								" AND exon.exon_id=exon_transcript.exon_id"+
								" AND transcript.transcript_id=transcript_stable_id.transcript_id"+
								" AND transcript_stable_id.stable_id=\""+ trans.getStableID()+ "\"";

				rs = stmt.executeQuery(query);	// TYPE and CONCUR value given by Statement

				if (rs== null) {
					System.err.println("No exons for transcript "+ trans.getStableID());
					continue;
				}

				while(rs.next()) {
					
					int start= Integer.parseInt(rs.getString(2));
					int end= Integer.parseInt(rs.getString(3));
					Exon e= trans.getGene().getExon(start, end); // exon exists?
					
					if (e== null) {						
						e= new Exon(trans, rs.getString(1), start, end);	// stableID
						// e.setGene(trans.getGene());
						e.setPhase(Byte.parseByte(rs.getString(4))); 
								//Byte.parseByte(rs.getString(5)))	// end phase
						
					} else
						e.addTranscript(trans);	
					
					trans.addExon(e);
//					if (!e.checkStrand(rs.getString(6)))
//						System.err.println("Strand mismatch of transcript "+trans.getStableID()
//							+" with exon "+ e.getStableID());
					
					tempExonMap.put(rs.getString(1), e);
				}
				
				if (pctBase== 0) {
					int k= (100/ transcr.length);
					for (int j = 0; j < k; j++) 
						System.out.print("*");						
				} else if (i/ pctBase> pctCtr) {
					pctCtr++;
					if (pctCtr% 2== 0) 
						System.out.print("*");
				}
			} catch (SQLException e) { 			// thrown by both methods
				System.err.println("Problem fetching exons (single): "+ e);
			}
		}
		System.out.println();
	}
	
	void retrieveTranslations(Connection con, Transcript[] transcr) {
		
		ResultSet rs= null;

		int pctCtr= 0;
		int pctBase= transcr.length/ 100;
		if (tempExonMap== null) 	// get or create transcript map
			tempExonMap= new HashMap(tempTranscriptMap.size()* 5, tempTranscriptMap.size()* 7);
		for (int i= 0; i < transcr.length; i++) {
			Transcript trans= transcr[i];
			try {
				Statement stmt = con.createStatement(
							ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
												// only one ResultSet per Statement possible
					// skippped: ,start_exon_id,end_exon_id
				String query= "SELECT translation_stable_id.stable_id,translation.seq_start,translation.seq_end,ex1.seq_region_start,ex2.seq_region_start "+
								"FROM translation,translation_stable_id,transcript,transcript_stable_id,exon as ex1,exon as ex2 "+	// ,exon_id,exon_stable_id 
								"WHERE translation_stable_id.translation_id=translation.translation_id"+
								" AND transcript.transcript_id=translation.transcript_id"+
								" AND transcript.transcript_id=transcript_stable_id.transcript_id"+
//								" AND exon.exon_id=exon_stable_id.exon_id"+
//								" AND exon.exon_id=translation.start_exon_id"+
//								" AND exon.exon_id=translation.end_exon_id"+
								" AND transcript_stable_id.stable_id=\""+ trans.getStableID()+ "\""+
								" AND translation.start_exon_id=ex1.exon_id AND translation.end_exon_id=ex2.exon_id";

				rs = stmt.executeQuery(query);	// TYPE and CONCUR value given by Statement

				if (rs== null) {
					System.err.println("No translations for transcript "+ trans.getStableID());
					continue;
				}

				while(rs.next()) {
					
					Translation l= new Translation(trans, rs.getString(1));	// stableID)
					//l.setTranscript(trans);
					l.setStart(Integer.parseInt(rs.getString(2))+ Integer.parseInt(rs.getString(4)));
					l.setEnd(Integer.parseInt(rs.getString(3))+ Integer.parseInt(rs.getString(5)));
					trans.addTranslation(l);
				}
				
				// do not close connection here - needed for exons
			} catch (SQLException e) { 			// thrown by both methods
				System.err.println("Problem fetching translations (single): "+ e);
			}
		}
	}	

		/**
		 * 
		 * @param species
		 * @param geneVec
		 * 
		 * @deprecated too slow to demand each gene in single at database
		 */
		void retrieveTranscripts_old(String species, Vector geneVec) {
			
			Connection con= connect(species+ "_core");	// reconnect to species-db
			ResultSet rs= null;
			System.out.println("getting transcripts for "+ species+"("+geneVec.size()+" genes)");
			System.out.println("0    10   20   30   40   50   60   70   80   90  100[%]");
			System.out.println("|----|----|----|----|----|----|----|----|----|----|");
			System.out.print("*");
			System.out.flush();
			int pctCtr= 0;
			int pctBase= geneVec.size()/ 100;
			long t0= System.currentTimeMillis();
			for (int i= 0; i < geneVec.size(); i++) {
				Gene gene= (Gene) geneVec.elementAt(i);
	//			long t1= System.currentTimeMillis();
				try {
					Statement stmt = con.createStatement(
								ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
													// only one ResultSet per Statement possible
					String query= "SELECT transcript_stable_id.stable_id,transcript.seq_region_start,transcript.seq_region_end,transcript.seq_region_strand "+
									"FROM gene,gene_stable_id,transcript,transcript_stable_id "+ 
									"WHERE transcript.gene_id=gene.gene_id"+
									" AND gene.gene_id=gene_stable_id.gene_id"+
									" AND gene_stable_id.stable_id=\""+ gene.getStableID()+ "\""+
									" AND transcript.transcript_id=transcript_stable_id.transcript_id ";
	
					rs = stmt.executeQuery(query);	// TYPE and CONCUR value given by Statement
				 	if (rs== null) {
				 		System.err.println("NO RESULT: "+ gene.getSpecies()+gene.getGeneID());
						continue;
					}
					
	//				long t2= System.currentTimeMillis();
				 
					while(rs.next()) {
						Transcript transcript= new Transcript(gene, rs.getString(1));	// stableID
						// transcript.setGene(gene);
						transcript.setStart(Integer.parseInt(rs.getString(2)));
						transcript.setEnd(Integer.parseInt(rs.getString(3)));
						transcript.checkStrand(rs.getString(4));
						gene.addTranscript(transcript);
					}
					
					if (i/ pctBase> pctCtr) {
						pctCtr++;
						if (pctCtr% 2== 0) {
							System.out.print("*");
							System.out.flush();
						}
					}
	//				long t3= System.currentTimeMillis();
					
	//				System.out.println((float) (t2-t1)+ ","+ (float) (t3-t2));
				} catch (SQLException e) { 			// thrown by both methods
					System.err.println(
						"!PROBLEM: "
							+ e
							+ "\n\r"
							+ "\tA database error occured during the query.\n\r"
							+ "\t\tvendor-specific ErrorCode: "
							+ e.getErrorCode()
							+ "\n\r"
							+ "--------------------------------------------------");
				}
			}
			
			System.out.println("done. ("+ ((System.currentTimeMillis()- t0)/ 1000) +" sec)");
			
		}
	
	Transcript[] retrieveTranscripts(Connection con, Species spec, Gene[] genes) {
			
			ResultSet rs= null;
		
				// build query
			StringBuffer sb= new StringBuffer("SELECT transcript_stable_id.stable_id,transcript.seq_region_start,transcript.seq_region_end,transcript.seq_region_strand,gene_stable_id.stable_id,transcript.biotype,transcript.status ");	//
			sb.append("FROM gene,gene_stable_id,transcript,transcript_stable_id "); 
			sb.append("WHERE transcript.gene_id=gene.gene_id");
			sb.append(" AND gene.gene_id=gene_stable_id.gene_id");
			sb.append(" AND transcript.transcript_id=transcript_stable_id.transcript_id ");
			Gene gene= genes[0];
			sb.append(" AND (gene_stable_id.stable_id=\""); 
			sb.append(gene.getStableID()); 
			sb.append("\"");
			for (int i= 1; i < genes.length; i++) {
				gene= genes[i];
				sb.append(" OR gene_stable_id.stable_id=\"");
				sb.append(gene.getStableID());
				sb.append("\"");
			}
			sb.append(")");
		
				// execute query
			Vector transVec= new Vector();
			try {
				Statement stmt = con.createStatement(
							ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
												// only one ResultSet per Statement possible
		
				rs = stmt.executeQuery(sb.toString());	// TYPE and CONCUR value given by Statement
				if (rs== null) {
					System.err.println("No result for Transcripts!");
					return null;
				}
					
					// store result
				if (tempTranscriptMap== null) {
					int ctr= countResults(rs);
					tempTranscriptMap= new HashMap(ctr, (int) (1.3* ctr));
					rs.beforeFirst();
				}
				while(rs.next()) {
					gene= spec.getGene(rs.getString(5));
					Transcript transcript= new Transcript(gene, rs.getString(1));	// stableID
					//transcript.setGene(gene);		// stableID
					transcript.setStart(Integer.parseInt(rs.getString(2)));
					transcript.setEnd(Integer.parseInt(rs.getString(3)));
					transcript.setType(rs.getString(6));
					transcript.setConfidence(rs.getString(7));
					
						
					if (gene== null) {
						System.err.println("EnsemblDBAdaptor.getTranscripts(): no gene found for "+ rs.getString(5));
						continue;
					}
					gene.addTranscript(transcript);
					
					if (!transcript.checkStrand(rs.getString(4)))
						System.err.println("Strand mismatch gene "+gene.getStableID()
							+"("+(transcript.getGene().isForward()?"fwd":"rev")
							+") with transcript "+ transcript.getStableID()
							+"("+(rs.getString(4).trim().equals("1")?"fwd":"rev"));
					
					tempTranscriptMap.put(rs.getString(1), transcript);
					transVec.add(transcript);
				}
				
				// do not close - connection is needed for translations and exons
			} catch (SQLException e) { 			// thrown by both methods
				System.err.println("Error fetching transcripts: "+ e);
			}
			
			return Transcript.toTranscriptArray(transVec);
		}

	static GTFObject[] getTranslations(Connection con, Species spec, Gene[] genes) {
		
		ResultSet rs= null;

			// build query
		StringBuffer sb= new StringBuffer("SELECT transcript_stable_id.stable_id,transcript.seq_region_start,transcript.seq_region_end,transcript.seq_region_strand,gene_stable_id.stable_id,transcript.biotype,transcript.status ");	//
		sb.append("FROM gene,gene_stable_id,transcript,transcript_stable_id "); 
		sb.append("WHERE transcript.gene_id=gene.gene_id");
		sb.append(" AND gene.gene_id=gene_stable_id.gene_id");
		sb.append(" AND transcript.transcript_id=transcript_stable_id.transcript_id ");
		Gene gene= genes[0];
		sb.append(" AND (gene_stable_id.stable_id=\""); 
		sb.append(gene.getStableID()); 
		sb.append("\"");
		for (int i= 1; i < genes.length; i++) {
			gene= genes[i];
			sb.append(" OR gene_stable_id.stable_id=\"");
			sb.append(gene.getStableID());
			sb.append("\"");
		}
		sb.append(")");
	
			// execute query
		Vector transVec= new Vector();
		try {
			Statement stmt = con.createStatement(
						ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
											// only one ResultSet per Statement possible
	
			rs = stmt.executeQuery(sb.toString());	// TYPE and CONCUR value given by Statement
			if (rs== null) {
				System.err.println("No result for Transcripts!");
				return null;
			}
				
				// store result
			if (tempTranscriptMap== null) {
				int ctr= countResults(rs);
				tempTranscriptMap= new HashMap(ctr, (int) (1.3* ctr));
				rs.beforeFirst();
			}
			while(rs.next()) {
				gene= spec.getGene(rs.getString(5));
				Transcript transcript= new Transcript(gene, rs.getString(1));	// stableID
				//transcript.setGene(gene);		// stableID
				transcript.setStart(Integer.parseInt(rs.getString(2)));
				transcript.setEnd(Integer.parseInt(rs.getString(3)));
				transcript.setType(rs.getString(6));
				transcript.setConfidence(rs.getString(7));
				
					
				if (gene== null) {
					System.err.println("EnsemblDBAdaptor.getTranscripts(): no gene found for "+ rs.getString(5));
					continue;
				}
				gene.addTranscript(transcript);
				
				if (!transcript.checkStrand(rs.getString(4)))
					System.err.println("Strand mismatch gene "+gene.getStableID()
						+"("+(transcript.getGene().isForward()?"fwd":"rev")
						+") with transcript "+ transcript.getStableID()
						+"("+(rs.getString(4).trim().equals("1")?"fwd":"rev"));
				
				tempTranscriptMap.put(rs.getString(1), transcript);
				transVec.add(transcript);
			}
			
			// do not close - connection is needed for translations and exons
		} catch (SQLException e) { 			// thrown by both methods
			System.err.println("Error fetching transcripts: "+ e);
		}
		
		return Transcript.toTranscriptArray(transVec);
	}
	


	/**
	 * get Genes in single: whole bulk suddenly did not run anymore!
	 */
	Gene[] retrieveGenes_single(Connection con, Gene[] genes) {
		
		try {
			checkVersion(Integer.parseInt(
				con.getCatalog().substring(
						con.getCatalog().indexOf("_core")+ 6,
						con.getCatalog().lastIndexOf('_'))));	// e.g. "homo_sapiens_core_31_35d" = Ensembl ver 31, based on NCBI genome 35
			String s= con.getCatalog().substring(
					con.getCatalog().lastIndexOf('_')+ 1,
					con.getCatalog().length()); // e.g. "homo_sapiens_core_31_35d"
			if (Character.isLetter(s.charAt(s.length()- 1)))
				s= s.substring(0, s.length()- 1); // quit last letter if present (e.g., the "d")
			genes[0].getSpecies().setBuildVersion(
				Integer.parseInt(s));
					
		} catch (SQLException e) {
			e.printStackTrace();
		}

			// build query
		for (int i= 1; i < genes.length; i++) {	// 1000
			StringBuffer sb= new StringBuffer("SELECT gene_stable_id.stable_id,gene.seq_region_start,gene.seq_region_end,gene.seq_region_strand,seq_region.name,gene.biotype,gene.status ");		// 
			sb.append("FROM gene,gene_stable_id,seq_region "); 
			sb.append("WHERE gene.gene_id=gene_stable_id.gene_id");
			sb.append(" AND gene.seq_region_id=seq_region.seq_region_id"); 
	
			sb.append(" AND gene_stable_id.stable_id=\""); 
			sb.append(genes[i].getStableID()); 
			sb.append("\"");
		
				// execute query
			ResultSet rs= null;
			try {
				Statement stmt = con.createStatement(
							ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
												// only one ResultSet per Statement possible
			
				rs = stmt.executeQuery(sb.toString());	// TYPE and CONCUR value given by Statement
				if (rs== null) {
					System.err.println("No result for Genes!");
					return genes;
				}
				
					// store result
				while(rs.next()) {
	
						// find gene
					int x= Arrays.binarySearch(genes, new Gene(genes[0].getSpecies(), rs.getString(1)), new Gene.StableIDComparator());	
					if (x< 0) {
						System.err.println("Gene not found: "+rs.getString(1));
						continue;
					}
					
					genes[x].setStart(Integer.parseInt(rs.getString(2)));
					genes[x].setEnd(Integer.parseInt(rs.getString(3)));
					genes[x].setStrand(Integer.parseInt(rs.getString(4)));
					genes[x].setChromosome(rs.getString(5));	// !! also with NT.. genes (see old method)
//					int chk= genes[x].getConfidence();
					genes[x].setType(rs.getString(6));
//					if (chk!= genes[x].getConfidence())
						//System.err.println("mismatching confidence in gene "+genes[x]);
					//chk= genes[x].getType();
					genes[x].setConfidence(rs.getString(7));
					//if (chk!= genes[x].getType())
						//System.err.println("mismatching type in gene "+genes[x]);
						
				}
				con.close();
			} catch (SQLException e) { 			// thrown by both methods
				System.err.println("Error fetching genes: "+ e);
			}
		}
				
		System.out.println("done.");
		return genes;
	}

	/**
	 * @deprecated: suddenly hanged
	 */
	Gene[] retrieveGenes(Connection con, Gene[] genes) {
		
			// build query
		StringBuffer sb= new StringBuffer("SELECT gene_stable_id.stable_id,gene.seq_region_start,gene.seq_region_end,gene.seq_region_strand,seq_region.name,gene.biotype,gene.status ");		// 
		sb.append("FROM gene,gene_stable_id,seq_region "); 
		sb.append("WHERE gene.gene_id=gene_stable_id.gene_id");
		sb.append(" AND gene.seq_region_id=seq_region.seq_region_id"); 

		sb.append(" AND (gene_stable_id.stable_id=\""); 
		sb.append(genes[0].getStableID()); 
		sb.append("\"");
		for (int i= 1; i < genes.length; i++) {	// 1000
			sb.append(" OR gene_stable_id.stable_id=\"");
			sb.append(genes[i].getStableID());
			sb.append("\"");
		}
		sb.append(")");
		
			// execute query
		try {
			checkVersion(Integer.parseInt(
				con.getCatalog().substring(
						con.getCatalog().indexOf("_core")+ 6,
						con.getCatalog().lastIndexOf('_'))));	// e.g. "homo_sapiens_core_31_35d" = Ensembl ver 31, based on NCBI genome 35
			String s= con.getCatalog().substring(
					con.getCatalog().lastIndexOf('_')+ 1,
					con.getCatalog().length()); // e.g. "homo_sapiens_core_31_35d"
			if (Character.isLetter(s.charAt(s.length()- 1)))
				s= s.substring(0, s.length()- 1); // quit last letter if present (e.g., the "d")
			genes[0].getSpecies().setBuildVersion(
				Integer.parseInt(s));
					
		} catch (SQLException e) {
			e.printStackTrace();
		}
		ResultSet rs= null;
		try {
			Statement stmt = con.createStatement(
						ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
											// only one ResultSet per Statement possible
		
			rs = stmt.executeQuery(sb.toString());	// TYPE and CONCUR value given by Statement
			if (rs== null) {
				System.err.println("No result for Genes!");
				return genes;
			}
			
				// store result
			while(rs.next()) {

					// find gene
				int x= Arrays.binarySearch(genes, new Gene(genes[0].getSpecies(), rs.getString(1)), new Gene.StableIDComparator());	
				if (x< 0) {
					System.err.println("Gene not found: "+rs.getString(1));
					continue;
				}
				
				genes[x].setStrand(Integer.parseInt(rs.getString(4)));
				genes[x].setStart(Integer.parseInt(rs.getString(2)));
				genes[x].setEnd(Integer.parseInt(rs.getString(3)));
				genes[x].setChromosome(rs.getString(5));	// !! also with NT.. genes (see old method)
				genes[x].setType(rs.getString(6));
				genes[x].setConfidence(rs.getString(7));
			}
			con.close();
		} catch (SQLException e) { 			// thrown by both methods
			System.err.println("Error fetching genes: "+ e);
		}
				
		System.out.println("done.");
		return genes;
	}

	/**
	 * @deprecated not in use
	 */
	Vector retrieveGenes(String species, Vector geneVec) {
				
		ResultSet rs= null;
		System.out.println("getting genes for "+ species+" ("+geneVec.size()+" genes)..");
	
			// build query
		StringBuffer sb= new StringBuffer("SELECT gene_stable_id.stable_id,gene.seq_region_start,gene.seq_region_end,gene.seq_region_strand,seq_region.name ");
		sb.append("FROM gene,gene_stable_id,seq_region "); 
		sb.append("WHERE gene.gene_id=gene_stable_id.gene_id");
		sb.append(" AND gene.seq_region_id=seq_region.seq_region_id"); 
		Gene gene= (Gene) geneVec.elementAt(0);
		sb.append(" AND (gene_stable_id.stable_id=\""); 
		sb.append(gene.getStableID()); 
		sb.append("\"");
		for (int i= 1; i < geneVec.size(); i++) {
			gene= (Gene) geneVec.elementAt(i);
			sb.append(" OR gene_stable_id.stable_id=\"");
			sb.append(gene.getStableID());
			sb.append("\"");
		}
		sb.append(")");
		
			// execute query
		Vector removeN2= new Vector();
		Connection con= connect(species+ "_core");	// reconnect to species-db
		try {
			checkVersion(Integer.parseInt(
				con.getCatalog().substring(
						con.getCatalog().indexOf("_core")+ 6,
						con.getCatalog().lastIndexOf('_'))));	// e.g. "homo_sapiens_core_31_35d"
			String s= con.getCatalog().substring(
					con.getCatalog().lastIndexOf('_')+ 1,
					con.getCatalog().length()); // e.g. "homo_sapiens_core_31_35d"
			if (Character.isLetter(s.charAt(s.length()- 1)))
				s= s.substring(0, s.length()- 1); // quit last letter if present (e.g., the "d")
			graph.getSpeciesByName(species).setBuildVersion(
				Integer.parseInt(s));
					
		} catch (SQLException e) {
			e.printStackTrace();
		}
		try {
			Statement stmt = con.createStatement(
						ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
											// only one ResultSet per Statement possible
		
			rs = stmt.executeQuery(sb.toString());	// TYPE and CONCUR value given by Statement
			if (rs== null) {
				System.err.println("No result for Genes!");
				return removeN2;
			}
			
				// store result
			while(rs.next()) {
				gene= graph.getGene(rs.getString(1));	// stableID
				if (gene== null)
					continue;
				
				String s= rs.getString(5).trim();	// filter off all NT... contigs
				try {
					Integer.parseInt(s);
				} catch (NumberFormatException e) {
					if (s!= "X"&& s!= "Y") {
						Gene[] removed= graph.removePurge(gene);		// eliminate NT..
						for (int i = 0; i < removed.length; i++) 
							if (!geneVec.remove(removed[i]))
								removeN2.add(removed[i]);	// remove in gene-Vector of the other species
						continue;
					}
				}
				
				gene.setStart(Integer.parseInt(rs.getString(2)));
				gene.setEnd(Integer.parseInt(rs.getString(3)));
				gene.setForward(rs.getString(4));
				gene.setChromosome(rs.getString(5));
			}
			con.close();
		} catch (SQLException e) { 			// thrown by both methods
			System.err.println("Error fetching genes: "+ e);
		}
				
		System.out.println("done.");
		return removeN2;
	}
	
		
		

	/**
	 * @deprecated not in use
	 */
	void retrieveExons(Connection con, Transcript[] transcr) {
		
		ResultSet rs= null;
		System.out.println("["+ new Date(System.currentTimeMillis())+"]: getting exons ("+transcr.length+" transcripts)");

			// build query
		StringBuffer sb= new StringBuffer(
			"SELECT exon_stable_id.stable_id,exon.seq_region_start,exon.seq_region_end,phase,end_phase,transcript_stable_id.stable_id ");
		sb.append("FROM exon,exon_stable_id,exon_transcript,transcript,transcript_stable_id "); 
		sb.append("WHERE exon_stable_id.exon_id=exon.exon_id");
		sb.append(" AND transcript.transcript_id=exon_transcript.transcript_id");
		sb.append(" AND exon.exon_id=exon_transcript.exon_id");
		sb.append(" AND transcript.transcript_id=transcript_stable_id.transcript_id");
		Transcript trans= (Transcript) transcr[0];
		sb.append(" AND (transcript_stable_id.stable_id=\""); 
		sb.append(trans.getStableID()); 
		sb.append("\"");
		for (int i= 1; i < transcr.length; i++) {
			trans= (Transcript) transcr[i];
			sb.append(" OR transcript_stable_id.stable_id=\"");
			sb.append(trans.getStableID());
			sb.append("\"");
		}
		sb.append(")");

			// execute query
		try {
			Statement stmt = con.createStatement(
						ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
											// only one ResultSet per Statement possible

			rs = stmt.executeQuery(sb.toString());	// TYPE and CONCUR value given by Statement
			System.out.println("-query executed.."+ new Date(System.currentTimeMillis()));
			if (rs== null) {
				System.err.println("No result for Exons!");
				return;
			}
	
				// store result			
			if (tempExonMap== null) {	// get or create transcript map
				int ctr= 0;
				while(rs.next()) ++ctr;
				tempExonMap= new HashMap(ctr, (int) (1.2* ctr));
				rs.beforeFirst();
			}
			while(rs.next()) {
				
				trans= (Transcript) tempTranscriptMap.get(rs.getString(6));	// stableID
				int start= Integer.parseInt(rs.getString(2));
				int end= Integer.parseInt(rs.getString(3));
				Exon e= trans.getGene().getExon(start, end); // exon exists?
				
				if (e== null) {
					e= new Exon(trans, rs.getString(1), start, end);	// stableID
					// e.setGene(trans.getGene());
					e.setFrame(trans.getDefaultTranslation(), new Phase(
							Byte.parseByte(rs.getString(4)), 
							Byte.parseByte(rs.getString(5))));
				}
				e.addTranscript(trans);	
				trans.addExon(e);
				
//				if (!e.checkStrand(rs.getString(6)))
//					System.err.println("Strand mismatch of transcript "+trans.getStableID()
//						+" with exon "+ e.getStableID());
				
				tempExonMap.put(rs.getString(1), e);
			}
			// connection closed in getTranscripts()
		} catch (SQLException e) { 			// thrown by both methods
			System.err.println(
				"Error fetching exons "+ e);
		}
	}	
	
	final static int countResults(ResultSet rs) {
		
		int ctr= 0;
		try {
			int row= rs.getRow();
			rs.beforeFirst();
			while (rs.next())
				ctr++;
			rs.beforeFirst();
			while(row-->0 )
				rs.next();

			rs.beforeFirst();	// reset ResultSet
		} catch (SQLException e) {
			e.printStackTrace();
		}

		return ctr;		
	}
	
	/**
	 * Retrieves all Genes that have a homolog with the given <code>refSpecies</code>.  
	 * 
	 * Uses <code>ensembl_mart</code> database, but could be extended by usage of 
	 * <code>ensembl_compara</code>.  
	 *  
	 * @param g <code>Graph</code> that provides the <code>Species</code> for which homolog
	 * genes have to be retrieved.
	 * @return <code>Graph</code> filled with <code>Gene</code> instances:
	 */
	GeneHomology[] retrieveHomologGenesReference_mart(Graph g, String refSpecName) {
		
		Species[] species= g.getSpecies();
		Vector allHomologies= new Vector(10000);	// bit bigger than default (16)
		
		// connect
		Connection con= connect("ensembl_mart");	// ensembl_mart_29
		try {
			checkVersion(Integer.parseInt(
				con.getCatalog().substring(con.getCatalog().lastIndexOf('_')+ 1)));	// "ensembl_mart_31"
		} catch (SQLException e) {
			System.err.println(e);;
		}

		Species refSpec= g.getSpeciesByName(refSpecName);
		for (int i= 0; i < species.length; i++) {
			
			if (species[i]== refSpec)
				continue;
			System.out.println(species[i]+" x "+ refSpec);
				// query
			ResultSet rs= null;
			String tableGmainI= species[i].getNameAbbrev()+"_gene_ensembl__gene__main";	
			String tableGmainJ= refSpec.getNameAbbrev()+"_gene_ensembl__gene__main";	
			String tableXlink= species[i].getNameAbbrev()+"_gene_ensembl__homologs_"+refSpec.getNameAbbrev()+"__dm";	
			try {
				Statement stmt = con.createStatement(
							ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY); // only one ResultSet per Statement possible
												
				String query= "SELECT "+tableXlink+".gene_stable_id,"+tableXlink+".homol_stable_id,"
								+ tableGmainI+".biotype,"+ tableGmainI+".confidence,"
								+ tableGmainI+".chr_name,"+ tableGmainI+".gene_chrom_start,"+ tableGmainI+".gene_chrom_end,"
								+ tableGmainJ+".biotype,"+ tableGmainJ+".confidence," 
								+ tableGmainJ+".chr_name,"+ tableGmainJ+".gene_chrom_start,"+ tableGmainJ+".gene_chrom_end"
								+ " FROM "+ tableGmainI+","+tableGmainJ+","+tableXlink
								+" WHERE homol_id is not NULL"	// does occur !!
								// +" AND description=\"UBRH\""	// we want all
								+" AND "+tableGmainI+".gene_stable_id="+tableXlink+".gene_stable_id"	// link tables
								+" AND "+tableGmainJ+".gene_stable_id="+tableXlink+".homol_stable_id";														
				
				rs = stmt.executeQuery(query);	// TYPE and CONCUR value given by Statement
				 
			} catch (SQLException e) { 			// thrown by both methods
				System.err.println(e);
			}
			
			if (rs== null) {		// table does not exist
				System.err.println("Error when trying to retrieve "+species[i]+"x"+refSpec+" homologs.");
				continue;
			}

			
				// iterate result set, Xlink genes
			Vector newGenes1= new Vector();
			Vector newGenes2= new Vector();	
			try {	
					// prepare species
				int ctr= countResults(rs);	// estimate gene nb (!! gene duplications)
				rs.first();	// jump
				String stableID1= rs.getString(1);
				String stableID2= rs.getString(2);
				Species spec1= g.getSpeciesByGeneID(stableID1);
				if (spec1== null)
					spec1= g.getSpeciesByEnsemblPrefix(Species.SP_NAMES_ENS_PFX[11]);	// guess tetraodon
				if (spec1.getGeneNb()< 1)
					spec1.setEstimatedGeneNb(ctr/2);
				Species spec2= g.getSpeciesByGeneID(stableID2);	// generate Species
				if (spec2== null)
					spec2= g.getSpeciesByEnsemblPrefix(Species.SP_NAMES_ENS_PFX[11]);	// guess tetraodon
				if (spec2.getGeneNb()< 1) 						
					spec2.setEstimatedGeneNb(ctr/2);

					// get genes
				rs.beforeFirst();
				while (rs.next()) {				// get body of table
					stableID1= rs.getString(1);	// renew
					Gene gene1= g.getGene(stableID1);
					if (gene1== null) {
						newGenes1.add(gene1= new Gene(spec1, stableID1));
						gene1.setSpecies(spec1);
						gene1.setType(rs.getString(3));
						gene1.setConfidence(rs.getString(4));
						gene1.setChromosome(rs.getString(5));
						gene1.setStart(Integer.parseInt(rs.getString(6)));
						gene1.setEnd(Integer.parseInt(rs.getString(7)));
						g.addGene(gene1);
					}
					stableID2= rs.getString(2);	// renew
					Gene gene2= g.getGene(stableID2);
					if (gene2== null) { 
						newGenes2.add(gene2= new Gene(spec2, stableID2));
						gene2.setSpecies(spec2);
						gene2.setType(rs.getString(8));
						gene2.setConfidence(rs.getString(9));
						gene2.setChromosome(rs.getString(10));
						gene2.setStart(Integer.parseInt(rs.getString(11)));
						gene2.setEnd(Integer.parseInt(rs.getString(12)));
						g.addGene(gene2);
					}
					
					GeneHomology hom= new GeneHomology(gene1, gene2);
					//hom.setDs();
					//hom.setDn();
					//hom.setDnDs();
					gene1.addHomology(hom);	// Xlink
					gene2.addHomology(hom);
					allHomologies.add(hom);
				}
				con.close();
			} catch (SQLException e) {
				e.printStackTrace();
			}				
			
		}

		try {con.close();
		BaseDriver.close(con);} 
		catch (SQLException e) {e.printStackTrace();}				
		return GeneHomology.toGeneHomologyArray(allHomologies); 
	}

	/**
	 * Retrieves all Genes that have at least one homolog in ALL <code>species</code>.
	 * (!!! check whether ALL makes sense !!!)  
	 * 
	 * Uses <code>ensembl_mart</code> database, but could be extended by usage of 
	 * <code>ensembl_compara</code>.  
	 *  
	 * @param g <code>Graph</code> that provides the <code>Species</code> for which homolog
	 * genes have to be retrieved.
	 * @return <code>Graph</code> filled with <code>Gene</code> instances:
	 */
	GeneHomology[] retrieveHomologGenesIterative_mart(Graph g, String refSpecName) {
		
		Species[] species= g.getSpecies();
		Vector allHomologies= new Vector(10000);	// bit bigger than default (16)
		
		// connect
		Connection con= connect("ensembl_mart");	// ensembl_mart_29
		try {
			checkVersion(Integer.parseInt(
				con.getCatalog().substring(con.getCatalog().lastIndexOf('_')+ 1)));	// "ensembl_mart_31"
		} catch (SQLException e) {
			System.err.println(e);;
		}
	
		Species refSpec= g.getSpeciesByName(refSpecName);
		for (int i= 0; i < species.length; i++) {
			if (refSpec== species[i]) continue;
			for (int j= (i+1); j < species.length; j++) {				
				if (refSpec== species[j]) continue;
				System.out.println(species[i]+" x "+ species[j]);
					// query
				ResultSet rs= null;
				String tableGmainI= species[i].getNameAbbrev()+"_gene_ensembl__gene__main";	
				String tableGmainJ= species[j].getNameAbbrev()+"_gene_ensembl__gene__main";	
				String tableXlink= species[i].getNameAbbrev()+"_gene_ensembl__homologs_"+species[j].getNameAbbrev()+"__dm";
				String query;
				try {
					Statement stmt = con.createStatement(
								ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY); // only one ResultSet per Statement possible
													
					Gene[] g1= species[i].getGenes();
					Gene[] g2= species[j].getGenes();
					query= "SELECT DISTINCT "+tableXlink+".gene_stable_id,"+tableXlink+".homol_stable_id,"
									+ tableGmainI+".biotype,"+ tableGmainI+".confidence,"
									+ tableGmainI+".chr_name,"+ tableGmainI+".gene_chrom_start,"+ tableGmainI+".gene_chrom_end,"
									+ tableGmainJ+".biotype,"+ tableGmainJ+".confidence," 
									+ tableGmainJ+".chr_name,"+ tableGmainJ+".gene_chrom_start,"+ tableGmainJ+".gene_chrom_end"
									+ " FROM "+ tableGmainI+","+tableGmainJ+","+tableXlink
									+" WHERE homol_id is not NULL"	// does occur !!
									// +" AND description=\"UBRH\""	// we want all
									+" AND "+tableGmainI+".gene_stable_id="+tableXlink+".gene_stable_id"	// link tables
									+" AND "+tableGmainJ+".gene_stable_id="+tableXlink+".homol_stable_id";
									
									query+=" AND (";		// filter for existing genes
									for (int k = 0; k < g1.length; k++) 
										query+= tableGmainI+".gene_stable_id=\""+ g1[k].getStableID()+ "\" OR ";
									for (int k = 0; k < g2.length; k++) {
										query+= tableGmainJ+".gene_stable_id=\""+ g2[k].getStableID()+ "\"";
										if (k+ 1< g2.length)
											query+= " OR ";
									}
									query+= ")";
					
					rs = stmt.executeQuery(query);	// TYPE and CONCUR value given by Statement
					 
				} catch (SQLException e) { 			// thrown by both methods
					System.err.println(e);
				}

				if (rs== null) {		// table does not exist
					System.err.println("Error when trying to retrieve "+species[i]+"x"+species[j]+" homologs.");
					continue;
				}
	
				
					// iterate result set, Xlink genes
				Vector newGenes1= new Vector();
				Vector newGenes2= new Vector();	
				try {	
						// prepare species
					int ctr= countResults(rs);	// estimate gene nb (!! gene duplications)
					rs.first();	// jump
					String stableID1= rs.getString(1);
					String stableID2= rs.getString(2);
					Species spec1= g.getSpeciesByGeneID(stableID1);
					if (spec1== null)
						spec1= g.getSpeciesByEnsemblPrefix(Species.SP_NAMES_ENS_PFX[11]);	// guess tetraodon
					if (spec1.getGeneNb()< 1)
						spec1.setEstimatedGeneNb(ctr/2);
					Species spec2= g.getSpeciesByGeneID(stableID2);	// generate Species
					if (spec2== null)
						spec2= g.getSpeciesByEnsemblPrefix(Species.SP_NAMES_ENS_PFX[11]);	// guess tetraodon
					if (spec2.getGeneNb()< 1) 						
						spec2.setEstimatedGeneNb(ctr/2);
	
						// get genes
					rs.beforeFirst();
					while (rs.next()) {				// get body of table
						stableID1= rs.getString(1);	// renew
						Gene gene1= g.getGene(stableID1);
						if (gene1== null) {
							newGenes1.add(gene1= new Gene(spec1, stableID1));
							gene1.setSpecies(spec1);
							gene1.setType(rs.getString(3));
							gene1.setConfidence(rs.getString(4));
							gene1.setChromosome(rs.getString(5));
							gene1.setStart(Integer.parseInt(rs.getString(6)));
							gene1.setEnd(Integer.parseInt(rs.getString(7)));
							g.addGene(gene1);
						}
						stableID2= rs.getString(2);	// renew
						Gene gene2= g.getGene(stableID2);
						if (gene2== null) { 
							newGenes2.add(gene2= new Gene(spec2, stableID2));
							gene2.setSpecies(spec2);
							gene2.setType(rs.getString(8));
							gene2.setConfidence(rs.getString(9));
							gene2.setChromosome(rs.getString(10));
							gene2.setStart(Integer.parseInt(rs.getString(11)));
							gene2.setEnd(Integer.parseInt(rs.getString(12)));
							g.addGene(gene2);
						}
						
						GeneHomology hom= new GeneHomology(gene1, gene2);
						//hom.setDs();
						//hom.setDn();
						//hom.setDnDs();
						gene1.addHomology(hom);	// Xlink
						gene2.addHomology(hom);
						allHomologies.add(hom);
					}
					con.close();
				} catch (SQLException e) {
					e.printStackTrace();
				}				
				
			}
		}
	
		try {con.close();
		BaseDriver.close(con);} 
		catch (SQLException e) {e.printStackTrace();}				
		return GeneHomology.toGeneHomologyArray(allHomologies); 
	}

	/**
	 * Retrieves all Genes that have at least one homolog in ALL <code>species</code>.
	 * (!!! check whether ALL makes sense !!!)  
	 * 
	 * Uses <code>ensembl_mart</code> database, but could be extended by usage of 
	 * <code>ensembl_compara</code>.  
	 *  
	 * @param g <code>Graph</code> that provides the <code>Species</code> for which homolog
	 * genes have to be retrieved.
	 * @return <code>Graph</code> filled with <code>Gene</code> instances:
	 */
	GeneHomology[] retrieveHomologGenesAll_mart(Graph g) {
		
		Species[] species= g.getSpecies();
		Vector allHomologies= new Vector(10000);	// bit bigger than default (16)
		
		// connect
		Connection con= connect("ensembl_mart");	// ensembl_mart_29
		try {
			checkVersion(Integer.parseInt(
				con.getCatalog().substring(con.getCatalog().lastIndexOf('_')+ 1)));	// "ensembl_mart_31"
		} catch (SQLException e) {
			System.err.println(e);;
		}
	
		for (int i= 0; i < species.length; i++) {
			for (int j= (i+1); j < species.length; j++) {
				
				System.out.println(species[i]+" x "+ species[j]);
					// query
				ResultSet rs= null;
				String tableGmainI= species[i].getNameAbbrev()+"_gene_ensembl__gene__main";	
				String tableGmainJ= species[j].getNameAbbrev()+"_gene_ensembl__gene__main";	
				String tableXlink= species[i].getNameAbbrev()+"_gene_ensembl__homologs_"+species[j].getNameAbbrev()+"__dm";	
				try {
					Statement stmt = con.createStatement(
								ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY); // only one ResultSet per Statement possible
													
					String query= "SELECT "+tableXlink+".gene_stable_id,"+tableXlink+".homol_stable_id,"
									+ tableGmainI+".biotype,"+ tableGmainI+".confidence,"
									+ tableGmainI+".chr_name,"+ tableGmainI+".chrom_strand," 
									+ tableGmainI+".gene_chrom_start,"+ tableGmainI+".gene_chrom_end,"
									+ tableGmainJ+".biotype,"+ tableGmainJ+".confidence," 
									+ tableGmainJ+".chr_name,"+ tableGmainJ+".chrom_strand,"
									+ tableGmainJ+".gene_chrom_start,"+ tableGmainJ+".gene_chrom_end, "
									+ tableXlink+".description"
									+ " FROM "+ tableGmainI+","+tableGmainJ+","+tableXlink
									+" WHERE homol_id is not NULL"	// does occur !!
									// +" AND description=\"UBRH\""	// we want all
									+" AND "+tableGmainI+".gene_stable_id="+tableXlink+".gene_stable_id"	// link tables
									+" AND "+tableGmainJ+".gene_stable_id="+tableXlink+".homol_stable_id";														
					
					for (int k = 0; k < species.length; k++) {	// ensure there are homologs in all other species
						if (k== i || k== j)
							continue;
						query+= " AND "+tableGmainI+"."+species[k].getNameAbbrev()+"_homolog_bool"
								+ " AND "+tableGmainJ+"."+species[k].getNameAbbrev()+"_homolog_bool";
					}
					query+= " AND "+tableXlink+".description='ortholog_one2one'";
					
					rs = stmt.executeQuery(query);	// TYPE and CONCUR value given by Statement
					 
				} catch (SQLException e) { 			// thrown by both methods
					System.err.println(e);
				}
				
				if (rs== null) {		// table does not exist
					System.err.println("Error when trying to retrieve "+species[i]+"x"+species[j]+" homologs.");
					continue;
				}
	
				
					// iterate result set, Xlink genes
				Vector newGenes1= new Vector();
				Vector newGenes2= new Vector();	
				try {	
						// prepare species
					int ctr= countResults(rs);	// estimate gene nb (!! gene duplications)
					rs.first();	// jump
					String stableID1= rs.getString(1);
					String stableID2= rs.getString(2);
					if (species[i].getGeneNb()< 1)
						species[i].setEstimatedGeneNb(ctr/2);
					if (species[j].getGeneNb()< 1) 						
						species[j].setEstimatedGeneNb(ctr/2);
	
						// get genes
					rs.beforeFirst();
					int x= 0;
					while (rs.next()) {				// get body of table
						if (DEBUG&& ++x>= 100)
							break;
						stableID1= rs.getString(1);	// renew
						Gene gene1= species[i].getGene(stableID1);
						if (gene1== null) {
							newGenes1.add(gene1= new Gene(species[i], stableID1));
							gene1.setSpecies(species[i]);
							gene1.setType(rs.getString(3));
							gene1.setConfidence(rs.getString(4));
							gene1.setChromosome(rs.getString(5));
							gene1.setStrand(Integer.parseInt(rs.getString(6)));
							gene1.setStart(Integer.parseInt(rs.getString(7)));
							gene1.setEnd(Integer.parseInt(rs.getString(8)));
							g.addGene(gene1);
						}
						stableID2= rs.getString(2);	// renew
						Gene gene2= species[j].getGene(stableID2);
						if (gene2== null) { 
							newGenes2.add(gene2= new Gene(species[j], stableID2));
							gene2.setSpecies(species[j]);
							gene2.setType(rs.getString(9));
							gene2.setConfidence(rs.getString(10));
							gene2.setChromosome(rs.getString(11));
							gene2.setStrand(Integer.parseInt(rs.getString(12)));
							gene2.setStart(Integer.parseInt(rs.getString(13)));
							gene2.setEnd(Integer.parseInt(rs.getString(14)));
							g.addGene(gene2);
						}
						
						GeneHomology hom= new GeneHomology(gene1, gene2);
						hom.setType(rs.getString(15));
						//hom.setDs();
						//hom.setDn();
						//hom.setDnDs();
						gene1.addHomology(hom);	// Xlink
						gene2.addHomology(hom);
						allHomologies.add(hom);
					}
					con.close();
				} catch (SQLException e) {
					e.printStackTrace();
				}				
				
			}
		}
	
		try {con.close();
		BaseDriver.close(con);} 
		catch (SQLException e) {e.printStackTrace();}				
		return GeneHomology.toGeneHomologyArray(allHomologies); 
	}

	Graph retrieveGenesAll_mart(Graph g, Species spec) {
		
		Species[] species= g.getSpecies();
		
		// connect
		Connection con= connect("ensembl_mart");	// ensembl_mart_29
		try {
			checkVersion(Integer.parseInt(
				con.getCatalog().substring(con.getCatalog().lastIndexOf('_')+ 1)));	// "ensembl_mart_31"
		} catch (SQLException e) {
			System.err.println(e);;
		}
	
		for (int i= 0; i < species.length; i++) {
				System.out.println(species[i]);
					// query
				ResultSet rs= null;
				String tableGmainI= species[i].getNameAbbrev()+"_gene_ensembl__gene__main";	
				try {
					Statement stmt = con.createStatement(
								ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY); // only one ResultSet per Statement possible
													
					String query= "SELECT "+tableGmainI+".gene_stable_id,"
									+ tableGmainI+".biotype,"+ tableGmainI+".confidence,"
									+ tableGmainI+".chr_name,"+ tableGmainI+ ".chrom_strand,"
									+ tableGmainI+".gene_chrom_start,"+ tableGmainI+".gene_chrom_end"
									+ " FROM "+ tableGmainI;
									//+ "  WHERE "+ tableGmainI+".biotype=\"protein_coding\"";
					
					rs = stmt.executeQuery(query);	// TYPE and CONCUR value given by Statement
					 
				} catch (SQLException e) { 			// thrown by both methods
					System.err.println(e);
				}
				
				if (rs== null) {		// table does not exist
					System.err.println("Error when trying to retrieve "+species[i]);
					continue;
				}
	
				
					// iterate result set, Xlink genes
				Vector newGenes1= new Vector();
				try {	
						// prepare species
					int ctr= countResults(rs);	// estimate gene nb (!! gene duplications)
					rs.first();	// jump
					String stableID1= rs.getString(1);
					if (spec.getGeneNb()< 1)
						spec.setEstimatedGeneNb(ctr/2);
	
						// get genes
					rs.beforeFirst();
					int x= 0;
					while (rs.next()) {				// get body of table
						if (DEBUG&& ++x>= 100)
							break;
						stableID1= rs.getString(1);	// renew
						Gene gene1= species[i].getGene(stableID1); // was g.getGene(stableID1);
						if (gene1== null) {
							newGenes1.add(gene1= new Gene(spec, stableID1)); 
							gene1.setSpecies(spec);
							gene1.setType(rs.getString(2));
							gene1.setConfidence(rs.getString(3));
							gene1.setChromosome(rs.getString(4));
							gene1.setStrand(Integer.parseInt(rs.getString(5)));
							gene1.setStart(Integer.parseInt(rs.getString(6)));
							gene1.setEnd(Integer.parseInt(rs.getString(7)));
							g.addGene(gene1);
						}
					}
					con.close();
				} catch (SQLException e) {
					e.printStackTrace();
				}				
				
		}
	
		try {con.close();
		BaseDriver.close(con);} 
		catch (SQLException e) {e.printStackTrace();}				
		return g; 
	}
	
	
	/**
	 * 
	 */
	GeneHomology[] retrieveHomologyInfo_single(GeneHomology[] homols) {
		
		

			// connect
		Connection con= connect("ensembl_compara");	// ensembl_compara_35
		try {
			checkVersion(Integer.parseInt(
				con.getCatalog().substring(con.getCatalog().lastIndexOf('_')+ 1)));	// "ensembl_mart_31"
		} catch (SQLException e) {
			System.err.println(e);;
		}
		System.out.println("0    10   20   30   40   50   60   70   80   90  100[%]");
		System.out.println("|----|----|----|----|----|----|----|----|----|----|");
		System.out.print("*");
		System.out.flush();
		int pctCtr= 0;
		int pctBase= homols.length/ 100;
		if (pctBase< 1)	// for debugging, little gene sets
			pctBase= 1;

		for (int i= 0; i < homols.length; i++) {

				// query
			ResultSet rs= null;
			try {
				Statement stmt = con.createStatement(
							ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY); // only one ResultSet per Statement possible
												
				String query= "SELECT DISTINCT "	// see method_link_species_set-Pfusch
								//+ "homology.homology_id, homology.description, "
								+ "homology.homology_id, "
								+ "homology.description, type, n, s, dn, ds, threshold_on_ds, lnl," 
								+ "perc_cov, perc_id, perc_pos,"// homology_member 
								+ "member.stable_id"			// unused just for strategy works with distinct: 
																// s. method_link_species_set-Pfusch
								+ " FROM homology, homology_member, member, method_link_species_set, method_link"
								+" WHERE (member.stable_id=\""+ homols[i].getGene1().getStableID()+"\""
								+" OR member.stable_id=\""+ homols[i].getGene2().getStableID()+"\")"
								+" AND homology_member.member_id=member.member_id"
								+" AND homology_member.homology_id=homology.homology_id"
								+" AND homology.method_link_species_set_id=method_link_species_set.method_link_species_set_id"
								+" AND method_link_species_set.method_link_id=method_link.method_link_id";
								
				rs = stmt.executeQuery(query);	// rs contains homologies between given 2 genes and all species in ensembl
												// however, the common relation appears twice (homology_id)
				 
			} catch (SQLException e) { 			// thrown by both methods
				System.err.println(e);
			}
			
			if (rs== null) {		// table does not exist
				System.err.println("Error when trying to retrieve "+homols[i]+" homology-info.");
				continue;
			}

			
			try {	
					// rs contains homologies for the two given genes, 
					// find the two that are linked by a common homol_id
				int[] homolID= new int[countResults(rs)];
				rs.first();
				for (int j = 0; j < homolID.length; j++, rs.next())		// fill in IDs 
					homolID[j]= Integer.parseInt(rs.getString(1));
				int x1= -1, x2= -1;										// find identical IDs
				for (x1 = 0; x1 < homolID.length; x1++) {
					for (x2 = (x1+1); x2 < homolID.length; x2++) 
						if (homolID[x1]== homolID[x2])
							break;
					if (x2< homolID.length)
						break;
				}
				
				rs.first();
				for (int j = 0; j < x1; j++) 
					rs.next();	// jump to homolog
				
				
				homols[i].setType(rs.getString(2));
				//homols[i].setSubtype(rs.getString(3));	// gibts net mehr
				homols[i].setMethod(rs.getString(3));
				homols[i].setN(rs.getDouble(4));
				homols[i].setS(rs.getDouble(5));
				homols[i].setDn(rs.getDouble(6));
				homols[i].setDs(rs.getDouble(7));
				homols[i].setThresholdOnDS(rs.getDouble(8));
				homols[i].setLnl(rs.getDouble(9));
				homols[i].setG1PercCov(rs.getInt(10));
				homols[i].setG1PercId(rs.getInt(11));
				homols[i].setG1PercPos(rs.getInt(12));
				homols[i].setG2PercCov(rs.getInt(10));
				homols[i].setG2PercId(rs.getInt(11));
				homols[i].setG2PercPos(rs.getInt(12));
				
				
				for (int j = 0; j < (x2- x1); j++) 
					rs.next();	// jump
				
			} catch (SQLException e) {
				e.printStackTrace();
			}				
			
			if (i/ pctBase> pctCtr) {
				pctCtr++;
				if (pctCtr% 2== 0) {
					System.out.print("*");
					System.out.flush();
				}
			}
		}

		System.out.println();
		return homols; 
	}

	/**
	 * @deprecated db connection times out!
	 * 
	 * @param homols
	 * @return
	 */
	GeneHomology[] retrieveHomologyInfo(GeneHomology[] homols) {
		
		
	
			// connect
		Connection con= connect("ensembl_compara");	// ensembl_compara_35
		try {
			checkVersion(Integer.parseInt(
				con.getCatalog().substring(con.getCatalog().lastIndexOf('_')+ 1)));	// "ensembl_compara_31"
		} catch (SQLException e) {
			System.err.println(e);;
		}

			// query
		ResultSet rs= null;
		try {
			Statement stmt = con.createStatement(
						ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY); // only one ResultSet per Statement possible
											
			String query= "SELECT DISTINCT "	// see method_link_species_set-Pfusch
							+ "homology.description, subtype, type, n, s, dn, ds, threshold_on_ds, lnl," 
							+ "hm1.perc_cov, hm1.perc_id, hm1.perc_pos,"// homology_member 
							+ "hm2.perc_cov, hm2.perc_id, hm2.perc_pos"// homology_member 
							+ ", m1.stable_id, m2.stable_id"	// test for same alphabet order
							+ " FROM homology FORCE INDEX (PRIMARY), homology_member AS hm1, homology_member AS hm2, "
							+ " member AS m1, member AS m2, method_link_species_set, method_link"
							+" WHERE hm1.homology_id= homology.homology_id"
							+" AND hm2.homology_id= homology.homology_id"
							+" AND hm1.member_id= m1.member_id"
							+" AND hm2.member_id= m2.member_id"
							+" AND homology.method_link_species_set_id=method_link_species_set.method_link_species_set_id"
							+" AND method_link_species_set.method_link_id=method_link.method_link_id";
			
//			query+= " AND ";
//			for (int i= 0; i < homols.length; i++) {
//				query+= "(m1.stable_id=\""+ homols[i].getGene1().getStableID()
//						+"\" AND m2.stable_id=\""+ homols[i].getGene2().getStableID()+"\") OR ";
//			}
//			query= query.substring(1, query.length()- 3);
			
			
				
			rs = stmt.executeQuery(query);	// rs contains homologies between given 2 genes and all species in ensembl
											// however, the common relation appears twice (homology_id)
			 
		} catch (SQLException e) { 			// thrown by both methods
			System.err.println(e);
		}
		
		if (rs== null) {		// table does not exist
			System.err.println("No homology-info.");
		}

		
		try {	

			Arrays.sort(homols, new GeneHomology.StableIDComparator());
			while (rs.next()) {
				
				Gene g1= new Gene(graph.getSpeciesByGeneID(rs.getString(16)), rs.getString(16));
				Gene g2= new Gene(graph.getSpeciesByGeneID(rs.getString(17)), rs.getString(17));
				GeneHomology test= new GeneHomology(g1, g2);
				int i= Arrays.binarySearch(homols, test, new GeneHomology.StableIDComparator());
				
				
				homols[i].setType(rs.getString(1));
				homols[i].setSubtype(rs.getString(2));
				homols[i].setMethod(rs.getString(3));
				homols[i].setN(rs.getDouble(4));
				homols[i].setS(rs.getDouble(5));
				homols[i].setDn(rs.getDouble(6));
				homols[i].setDs(rs.getDouble(7));
				homols[i].setThresholdOnDS(rs.getDouble(8));
				homols[i].setLnl(rs.getDouble(9));
				homols[i].setG1PercCov(rs.getInt(10));
				homols[i].setG1PercId(rs.getInt(11));
				homols[i].setG1PercPos(rs.getInt(12));
				homols[i].setG2PercCov(rs.getInt(13));
				homols[i].setG2PercId(rs.getInt(14));
				homols[i].setG2PercPos(rs.getInt(15));
			}
			
			
		} catch (SQLException e) {
			e.printStackTrace();
		}				
		
		return homols; 
	}	
	
	Graph retrieveGeneInfo(Graph g) {
		
		Species[] spec= g.getSpecies();
		for (int i = 0; i < spec.length; i++) {		// each species is one connect!
								
			String specName= spec[i].getBinomialName();
			Gene[] specGenes= spec[i].getGenes();
			Arrays.sort(specGenes, new Gene.StableIDComparator());	// binary search in retriveGenes()
			
			
				// connect
			Connection con= connect(specName+ "_core");	// reconnect to species-db
			try {
				checkVersion(Integer.parseInt(
					con.getCatalog().substring(
							con.getCatalog().indexOf("_core")+ 6,
							con.getCatalog().lastIndexOf('_'))));	// e.g. "homo_sapiens_core_31_35d"
			} catch (SQLException e) {
				e.printStackTrace();
			}

				// requests
			System.out.println("["+ new Date(System.currentTimeMillis())+"]: getting gene info for "+ 
					specName+"("+specGenes.length+" genes)");
			specGenes= retrieveGenes_single(con, specGenes);	// genes

			System.out.println("["+ new Date(System.currentTimeMillis())+"]: getting transcripts for "+ 
					specName+"("+specGenes.length+" genes)");
			Transcript[] specTrans= retrieveTranscripts(con, spec[i], specGenes);	// transcripts
			
			System.out.println("["+ new Date(System.currentTimeMillis())+ "]: getting translations ("+
					specTrans.length+" transcripts)");
			retrieveTranslations(con, specTrans);								// proteins
			
			System.out.println("getting exons ("+specTrans.length+" transcripts)");		// exons
			System.out.println("0    10   20   30   40   50   60   70   80   90  100[%]");
			System.out.println("|----|----|----|----|----|----|----|----|----|----|");
			System.out.print("*");
			retrieveExonsSingle(con, specTrans);	// strictly after adding transcripts to gene
											// necessary for initing AS types correctly
											// see: Gene.addExon()
			

				// close
			try { con.close(); } 
			catch (SQLException e1) { e1.printStackTrace(); } 
		}
		
		return g;
	}

	/**
	 * @deprecated not in use
	 */
	ResultSet retrieveJoinSpecies(Connection con, String species1, String species2) {
		
		
		ResultSet rs= null;
		String table= species1+"_gene_ensembl__homologs_"+species2+"__dm";	// for human also VEGA!
		try {
			Statement stmt = con.createStatement(
						ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
											// only one ResultSet per Statement possible
			String query= "SELECT gene_stable_id,homol_stable_id FROM "+table
							+" WHERE homol_id is not NULL ";	// does occur !!
							// +" AND description=\"UBRH\"";	// we want all

			rs = stmt.executeQuery(query);	// TYPE and CONCUR value given by Statement
			 
		} catch (SQLException e) { 			// thrown by both methods
			System.err.println(e);
		}
		
		return rs;
	}

	void outputResultSet(ResultSet rs) {

		try {		
			ResultSetMetaData rsmd= rs.getMetaData();
			int cols= rsmd.getColumnCount();
	
			for (int i= 1; i<= cols; i++)	// get header of table
				System.out.print(rsmd.getColumnName(i)+ "\t");
			System.out.println("\n--");
			
			while (rs.next()) {				// get body of table
				for (int i= 1; i<= cols; i++)
					System.out.println(rs.getString(i)+ "\t");
			}		
		} catch (SQLException e) {
			e.printStackTrace();
		}
	}
	
	
	public boolean testGraphEncode(Graph g) {
		
		FASTAWrapper fasta= new FASTAWrapper(
				Constants.HOME_DIR+ File.separator+ GraphHandler.GRAPH_SUBDIR+ File.separator+ Graph.GRAPH_ENCODE_REF_GENES);
		try {
			fasta.read();
		} catch (Exception e) {
			e.printStackTrace();
		}
		String[] names= fasta.getSeqNames();
		String[] seqs= fasta.getSequences();
		
		Species[] spec= g.getSpecies();
		boolean tst= true;
		for (int i = 0; i < spec.length; i++) {
			Gene[] ge= spec[i].getGenes();
			
			String nm= ge[0].getStableID();				// check first gene
			int pos;
			for (pos = 0; pos< names.length; pos++) 
				if (nm.equalsIgnoreCase(names[pos]))
					break;
			if (pos>= names.length) { 
				System.err.println(nm+ " not found in reference data set.");
				tst= false;
			} else {
				try {
					String test= g.readSequence(ge[0]);
					if (!seqs[pos].equalsIgnoreCase(test)) { 
						System.err.println(nm+ "mismatches the reference:\n\t"
								+ seqs[pos]+ "\n\t"+ test);
						tst= false;
					} else
						System.out.println(nm+ " passed test.");
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
			
			
			nm= ge[ge.length- 1].getStableID();				// check last gene
			for (pos = 0; pos< names.length; pos++) 
				if (nm.equalsIgnoreCase(names[pos]))
					break;
			if (pos>= names.length) { 
				System.err.println(nm+ " not found in reference data set.");
				tst= false;
			} else {
				try{
					String test= g.readSequence(ge[ge.length- 1]);
					if (!seqs[pos].equalsIgnoreCase(test)) { 
						System.err.println(nm+ " mismatches the reference:\n\t"
								+ seqs[pos]+ "\n\t"+ test);
						tst= false;
					} else
						System.out.println(nm+ " passed test.");
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
			
		}
		
		return tst;
	}
	
	public Graph getGraphEncode() {
		
		if (graph == null) {
			System.out.println(Constants.getDateString()+ " loading Graph");
			graph= GraphHandler.readIn(
					Constants.HOME_DIR+ File.separator+ GraphHandler.GRAPH_SUBDIR+ File.separator+ GraphHandler.GRAPH_ENCODE_FNAME);
			if (graph!= null) 
				return graph;
			
				// rebuild
			if (SPECIES== null || SPECIES.length< 1)	// or (< 2) ?
				return null;
			
			System.out.println("["+new Date(System.currentTimeMillis())+"] -- start building graph --");
			Species[] spec= new Species[SPECIES.length];
			for (int i = 0; i < spec.length; i++) 
				spec[i]= new Species(SPECIES[i]);
				
			graph= new Graph(spec);
			
			System.out.println("["+new Date(System.currentTimeMillis())+"] "+graph.getSpecies().length+ " species.");
			GeneHomology[] homols= retrieveHomologGenesReference_mart(graph, Species.SP_NAMES_BINOMIAL[0]);
			System.out.println("retrieved "+homols.length+" homologies.");
			System.out.println("["+new Date(System.currentTimeMillis())+"] filtering graph.");
			homols= graph.filterForRegions(readEncodeRegions(graph.getSpeciesByName("homo_sapiens"), EncodeRegion.ENC_REG_HUMAN));	// reference
			GeneHomology[] homOld= homols;
			GeneHomology[] homNew= retrieveHomologGenesIterative_mart(graph, Species.SP_NAMES_BINOMIAL[0]);
			homols= new GeneHomology[homOld.length+ homNew.length];
			for (int i = 0; i < homOld.length; i++) 
				homols[i]= homOld[i];
			for (int i = 0; i < homNew.length; i++) 
				homols[i+ homOld.length]= homNew[i];
			System.out.println("retrieved "+homols.length+" homologies.");
			
			
			System.out.println("["+new Date(System.currentTimeMillis())+"] retrieving info for "+homols.length+" homologies.");
			homols= retrieveHomologyInfo_single(homols);			// cannot join two databases (mart, compara)
			
			System.out.println("["+new Date(System.currentTimeMillis())+"] retrieving info for "+graph.getGenes().length+" genes.");
			graph= retrieveGeneInfo(graph);
			System.out.println("["+new Date(System.currentTimeMillis())+"] # finished building graph #");
			

			
			System.out.println(Constants.getDateString()+ " Graph downloaded");
			tempExonMap= null;	// remove temporal references
			tempTranscriptMap= null;
			System.gc();

			GraphHandler.writeOut(graph,
					Constants.HOME_DIR+ File.separator+ GraphHandler.GRAPH_SUBDIR+ File.separator+ GraphHandler.GRAPH_ENCODE_FNAME); 		// writeGraph();
			System.out.println(Constants.getDateString()+ " Graph written");
		}
		
//		Gene g0= (Gene) graph.getSpecies()[0].getGeneIterator().next();
//		boolean init= true;
//		for (int i = 0; i < g0.getExons().length; i++) 
//			if (g0.getExons()[i].getHomologs()!= null) {
//				init= false;
//				break;
//			}
//		if (init) {
//			graph.init();	// init graph
//			System.out.println("Graph inited ---");
//			graph.writeOut();
//		}

		return graph;
	}
	
	
	EncodeRegion[] readEncodeRegions(Species spe, String fName) {
		
		Vector regionsVec= new Vector(); 
		try {
				// ENm001  chr7    115404472       117281897
			BufferedReader buffy= new BufferedReader(new FileReader(fName));
			buffy.readLine();	// descriptor: #name...
			while (buffy.ready()) {
				String line= buffy.readLine();				
				StringTokenizer toki= new StringTokenizer(line);
				EncodeRegion reg= new EncodeRegion(spe);
				reg.setSilentAndNumber(toki.nextToken());
				reg.setChromosome(toki.nextToken().trim().substring(3));	// eliminate "chr"
				reg.setStart(Integer.parseInt(toki.nextToken()));
				reg.setEnd(Integer.parseInt(toki.nextToken()));
				regionsVec.add(reg);
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return EncodeRegion.toArray(regionsVec);
		
	}
	
	public Graph getGraphAllGenes(String specName) {
		return getGraphAllGenes(new Species(specName));
	}
	
	public Graph getGraphAllGenes(Species spec) {
		
		System.out.println(Constants.getDateString()+ " loading Graph");
		graph= GraphHandler.readIn(GraphHandler.getGraphAbsPath(spec));
		if (graph== null) {
				// rebuild
			if (spec== null)
				return null;
			
			System.out.println("["+new Date(System.currentTimeMillis())+"] -- start building graph --");
			Species[] speci= new Species[] {spec};
			graph= new Graph(speci);
			
			System.out.println("["+new Date(System.currentTimeMillis())+"] "+graph.getSpecies().length+ " species.");
			for (int i = 0; i < speci.length; i++) 
				graph= retrieveGenesAll_mart(graph, speci[i]);
			
			
			System.out.println("["+new Date(System.currentTimeMillis())+"] retrieving info for "+graph.getGenes().length+" genes.");
			graph= retrieveGeneInfo(graph);
			System.out.println("["+new Date(System.currentTimeMillis())+"] # finished building graph #");
			

			
			System.out.println(Constants.getDateString()+ " Graph downloaded");
			tempExonMap= null;	// remove temporal references
			tempTranscriptMap= null;
			System.gc();

			GraphHandler.writeOut(graph, GraphHandler.getGraphAbsPath(spec)+ "_download"); 		// writeGraph();
			System.out.println(Constants.getDateString()+ " Graph written");

			filter(graph);
			
			GraphHandler.writeOut(graph, GraphHandler.getGraphAbsPath(spec)+ "filtered"); 
			System.out.println(Constants.getDateString()+ " Graph written");
			
			GraphHandler.writeOut(graph, GraphHandler.getGraphAbsPath(spec)); 		// working copy
			System.out.println(Constants.getDateString()+ " Graph written");
		}
		return graph;		
	}

	static void filter(Graph g) {
		System.out.println("pw as variations: "+g.countASVariations());
		g.filterNonsense();	// first remove transcripts, --> later genes have less transcripts
		g.filterSingleTranscriptGenes();
		System.out.println("pw as variations: "+g.countASVariations());
		GraphHandler.writeOut(g);
	}
	
	public Graph getGraphAllHomologs(String[] specNames) {
		
		if (graph == null) {
			System.out.println(Constants.getDateString()+ " loading Graph");
			graph= GraphHandler.readIn(GraphHandler.getGraphAbsPath());
			if (graph!= null) 
				return graph;
			
				// rebuild
			if (specNames== null || specNames.length< 1)	// or (< 2) ?
				return null;
			
			System.out.println("["+new Date(System.currentTimeMillis())+"] -- start building graph --");
			Species[] spec= new Species[specNames.length];
			for (int i = 0; i < spec.length; i++) 
				spec[i]= new Species(specNames[i]);
				
			graph= new Graph(spec);
			
			System.out.println("["+new Date(System.currentTimeMillis())+"] "+graph.getSpecies().length+ " species.");
			GeneHomology[] homols= retrieveHomologGenesAll_mart(graph);
			System.out.println("["+new Date(System.currentTimeMillis())+"] retrieving info for "+homols.length+" homologies.");
			//homols= retrieveHomologyInfo_single(homols);			// cannot join two databases (mart, compara)
			
			System.out.println("["+new Date(System.currentTimeMillis())+"] retrieving info for "+graph.getGenes().length+" genes.");
			graph= retrieveGeneInfo(graph);
			System.out.println("["+new Date(System.currentTimeMillis())+"] # finished building graph #");
			

			
			System.out.println(Constants.getDateString()+ " Graph downloaded");
			tempExonMap= null;	// remove temporal references
			tempTranscriptMap= null;
			System.gc();

			GraphHandler.writeOut(graph, GraphHandler.getGraphAbsPath()); 		// writeGraph();
			System.out.println(Constants.getDateString()+ " Graph written");
		}
		
//		Gene g0= (Gene) graph.getSpecies()[0].getGeneIterator().next();
//		boolean init= true;
//		for (int i = 0; i < g0.getExons().length; i++) 
//			if (g0.getExons()[i].getHomologs()!= null) {
//				init= false;
//				break;
//			}
//		if (init) {
//			graph.init();	// init graph
//			System.out.println("Graph inited ---");
//			graph.writeOut();
//		}

		return graph;
	}

	static void outputVariantHomologs(Graph g) {
		
		PrintStream p= System.out;
		
		Gene[] ge= g.getSpecies()[0].getGenes();	// does not matter which spec is base
		for (int i = 0; i < ge.length; i++) {
			int grad= 0;
			for (int j = 0; j < g.getSpecies().length; j++) {
				Gene ge1= null;
				if (g.getSpecies()[j]== ge[i].getSpecies())
					ge1= ge[i];
				else
					ge1= ge[i].getHomologies(g.getSpecies()[j])[0].getOtherGene(ge[i]);
				for (int k = (j+1); k < g.getSpecies().length; k++) {
					Gene ge2= ge[i].getHomologies(g.getSpecies()[k])[0].getOtherGene(ge[i]);
					if (!ge1.toStringSSPattern().equals(ge2.toStringSSPattern()))
						++grad;
				}
			}
			
			if (grad> 0&& ge[i].toStringSSPattern().length()<= 18) {
				p.println(grad);
				for (int j = 0; j < g.getSpecies().length; j++) {
					if (g.getSpecies()[j]== ge[i].getSpecies()) {
						p.print(ge[i].toStringSSPattern()+"\t"+ ge[i].getSpecies().getCommonName()
								+"("+ ge[i].getGeneID()+ ") "
								+ ge[i].getTranscriptCount()+ " [");
						for (int k = 0; k < ge[i].getTranscriptCount(); k++) 
							p.print(ge[i].getTranscripts()[k].getExons().length);
						p.println("]");
					} else {
						Gene gx= ge[i].getHomologies(g.getSpecies()[j])[0].getOtherGene(ge[i]);
						p.print(gx.toStringSSPattern()+"\t"+ gx.getSpecies().getCommonName()
								+"("+ gx.getGeneID()+ ")"
								+ gx.getTranscriptCount()+ " [");
						for (int k = 0; k < gx.getTranscriptCount(); k++) 
							p.print(gx.getTranscripts()[k].getExons().length);
						p.println("]");
					}
				}
			}
		}

	}
	
	static void removeNotAllHomologGenes(Graph g) {
		Vector v= new Vector();
		for (int x = 0; x < g.getSpecies().length; x++) {
			Gene[] ge= g.getSpecies()[x].getGenes();
			for (int i = 0; i < ge.length; i++) {
				int j;
				for (j = 0; j < g.getSpecies().length; j++) {
					if (j== x)
						continue;
					if (ge[i].getHomologies(g.getSpecies()[j])== null)
						break;
				}
				if (j< g.getSpecies().length) {
					v.add(ge[i]);
					for (int k = 0; k < g.getSpecies().length; k++) {
						GeneHomology[] hgg= ge[i].getHomologies(g.getSpecies()[k]);
						for (int m = 0; hgg!= null&& m < hgg.length; m++) {
							Gene hg= hgg[m].getOtherGene(ge[i]);
//							if (hg.getHomologies(ge[i].getSpecies())== null||
//									hg.getHomologies(ge[i].getSpecies()).length< 1)
								v.add(hg);
						}
					}
				}
			}
		}
		for (int i = 0; i < v.size(); i++) 
			g.removeKill((Gene) v.elementAt(i)); 
		//GraphHandler.writeOut(g, GraphHandler.getGraphAbsPath()); 		// writeGraph();
	}
	
	public static void main(String[] args) {
			
			EnsemblDBAdaptor adaptor= new EnsemblDBAdaptor();
	//		con= adaptor.connect("homo_sapiens_core");	// homo_sapiens_core_29_35b
	//		adaptor.testStatement(con);
			
	//		adaptor.testGraphEncode(adaptor.getGraphEncode());
			
			
			Graph g= adaptor.getGraphAllHomologs(SPECIES_ISMB);
			g.filterNonCodingTranscripts();

			
			try { 
				PrintStream p= System.out;
				Gene tst= g.getSpecies()[0].getGenes()[0];
				p.println(tst.getSpecies().getCommonName()+"\t"+tst.getGeneID());
				for (int i = 0; i < tst.getExons().length; i++) {
					p.println(">"+ tst.getExons()[i].getExonID());
					p.println(Graph.readSequence(tst.getExons()[i]));
				}
				p.flush(); p.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
//			try {
//				PrintStream p= new PrintStream("vars_out");
//				for (int i = 0; i < g.getSpecies().length; i++) {
//					p.println(g.getSpecies()[i].getCommonName());
//					ASVariation[][] classes= g.getSpecies()[i].getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);
//					Method m = classes[0][0].getClass().getMethod("isTrue", null);
//					ASVariation[][] filtClasses= ASAnalyzer.filter(classes, m);
//					gphase.tools.Arrays.sort2DFieldRev(filtClasses);
//					ASAnalyzer.outputVariations(filtClasses, false, false, p);
//					p.println("\n");
//				}
//				p.flush(); p.close();
//			} catch (Exception e) {
//				e.printStackTrace();
//			}
			
//			if (args.length< 1) {
//				System.err.println("at least one species needed");
//				System.exit(-1);
//			}
//			
//			Graph g= null;
//			if (args.length== 1) {
//				g= adaptor.getGraphAllGenes(args[0]);
//				g.filterNonCodingTranscripts();
//				try {
//					PrintStream p= new PrintStream("structures");
//					ASAnalyzer.test04_determineVariations(g, p);
//					p.flush(); p.close();
//				} catch (Exception e) {
//					e.printStackTrace();
//				}
//			} else {
//				g= adaptor.getGraphAllHomologs(args);
//			}
//			
//			g.initSpliceSiteHomology();
	}
}
