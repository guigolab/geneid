package org.ensembl.driver.plugin.compara;

import java.sql.PreparedStatement;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.ensembl.datamodel.Persistent;
import org.ensembl.datamodel.compara.DnaFragment;
import org.ensembl.datamodel.compara.GenomeDB;
import org.ensembl.datamodel.compara.GenomicAlign;
import org.ensembl.datamodel.compara.MethodLink;
import org.ensembl.datamodel.impl.PersistentImpl;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.compara.DnaFragmentAdaptor;
import org.ensembl.driver.compara.GenomeDBAdaptor;
import org.ensembl.driver.compara.GenomicAlignAdaptor;
import org.ensembl.driver.compara.MethodLinkAdaptor;


/**
 * Fetches DNA-DNA alignment information from the compara database.
 * @author <a href="maito:vvi@sanger.ac.uk">Vivek Iyer</a>
**/
public class 
  MySQLGenomicAlignAdaptor 
extends 
  ComparaBaseAdaptor
implements 
  GenomicAlignAdaptor
{
  public static String TABLE_NAME = "genomic_align_block";
  public static String CONSENSUS_DNA_FRAG_ID = TABLE_NAME+"."+"consensus_dnafrag_id";
  public static String CONSENSUS_START = TABLE_NAME+"."+"consensus_start";
  public static String CONSENSUS_END = TABLE_NAME+"."+"consensus_end";
  public static String QUERY_DNA_FRAG_ID = TABLE_NAME+"."+"query_dnafrag_id";
  public static String QUERY_START = TABLE_NAME+"."+"query_start";
  public static String QUERY_END = TABLE_NAME+"."+"query_end";
  public static String SCORE = TABLE_NAME+"."+"score";
  public static String PERC_ID = TABLE_NAME+"."+"perc_id";
  public static String CIGAR_LINE = TABLE_NAME+"."+"cigar_line";
  public static String METHOD_LINK_ID = TABLE_NAME+"."+"method_link_id";
  
  public static int DEFAULT_MAX_ALIGNMENT = 20000;

  public MySQLGenomicAlignAdaptor(ComparaMySQLDriver driver) {
    super(driver);
  }//end MySQLGenomicAlignAdaptor

  public String getType(){
    return TYPE;
  }//end getType
    
  public int store(GenomicAlign genomicAlign) throws  AdaptorException {
    super.store(genomicAlign);
    return 0;
  }//end store

  public GenomicAlign fetch(long internalID ) throws AdaptorException {
    return (GenomicAlign) super.fetch(new Long(internalID));
  }//end fetch
  
  /**
   * Fetch all genomic aligns matching the parameters.
   *
   * @param dnaFragment the input dnaFragment (implicitly
   * carrying the concensus genomeDB and the targetSpecies)
   * @param targetGenome genome to find alignments in
   * @param start bounded start coord
   * @param end bounded end coord
   * @param methodLinkType input method link type. if null, fetch regardless of method link type.
  **/
  public List fetch(
    DnaFragment dnaFragment, //implicit source species
    GenomeDB targetGenome, 
    int start,
    int end,
    String methodLinkType
  ) throws AdaptorException {

    StringBuffer select = new StringBuffer();
    StringBuffer whereClause = null;
    List arguments = new ArrayList();
    GenomeDBAdaptor genomeDBAdaptor = 
      (GenomeDBAdaptor)getDriver().getAdaptor(GenomeDBAdaptor.TYPE);

    MethodLinkAdaptor methodLinkAdaptor = null;
    
    int lowerBound = 0;
    int upperBound = 0;
    PreparedStatement statement;
    List returnList = new ArrayList();
    boolean inputDnaFragBelongsToConsensusAndTargetGenomeDBBelongsToQuery = true; 
    boolean inputDnaFragBelongsToQueryAndTargetGenomeDBBelongsToConsensus = false;
    List temporaryList;
    MethodLink methodLink = null;
    
    if(dnaFragment == null){
      throw new AdaptorException("Must supply dnaFragment for Genomic Align query");
    }//end if
    
    if(methodLinkType != null){
      methodLinkAdaptor = (MethodLinkAdaptor)getDriver().getAdaptor(MethodLinkAdaptor.TYPE);
      methodLink = methodLinkAdaptor.fetch(methodLinkType);
    }//end if
    
    //
    //If no target, prepare two selects - one which grabs all align rows where
    //the dnaFrag's genome appears as a consensus, and the other where it appears as a query.
    if(targetGenome == null){

      //
      //CONSENSUS-side select
      whereClause = addEqualsClause(CONSENSUS_DNA_FRAG_ID, new StringBuffer());
      arguments.add(new Long(dnaFragment.getInternalID()));
      
      if(methodLink != null){
        addEqualsClause(METHOD_LINK_ID, whereClause);
        arguments.add(new Long(methodLink.getInternalID()));
      }//end if
      
      if(start > 0 && end > 0){
        
        lowerBound = start - DEFAULT_MAX_ALIGNMENT;
        
        /*
               $sql .= ( " AND gab.consensus_start <= $end
                         AND gab.consensus_start >= $lower_bound
                         AND gab.consensus_end >= $start" );
        */
        
        addLEClause(CONSENSUS_START, whereClause);
        arguments.add(new Integer(end));
        
        addGEClause(CONSENSUS_START, whereClause);
        arguments.add(new Integer(lowerBound));
        
        addGEClause(CONSENSUS_END, whereClause);
        arguments.add(new Integer(start));

      }//end if
      
      statement = prepareSelectWithWhereClause(whereClause.toString());
      returnList = executeStatementAndConvertResultToPersistent(statement, arguments);
      inflateGenomicAlignBlocks(returnList, false);//DONT invert the genomic align information!

      //
      //Plus result of query-side select.
      whereClause = addEqualsClause(QUERY_DNA_FRAG_ID, new StringBuffer());
      arguments.clear();
      arguments.add(new Long(dnaFragment.getInternalID()));
      
      if(start > 0 && end > 0){
        
        addLEClause(QUERY_START, whereClause);
        arguments.add(new Integer(end));
        
        addGEClause(QUERY_START, whereClause);
        arguments.add(new Integer(lowerBound));
        
        addGEClause(QUERY_END, whereClause);
        arguments.add(new Integer(start));

      }//end fi
      
      statement = prepareSelectWithWhereClause(whereClause.toString());
      temporaryList = executeStatementAndConvertResultToPersistent(statement, arguments);
      inflateGenomicAlignBlocks(temporaryList, true);//DO INVERT the genomic align information!
      if(temporaryList.size() > 0){
        returnList.addAll(temporaryList);
      }//end if
      
    }else{
      //
      //If there IS a target, then establish whether the dnaFrag's genome is related
      //to the target as consensus->query or query->consensus. Then phrase the sql
      //accordingly.
      inputDnaFragBelongsToConsensusAndTargetGenomeDBBelongsToQuery = 
        genomeDBAdaptor.firstArgumentIsKnownConsensusAndSecondIsKnowQuery(
          dnaFragment.getGenomeDB(), 
          targetGenome
        );

      inputDnaFragBelongsToQueryAndTargetGenomeDBBelongsToConsensus = 
        genomeDBAdaptor.firstArgumentIsKnownConsensusAndSecondIsKnowQuery(
          targetGenome, 
          dnaFragment.getGenomeDB()
        );

      if(inputDnaFragBelongsToConsensusAndTargetGenomeDBBelongsToQuery){
        //
        //create join so as to pin the QUERY dnafrag id's by having to 
        //belong to the target genome db id
        select = 
          createJoinedSelect(
            new String[]{TABLE_NAME, MySQLDnaFragmentAdaptor.TABLE_NAME},
            new String[]{QUERY_DNA_FRAG_ID, MySQLDnaFragmentAdaptor.DNAFRAG_ID}
          );

        if(methodLink != null){
          addEqualsClause(METHOD_LINK_ID, select);
          arguments.add(new Long(methodLink.getInternalID()));
        }//end if
          
        //
        //Pin the consensus dnafrag id by the input dnafrag id
        addEqualsClause(CONSENSUS_DNA_FRAG_ID, select);
        arguments.add(new Long(dnaFragment.getInternalID()));

        //
        //pin the QUERY dnafrag id's by having to belong to the target genome db id
        addEqualsClause(MySQLDnaFragmentAdaptor.GENOME_DB_ID, select);
        arguments.add(new Long(targetGenome.getInternalID()));

        if(start >0 && end >0){
          lowerBound = start - DEFAULT_MAX_ALIGNMENT;

          addLEClause(CONSENSUS_START, select);
          arguments.add(new Integer(end));

          addGEClause(CONSENSUS_START, select);
          arguments.add(new Integer(lowerBound));

          addGEClause(CONSENSUS_END, select);
          arguments.add(new Integer(start));
        }//end if
          
      }else if(inputDnaFragBelongsToQueryAndTargetGenomeDBBelongsToConsensus){
        //
        //create join so as to pin the CONSENSUS dnafrag id's by having to 
        //belong to the target genome db id
        select = 
          createJoinedSelect(
            new String[]{TABLE_NAME, MySQLDnaFragmentAdaptor.TABLE_NAME},
            new String[]{CONSENSUS_DNA_FRAG_ID, MySQLDnaFragmentAdaptor.DNAFRAG_ID}
          );
          
        if(methodLink != null){
          addEqualsClause(METHOD_LINK_ID, select);
          arguments.add(new Long(methodLink.getInternalID()));
        }//end if
          
        //
        //Pin the QUERY dnafrag id by the input dnafrag id
        addEqualsClause(QUERY_DNA_FRAG_ID, select);
        arguments.add(new Long(dnaFragment.getInternalID()));
        
        //
        //pin the CONSENSUS dnafrag id's by having to belong to the target genome db id
        addEqualsClause(MySQLDnaFragmentAdaptor.GENOME_DB_ID, select);
        arguments.add(new Long(targetGenome.getInternalID()));
          
        if(start >0 && end >0){
          lowerBound = start - DEFAULT_MAX_ALIGNMENT;
          
          addLEClause(QUERY_START, select);
          arguments.add(new Integer(end));

          addGEClause(QUERY_START, select);
          arguments.add(new Integer(lowerBound));

          addGEClause(QUERY_END, select);
          arguments.add(new Integer(start));
        }//end if
      }else{
        throw new AdaptorException("The compara data doesnt have a direct comparison of "+dnaFragment.getGenomeDB().getName()+" against "+targetGenome.getName());
      }//end if (was source genome stored as query or consensus?)
      
      statement = prepareStatement(select.toString());
      returnList = executeStatementAndConvertResultToPersistent(statement, arguments);
      
      if(inputDnaFragBelongsToQueryAndTargetGenomeDBBelongsToConsensus){
        inflateGenomicAlignBlocks(returnList, true);//invert the genomic align information!
      }else{
        inflateGenomicAlignBlocks(returnList, false);//DONT invert the genomic align information!
      }
    }//end if (was target genome not defined?)
    
    System.out.println("SIZE OF GENOMIC ALIGNS RETURNED: "+returnList.size());
    return returnList;
  }//end fetch
  
  /**
   * Input - uninflated GenomicAlignBlocks
   * Output - genomic align blocks with dnafrags inflated
  **/
  private void inflateGenomicAlignBlocks(List genomicAlignBlocks, boolean invertAlign) throws AdaptorException{
    //
    //We cache retrieved dnafrags to aid their 'fluffing up'. Oh, for a real
    //object cache...sigh...
    long dnaFragInternalId = 0;
    HashMap dnaFragmentCache = new HashMap();
    DnaFragmentAdaptor dnaFragAdaptor = 
      (DnaFragmentAdaptor)getDriver().getAdaptor(DnaFragmentAdaptor.TYPE);
    DnaFragment dnaFragment;
    Object object = null;
    GenomicAlign genomicAlign = null;
    MethodLinkAdaptor methodLinkAdaptor = (MethodLinkAdaptor)getDriver().getAdaptor(MethodLinkAdaptor.TYPE);
    
    //
    //Fluff up the returned objects with the DnaFrags they have a handle to. 
    for(int i=0; i<genomicAlignBlocks.size(); i++){
      genomicAlign = (GenomicAlign)genomicAlignBlocks.get(i);
      
      dnaFragInternalId = genomicAlign.getConsensusDnaFragmentId();
      if(genomicAlign.getConsensusDnaFragment() == null){
        dnaFragment = (DnaFragment)dnaFragmentCache.get(new Long(dnaFragInternalId));
        if(dnaFragment == null){
          dnaFragment = dnaFragAdaptor.fetch(dnaFragInternalId);
          dnaFragmentCache.put(new Long(dnaFragInternalId), dnaFragment);
        }//end if
        genomicAlign.setConsensusDnaFragment(dnaFragment);
      }//end if
      
      if(genomicAlign.getQueryDnaFragment() == null){
        dnaFragInternalId = 
          genomicAlign.getQueryDnaFragmentId();
        dnaFragment = (DnaFragment)dnaFragmentCache.get(new Long(dnaFragInternalId));
        if(dnaFragment == null){
          dnaFragment = dnaFragAdaptor.fetch(dnaFragInternalId);
          dnaFragmentCache.put(new Long(dnaFragInternalId), dnaFragment);
        }//end if
        genomicAlign.setQueryDnaFragment(dnaFragment);
      }//end if

      if(genomicAlign.getMethodLinkInternalId() >0){
        genomicAlign.setMethodLink(
          (MethodLink)
          methodLinkAdaptor.fetch(genomicAlign.getMethodLinkInternalId())
        );
      }
      
      if(invertAlign){
        DnaFragment frag = genomicAlign.getConsensusDnaFragment();
        int dnafragid = genomicAlign.getConsensusDnaFragmentId();
        int start = genomicAlign.getConsensusStart();
        int end = genomicAlign.getConsensusEnd();
        genomicAlign.setConsensusDnaFragment(genomicAlign.getQueryDnaFragment());
        genomicAlign.setConsensusDnaFragmentId(genomicAlign.getQueryDnaFragmentId());
        genomicAlign.setConsensusStart(genomicAlign.getQueryStart());
        genomicAlign.setConsensusEnd(genomicAlign.getQueryEnd());
        genomicAlign.setQueryDnaFragment(frag);
        genomicAlign.setQueryDnaFragmentId(dnafragid);
        genomicAlign.setQueryStart(start);
        genomicAlign.setQueryEnd(end);
      }//end if
      
    }//end for
  }//end inflateGenomicAlignBlocks
    
  protected String getTableName(){
    return TABLE_NAME;
  }//end getTableName
  
  protected PersistentImpl createNewObject() {
    return (PersistentImpl)getFactory().createGenomicAlign();
  }//end createNewObject() 
  
  /**
   * Populate object attributes using the columns map.
   * @param columns Hashmap of column names and values. 
   * @param object PersistentImpl passed in to have its attributes populated by these values.
  **/
  protected void mapColumnsToObject(HashMap columns, Persistent object) {
    GenomicAlign genomicAlign = (GenomicAlign)object;
    genomicAlign.setConsensusDnaFragmentId(((Integer)columns.get(CONSENSUS_DNA_FRAG_ID)).intValue());
    genomicAlign.setConsensusStart(((Integer)columns.get(CONSENSUS_START)).intValue());
    genomicAlign.setConsensusEnd(((Integer)columns.get(CONSENSUS_END)).intValue());
    genomicAlign.setQueryDnaFragmentId(((Integer)columns.get(QUERY_DNA_FRAG_ID)).intValue());
    genomicAlign.setQueryStart(((Integer)columns.get(QUERY_START)).intValue());
    genomicAlign.setQueryEnd(((Integer)columns.get(QUERY_END)).intValue());
    genomicAlign.setScore(((Double)columns.get(SCORE)).doubleValue());
    if(columns.get(PERC_ID) != null){
      genomicAlign.setPercentageId(((Integer)columns.get(PERC_ID)).intValue());
    }
    genomicAlign.setCigarString((String)columns.get(CIGAR_LINE));
    genomicAlign.setMethodLinkInternalId(((Integer)columns.get(METHOD_LINK_ID)).intValue());
  }//end mapColumnsToObject
  
  /**
   * Input; genomeDB with populated attributes. Output - hashmap of columns
   * and their values appropriate for an insert/update.
  **/
  protected HashMap mapObjectToColumns(Persistent object) {
    HashMap columns = new HashMap();
    GenomicAlign genomicAlign = (GenomicAlign)object;
    
    columns.put(CONSENSUS_DNA_FRAG_ID, new Integer(genomicAlign.getConsensusDnaFragmentId()));
    columns.put(CONSENSUS_START, new Integer(genomicAlign.getConsensusDnaFragment().getLocation().getStart()));
    columns.put(CONSENSUS_END, new Integer(genomicAlign.getConsensusDnaFragment().getLocation().getEnd()));
    
    columns.put(QUERY_DNA_FRAG_ID, new Integer(genomicAlign.getQueryDnaFragmentId()));
    columns.put(QUERY_START, new Integer(genomicAlign.getQueryDnaFragment().getLocation().getStart()));
    columns.put(QUERY_END, new Integer(genomicAlign.getQueryDnaFragment().getLocation().getEnd()));
    
    columns.put(SCORE, new Double(genomicAlign.getScore()));
    columns.put(PERC_ID, new Integer(genomicAlign.getPercentageId()));
    columns.put(CIGAR_LINE, genomicAlign.getCigarString());
    columns.put(METHOD_LINK_ID, new Integer(genomicAlign.getMethodLinkInternalId()));
    return columns;
  }//end mapObjectToColumns
  
  public HashMap getLogicalKeyPairs(Persistent genomeDB) throws AdaptorException{
    HashMap logicalKeyPairs = new HashMap();
    return logicalKeyPairs;
  }//end getLogicalKeyPairs
  
  public void validate(Persistent object) throws AdaptorException{
    GenomeDB genomeDB = (GenomeDB)object;
    
    if(genomeDB.getInternalID() <= 0){
      throw new AdaptorException("Attempt to store genomeDB "+genomeDB.getName()+" with missing id");
    }//end if
    
    if(genomeDB.getName() == null){
      throw new AdaptorException("Attempt to store genomeDB "+genomeDB.getInternalID()+" with missing name");
    }//end if
    
    if(genomeDB.getAssembly() == null){
      throw new AdaptorException("Attempt to store genomeDB "+genomeDB.getInternalID()+" with missing assembly");
    }//end if
  }//end validate
}
