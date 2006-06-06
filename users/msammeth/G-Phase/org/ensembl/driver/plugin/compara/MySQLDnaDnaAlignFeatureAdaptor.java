package org.ensembl.driver.plugin.compara;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import org.ensembl.datamodel.DnaDnaAlignFeature;
import org.ensembl.datamodel.Feature;
import org.ensembl.datamodel.Location;
import org.ensembl.datamodel.Persistent;
import org.ensembl.datamodel.compara.DnaFragment;
import org.ensembl.datamodel.compara.GenomeDB;
import org.ensembl.datamodel.compara.GenomicAlign;
import org.ensembl.datamodel.impl.PersistentImpl;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.compara.DnaDnaAlignFeatureAdaptor;
import org.ensembl.driver.compara.DnaFragmentAdaptor;
import org.ensembl.driver.compara.GenomeDBAdaptor;
import org.ensembl.driver.compara.GenomicAlignAdaptor;
import org.ensembl.util.NotImplementedYetException;


/**
 * Fetches DnaDnaAlignFeature objects.
 * @author <a href="maito:vvi@sanger.ac.uk">Vivek Iyer</a>
 * Need fetch(Location), fetch(internalId), and some mess involving analysis_id_conditions (!!)
**/
public class 
  MySQLDnaDnaAlignFeatureAdaptor 
extends 
  ComparaBaseAdaptor
implements 
  DnaDnaAlignFeatureAdaptor
{
  private HashMap dnaFragmentHash = new HashMap();

  public MySQLDnaDnaAlignFeatureAdaptor(ComparaMySQLDriver driver) {
    super(driver);
  }//end MySQLDnaDnaAlignFeatureAdaptor

  public String getType(){
    return TYPE;
  }
    
  /**
   * Not implemented yet!
  **/
  public void store( Feature feature ) throws  AdaptorException {
    throw new NotImplementedYetException("Store not implemented on DnaDnaAlignFeatureAdaptor");
  }//end store

  public List fetch(
      String consensusSpecies, 
      Location consensusLocation, 
      String querySpecies
  ) throws AdaptorException {
    
    return fetch(
      consensusSpecies, 
      consensusLocation, 
      querySpecies,
      null
    );
    
  }
  
  /**
   * <p>Returns a list of Dna-Dna alignment features for the input
   * query species, chromosome location and hit species. Namely, </p>
   * <ol>
   * <li> Fetches all DnaFragments for the input location and species </li>
   * <li> For each DnaFragment - retrieve the genomic aligns associated to that fragment
   * and consensus/query species</li>
   * <li> For each retrieved genomic align, convert to a dnadnaalignfeature, and pass out.</li>
   * </ol>
  **/
  public List fetch(
      String consensusSpecies, 
      Location consensusLocation, 
      String querySpecies,
      String methodLinkType
  ) throws AdaptorException {
    
    List returnList = new ArrayList();
    List dnaFragments;
    Iterator dnaFragmentIterator;
    DnaFragment currentFragment;
    List genomicAlignBlocksForDnaFragment;
    List dnaDnaAlignFeaturesForAlignBlockList;
    int startRelativeToDnaFrag = 0;
    int endRelativeToDnaFrag = 0;
    int lengthOfRequestedSegment = 0;
    
    GenomeDBAdaptor genomeDBAdaptor = 
      (GenomeDBAdaptor)getDriver().getAdaptor(MySQLGenomeDBAdaptor.TYPE);

    DnaFragmentAdaptor dnaFragmentAdaptor = 
      (DnaFragmentAdaptor)getDriver().getAdaptor(MySQLDnaFragmentAdaptor.TYPE);
    
    GenomicAlignAdaptor genomicAlignAdaptor =
      (GenomicAlignAdaptor)getDriver().getAdaptor(MySQLGenomicAlignAdaptor.TYPE);
      
    GenomeDB consensusGenome = genomeDBAdaptor.fetch(consensusSpecies);
    GenomeDB targetGenome = null;
    
    if(querySpecies != null){
      targetGenome = genomeDBAdaptor.fetch(querySpecies);
    }
    
    getLogger().fine("Fetching dna frags:");
    dnaFragments = 
      dnaFragmentAdaptor.fetch(
        consensusGenome, 
        "Chromosome", 
        consensusLocation.getSeqRegionName(),
        consensusLocation.getStart(),
        consensusLocation.getEnd()
      );
    
    getLogger().fine("Finished fetching dna frags:"+dnaFragments.size());
    dnaFragmentIterator = dnaFragments.iterator();

    while(dnaFragmentIterator.hasNext()){
      currentFragment = (DnaFragment)dnaFragmentIterator.next();

      //
      //It's understood here that start/end for the currentFragment are GLOBAL
      //starts and ends - that is, they are in chromosome coordinates already.
      startRelativeToDnaFrag = 
        consensusLocation.getStart() - currentFragment.getLocation().getStart() + 1;
      
      endRelativeToDnaFrag = 
        consensusLocation.getEnd() - currentFragment.getLocation().getStart() + 1;
      
      lengthOfRequestedSegment = endRelativeToDnaFrag - startRelativeToDnaFrag;

      if(startRelativeToDnaFrag < 1){
        startRelativeToDnaFrag = 1;
      }//end if

      /*
      if(endRelativeToDnaFrag > lengthOfRequestedSegment){
        endRelativeToDnaFrag = lengthOfRequestedSegment;
      }//end if
       */
      
      //
      //Cache the fragment we retrieve - this will be handy after the align blocks are
      //retrieved.
      getLogger().fine("fetching genomic aligns for dna frag");

      genomicAlignBlocksForDnaFragment = 
        genomicAlignAdaptor.fetch(
          currentFragment,
          targetGenome,
          startRelativeToDnaFrag,
          endRelativeToDnaFrag,
          methodLinkType
        );

      getLogger().fine("fetched "+genomicAlignBlocksForDnaFragment.size()+" blocks");

      dnaDnaAlignFeaturesForAlignBlockList = 
        createDnaDnaAlignFeaturesForAlignBlocks(
          consensusSpecies,
          consensusLocation.getStart(),
          consensusLocation.getEnd(),
          querySpecies,
          genomicAlignBlocksForDnaFragment,
          genomeDBAdaptor
        );

      returnList.addAll(dnaDnaAlignFeaturesForAlignBlockList);
    }//end while
    
    return returnList;
  }//end fetch
  
  /**
   * Copy each genomic align directly into a dna-dna-align-feature.
  **/
  private List createDnaDnaAlignFeaturesForAlignBlocks(
    String consensusSpecies,
    int consensusChromosomeStart,
    int consensusChromosomeEnd,
    String querySpecies,
    List genomicAlignBlocks,
    GenomeDBAdaptor genomeDBAdaptor
  ) throws AdaptorException {
    
    List returnList = new ArrayList();
    String tempQuerySpecies = null;
    GenomicAlign genomicAlign;

    DnaDnaAlignFeature dnaDnaAlign;
    Iterator iterator = genomicAlignBlocks.iterator();

    while(iterator.hasNext()){
      genomicAlign = (GenomicAlign)iterator.next();
      dnaDnaAlign = getFactory().createDnaDnaAlignFeature();
      if(querySpecies == null){
        dnaDnaAlign.setHitSpecies(genomicAlign.getQueryDnaFragment().getGenomeDB().getName());
      }else{
        dnaDnaAlign.setHitSpecies(querySpecies);
      }
      dnaDnaAlign.setScore(genomicAlign.getScore());
      dnaDnaAlign.setCigarString(genomicAlign.getCigarString());
      dnaDnaAlign.setSpecies(consensusSpecies);
      dnaDnaAlign.setDescription(null);
      dnaDnaAlign.setDisplayName(null);
      dnaDnaAlign.setMethodLinkType(genomicAlign.getMethodLink().getType());
      dnaDnaAlign.setPercentageIdentity(genomicAlign.getPercentageId());
      
      try{
        dnaDnaAlign.setLocation(
          new Location(
            "chromosome:"+
            genomicAlign.getConsensusDnaFragment().getName()+":"+
            genomicAlign.getConsensusStart()+"-"+
            genomicAlign.getConsensusEnd()+
            ":1"
          )
        );

      }catch(java.text.ParseException exception){
        //I find this use of RuntimeException distasteful, but 
        //we don't have a common FatalException subclass of RuntimeException!
        throw new RuntimeException(
          "I cant convert the string chromosome:"+
          genomicAlign.getConsensusDnaFragment().getName()+":"+
          genomicAlign.getConsensusStart()+"-"+
          genomicAlign.getConsensusEnd()+" to a location",
          exception
        );
      }
      
      try{
        dnaDnaAlign.setHitLocation(
          new Location(
            "I cant convert the string chromosome:"+
            genomicAlign.getQueryDnaFragment().getName()+":"+
            genomicAlign.getQueryStart()+"-"+
            genomicAlign.getQueryEnd()+":"+
            genomicAlign.getQueryStrand()
          )
        );
      }catch(java.text.ParseException exception){
        //I find this use of RuntimeException distasteful, but 
        //we don't have a common FatalException subclass of RuntimeException!
        throw new RuntimeException(
          "I cant convert the string chromosome:"+
          genomicAlign.getQueryDnaFragment().getName()+":"+
          genomicAlign.getQueryStart()+":"+
          genomicAlign.getQueryEnd()+":"+
          genomicAlign.getQueryStrand()+" to a location",
          exception
        );
      }

      returnList.add(dnaDnaAlign);
    }//end while
    
    return returnList;
  }//end createDnaDnaAlignFeaturesForAlignBlocks
  
  /**
   * A hash to cache DnaFragments based on internal id. As a stateful property,
   * this is harmless - very unlikely that there's a problem in it getting stored across
   * successive invocations.
  **/
  private HashMap getDnaFragmentHash(){
    return  dnaFragmentHash;
  }//end getDnaFragmentHash

  protected void configure(){
    //do nothing
  }
  
  protected  String getTableName(){
    if(true){
      throw new IllegalStateException("This adaptor should not be writing this table");
    }
    return null;
  }
  
  protected HashMap mapObjectToColumns(Persistent object){
    if(true){
      throw new IllegalStateException("This adaptor should not be writing this table");
    }
    return null;
  }
  
  protected void mapColumnsToObject(HashMap columns, Persistent object){
    if(true){
      throw new IllegalStateException("This adaptor should not be writing this table");
    }
  }
  
  protected PersistentImpl createNewObject(){
    if(true){
      throw new IllegalStateException("This adaptor should not be writing this table");
    }
    return null;
  }//end createNewObject
  
  protected void validate(Persistent object) throws AdaptorException{
    if(true){
      throw new IllegalStateException("This adaptor should not be writing this table");
    }
  }
  
  public HashMap getLogicalKeyPairs(Persistent genomeDB) throws AdaptorException{
    if(true){
      throw new IllegalStateException("This adaptor should not be writing this table");
    }
    return null;
  }
}//end MySQLDnaDnaAlignFeatureAdaptor 
