package org.ensembl.driver.plugin.compara;

import java.sql.PreparedStatement;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.ensembl.datamodel.Location;
import org.ensembl.datamodel.Persistent;
import org.ensembl.datamodel.compara.DnaFragment;
import org.ensembl.datamodel.compara.GenomeDB;
import org.ensembl.datamodel.impl.PersistentImpl;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.compara.DnaFragmentAdaptor;
import org.ensembl.driver.compara.GenomeDBAdaptor;

/**
 * @author <a href="maito:vvi@sanger.ac.uk">Vivek Iyer</a>
**/
public class MySQLDnaFragmentAdaptor extends ComparaBaseAdaptor implements DnaFragmentAdaptor{

  public static String TABLE_NAME = "dnafrag";
  
  public static String DNAFRAG_ID = TABLE_NAME+"."+"dnafrag_id";
  public static String START = TABLE_NAME+"."+"start";
  public static String END = TABLE_NAME+"."+"end";
  public static String NAME = TABLE_NAME+"."+"name";
  public static String GENOME_DB_ID = TABLE_NAME+"."+"genome_db_id";
  public static String DNAFRAG_TYPE = TABLE_NAME+"."+"dnafrag_type";
  

  public MySQLDnaFragmentAdaptor(ComparaMySQLDriver driver) {
    super(driver);
  }//end MySQLAlignAdaptor

  public String getType(){
    return TYPE;
  }//end getType

  public int store(DnaFragment object) throws AdaptorException{
    super.store(object);
    return 0;
  }
  
  public DnaFragment fetch( long internalID ) throws AdaptorException {
    Long id = new Long(internalID);
    DnaFragment fragment = (DnaFragment)super.fetch(id);
    inflateDnaFragment(fragment);
    return fragment;
  }//end 
  
  /**
   * Fetch genome db by species name
  **/
  public List fetch(
    GenomeDB genomeDB, 
    String dnaFragType, 
    String name, 
    int start, 
    int end
  )throws AdaptorException{
    StringBuffer whereClause = null;
    PreparedStatement statement = null;
    List arguments = new ArrayList();
    GenomeDBAdaptor genomeDBAdaptor = null;
    List returnList;
    DnaFragment dnaFragment;
    
    if(genomeDB == null){
      throw new AdaptorException("Genome DB must be specified for fetch");
    }//end if
    
    if(dnaFragType == null){
      throw new AdaptorException("DnaFrag type must be specified for fetch");
    }//end if
    
    whereClause = 
      addEqualsClause(GENOME_DB_ID, (new StringBuffer()));
    
    whereClause = 
      addEqualsClause(DNAFRAG_TYPE, whereClause);
    
    arguments.add(new Long(genomeDB.getInternalID()));
    arguments.add(dnaFragType);
      
    if(name != null){
      addEqualsClause(NAME, whereClause);
      arguments.add(name);
    }//end if
    
    if(start > 0){
      addGEClause(END, whereClause);
      arguments.add(new Integer(start));
    }//end if
    
    if(end > 0){
      addLEClause(START, whereClause);
      arguments.add(new Integer(end));
    }//end if
    
    statement = prepareSelectWithWhereClause(whereClause.toString());
    
    returnList = executeStatementAndConvertResultToPersistent(statement, arguments);

    for(int i=0; i< returnList.size(); i++){
      dnaFragment = (DnaFragment)returnList.get(i);
      dnaFragment.setGenomeDB(genomeDB);
    }//end for
    
    return returnList;
  }//end fetch

  /**
   * Input - dnafragment with no genome db object
   * output - dnafragment with genomeDb attached
  **/
  private void inflateDnaFragment(DnaFragment fragment) throws AdaptorException{
    GenomeDBAdaptor genomeDBAdaptor = 
      (GenomeDBAdaptor)getDriver().getAdaptor(GenomeDBAdaptor.TYPE);
    
    if(fragment.getGenomeDB() == null){
      fragment.setGenomeDB(genomeDBAdaptor.fetch(fragment.getGenomeDbInternalId()));
    }//end if
  }
  
  protected String getTableName(){
    return TABLE_NAME;
  }//end getTableName
  
  protected PersistentImpl createNewObject() {
    return (PersistentImpl)getFactory().createDnaFragment();
  }//end createNewObject() 
  
  /**
   * Does NOT gather the GenomeDB object - that's left to the fetch() methods.
  **/
  protected void mapColumnsToObject(HashMap columns, Persistent object) {
    DnaFragment dnaFragment = (DnaFragment)object;
    int start;
    int end;
    Location location; 
    
    dnaFragment.setInternalID(((Integer)columns.get(DNAFRAG_ID)).intValue());
    dnaFragment.setGenomeDbInternalId(((Integer)columns.get(GENOME_DB_ID)).intValue());
    dnaFragment.setType((String)columns.get(DNAFRAG_TYPE));
    dnaFragment.setName((String)columns.get(NAME));
    
    start = ((Integer)columns.get(START)).intValue();
    end = ((Integer)columns.get(END)).intValue();
    
    try{
      location = 
        new Location(
          "chromosome"+
          dnaFragment.getName()+":"+
          start+":"+
          end+":0"
        );
    }catch(java.text.ParseException exception){
      //I find this use of RuntimeException distasteful, but 
      //we don't have a common FatalException subclass of RuntimeException!
      throw new RuntimeException(
        "Fatal problem parsing location: chromosome"+
        dnaFragment.getName()+":"+
        start+":"+
        end+":0",
        exception
       );
    }

    dnaFragment.setLocation(location);
  }//end mapColumnsToObject
  
  /**
   * Input; genomeDB with populated attributes. Output - hashmap of columns
   * and their values appropriate for an insert/update.
  **/
  protected HashMap mapObjectToColumns(Persistent object) {
    HashMap columns = new HashMap();
    DnaFragment dnaFragment = (DnaFragment)object;
    columns.put(DNAFRAG_ID, new Long(dnaFragment.getInternalID()));
    columns.put(TYPE, new Integer(dnaFragment.getType()));
    columns.put(NAME, dnaFragment.getName());
    columns.put(START, new Integer(dnaFragment.getLocation().getStart()));
    columns.put(END, new Integer(dnaFragment.getLocation().getEnd()));
    columns.put(GENOME_DB_ID, new Integer(dnaFragment.getGenomeDbInternalId()));
    return columns;
  }//end mapObjectToColumns
  
  public HashMap getLogicalKeyPairs(Persistent object) throws AdaptorException{
    HashMap logicalKeyPairs = new HashMap();
    DnaFragment dnaFragment = (DnaFragment)object;
    logicalKeyPairs.put(TYPE, new Integer(dnaFragment.getType()));
    logicalKeyPairs.put(NAME, dnaFragment.getName());
    logicalKeyPairs.put(START, new Integer(dnaFragment.getLocation().getStart()));
    logicalKeyPairs.put(END, new Integer(dnaFragment.getLocation().getEnd()));
    return logicalKeyPairs;
  }//end getLogicalKeyPairs
  
  public void validate(Persistent object) throws AdaptorException{
    DnaFragment dnaFragment = (DnaFragment)object;
    if(dnaFragment.getInternalID() <= 0){
      throw new AdaptorException("Attempt to store dnaFragment "+dnaFragment.getName()+" with missing id");
    }//end if
    
    if(dnaFragment.getName() == null){
      throw new AdaptorException("Attempt to store dnaFragment "+dnaFragment.getInternalID()+" with missing name");
    }//end if   
    
    if(dnaFragment.getType() == null){
      throw new AdaptorException("Attempt to store dnaFragment "+dnaFragment.getName()+" with missing type");
    }//end if   
    
    if(dnaFragment.getLocation() == null){
      throw new AdaptorException("Attempt to store dnaFragment "+dnaFragment.getName()+" with missing start/end");
    }//end if   
    
    if(dnaFragment.getGenomeDB() == null){
      throw new AdaptorException("Attempt to store dnaFragment "+dnaFragment.getName()+" with missing GenomeDB");
    }//end if   
  }//end validate
}//end MySQLDnaFragmentAdaptor 
