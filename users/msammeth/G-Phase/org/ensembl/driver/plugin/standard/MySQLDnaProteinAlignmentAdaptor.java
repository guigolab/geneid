package org.ensembl.driver.plugin.standard;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;

import org.ensembl.datamodel.CoordinateSystem;
import org.ensembl.datamodel.DnaProteinAlignment;
import org.ensembl.datamodel.InvalidLocationException;
import org.ensembl.datamodel.Location;
import org.ensembl.datamodel.impl.DnaProteinAlignmentImpl;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.ConfigurationException;
import org.ensembl.driver.DnaProteinAlignmentAdaptor;
import org.ensembl.util.NotImplementedYetException;

// TODO Implement fetch/store/delete for new schema/API

public class
      MySQLDnaProteinAlignmentAdaptor
extends
      MySQLBaseFeatureAdaptor
implements
      DnaProteinAlignmentAdaptor
{
  private static final Logger logger = Logger.getLogger(MySQLDnaProteinAlignmentAdaptor.class.getName());


  public MySQLDnaProteinAlignmentAdaptor(MySQLDriver driver) {
    super(driver, TYPE);
  }

  public MySQLDnaProteinAlignmentAdaptor(MySQLDriver driver, String logicName, String type) {
    super(driver, logicName, type);
  }

  public MySQLDnaProteinAlignmentAdaptor(MySQLDriver driver, String[] logicNames, String type) {
    super(driver, logicNames, type);
  }

  void configure() throws ConfigurationException, AdaptorException {
    super.configure();
  }

  public String[] columns() {
    return new String[]{
             "protein_align_feature.protein_align_feature_id",
             "protein_align_feature.seq_region_id",
             "protein_align_feature.seq_region_start",
             "protein_align_feature.seq_region_end",
             "protein_align_feature.seq_region_strand",
             "protein_align_feature.hit_start",
             "protein_align_feature.hit_end",
             "protein_align_feature.hit_name",
             "protein_align_feature.analysis_id",
             "protein_align_feature.score",
             "protein_align_feature.evalue",
             "protein_align_feature.perc_ident",
             "protein_align_feature.cigar_line"
           };
  } // columns

  public String[][] tables() {
    return new String[][] {{"protein_align_feature","protein_align_feature"}};
  }

  public Object createObject(ResultSet resultSet) throws AdaptorException {
    DnaProteinAlignmentImpl align = null;
    Location location;
    CoordinateSystem myCoordinateSystem;

    try {
      if (resultSet.next()) {
        align = new DnaProteinAlignmentImpl(getDriver());
        align.setInternalID(resultSet.getLong("protein_align_feature_id"));
        location =
          locationConverter.idToLocation(
            resultSet.getLong("protein_align_feature.seq_region_id"),
            resultSet.getInt("protein_align_feature.seq_region_start"),
            resultSet.getInt("protein_align_feature.seq_region_end"),
            resultSet.getInt("protein_align_feature.seq_region_strand")
          );

        myCoordinateSystem = location.getCoordinateSystem();

        align.setLocation(location);
        align.setScore(resultSet.getDouble("protein_align_feature.score"));
        align.setAnalysis(
          getAnalysisAdaptor().fetch(resultSet.getLong("protein_align_feature.analysis_id"))
        );

        align.setHitAccession(resultSet.getString("protein_align_feature.hit_name"));
        align.setHitDisplayName(align.getHitAccession());
        align.setHitDescription(align.getHitAccession());
        align.setHitLocation(
          new Location(
            myCoordinateSystem, 
            align.getHitAccession(),
            resultSet.getInt("protein_align_feature.hit_start"),
            resultSet.getInt("protein_align_feature.hit_end"),
            0
          )
        );

        align.setPercentageIdentity(resultSet.getDouble("protein_align_feature.perc_ident"));
        align.setEvalue(resultSet.getDouble("protein_align_feature.evalue"));
        align.setScore(resultSet.getDouble("protein_align_feature.score"));
        align.setCigarString(resultSet.getString("protein_align_feature.cigar_line"));
        align.setDriver(getDriver());
      }

    } catch (InvalidLocationException exception) {
      throw new AdaptorException("Error when building Location", exception);
    } catch (SQLException exception) {
      throw new AdaptorException("SQL error when building object", exception);
    }

    return align;
  }

  /**
   * @param internalID The id of the DnaProteinAlignment to be fetched
   * @return A feature matching the internalID, or null if non found.
   */
  public DnaProteinAlignment fetch(long internalID) throws AdaptorException {

    return (DnaProteinAlignment) fetchByInternalID(internalID);

  }

  /**
   * Warning: this is potentially a very slow query.
   * @return list of DnaProteinAlignments corresponding to the specified analysis
   */
  public List fetch(String logicalName) throws AdaptorException {
    return fetchByNonLocationConstraint( getAnalysisIDCondition(logicalName) );
  }
  
  public Iterator fetchIterator(String logicalName) throws AdaptorException {
  	return fetchIteratorBySQL( "SELECT protein_align_feature_id " +
  			"FROM protein_align_feature paf, analysis a " +
  			"WHERE paf.analysis_id=a.analysis_id AND logic_name=\""+logicalName+"\"" );
  }
  
  public long store(DnaProteinAlignment feature) throws AdaptorException {
    throw new NotImplementedYetException("Not yet implemented for new API");
  }

  public void delete(DnaProteinAlignment feature) throws AdaptorException {

    throw new NotImplementedYetException("Not yet implemented for new API");

  }

  public void delete(long internalID) throws AdaptorException {

    throw new NotImplementedYetException("Not yet implemented for new API");

  }

  
}
