package org.ensembl.driver.plugin.standard;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.List;
import java.util.logging.Logger;

import org.ensembl.datamodel.CoordinateSystem;
import org.ensembl.datamodel.DnaDnaAlignment;
import org.ensembl.datamodel.InvalidLocationException;
import org.ensembl.datamodel.Location;
import org.ensembl.datamodel.impl.DnaDnaAlignmentImpl;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.DnaDnaAlignmentAdaptor;
import org.ensembl.util.NotImplementedYetException;

public class
      MySQLDnaDnaAlignmentAdaptor
      extends
      MySQLBaseFeatureAdaptor
      implements
      DnaDnaAlignmentAdaptor
{

  private static final Logger logger = Logger.getLogger(MySQLDnaDnaAlignmentAdaptor.class.getName());


  /**
   * Constructs an adaptor capable of retrieving all DnaDnaAlignments from
   * the db table "dna_align_feature".
   */
  public MySQLDnaDnaAlignmentAdaptor(MySQLDriver driver) {
    super(driver, TYPE);
  }

  /**
   * Constructs an adaptor capable of retrieving all DnaDnaAlignments from
   * the db table "dna_align_feature" with the specified analysis logicName.
   */
  public MySQLDnaDnaAlignmentAdaptor(MySQLDriver driver, String logicName, String type) {
    super(driver, logicName, type);
  }

  /**
   * Constructs an adaptor capable of retrieving all DnaDnaAlignments from
   * the db table "dna_align_feature" with the specified analysis logicNames.
   */
  public MySQLDnaDnaAlignmentAdaptor(MySQLDriver driver, String[] logicNames, String type) {
    super(driver, logicNames, type);
  }

  public String[] columns() {

    return
      new String[]{
        "dna_align_feature.dna_align_feature_id",
        "dna_align_feature.seq_region_id",
        "dna_align_feature.seq_region_start",
        "dna_align_feature.seq_region_end",
        "dna_align_feature.seq_region_strand",
        "dna_align_feature.hit_start",
        "dna_align_feature.hit_end",
        "dna_align_feature.hit_strand",
        "dna_align_feature.hit_name",
        "dna_align_feature.analysis_id",
        "dna_align_feature.score",
        "dna_align_feature.evalue",
        "dna_align_feature.perc_ident",
        "dna_align_feature.cigar_line"
      };

  } // columns

  public String[][] tables() {
    return new String[][] {{"dna_align_feature","dna_align_feature"}};
  }

  public Object createObject(ResultSet resultSet) throws AdaptorException {

    DnaDnaAlignmentImpl align = null;
    Location location;
    CoordinateSystem myCoordinateSystem;

    try {
      if (resultSet.next()) {

        align = new DnaDnaAlignmentImpl(getDriver());

        align.setInternalID(resultSet.getLong("dna_align_feature.dna_align_feature_id"));
        
        location =
          locationConverter.idToLocation(
            resultSet.getLong("dna_align_feature.seq_region_id"),
            resultSet.getInt("dna_align_feature.seq_region_start"),
            resultSet.getInt("dna_align_feature.seq_region_end"),
            resultSet.getInt("dna_align_feature.seq_region_strand")
          );

        myCoordinateSystem = location.getCoordinateSystem();

        align.setLocation(location);
        align.setScore(resultSet.getDouble("dna_align_feature.score"));
        align.setAnalysis(
          getAnalysisAdaptor().fetch(resultSet.getLong("dna_align_feature.analysis_id"))
        );

        location =
          locationConverter.idToLocation(
            resultSet.getLong("dna_align_feature.seq_region_id"),
            resultSet.getInt("dna_align_feature.hit_start"),
            resultSet.getInt("dna_align_feature.hit_end"),
            resultSet.getInt("dna_align_feature.hit_strand")
          );

        align.setHitAccession(resultSet.getString("dna_align_feature.hit_name"));
        align.setHitDisplayName(align.getHitAccession());
        align.setHitDescription(align.getHitAccession());
        align.setHitLocation(
          new Location(
            myCoordinateSystem,
            align.getHitAccession(),
            resultSet.getInt("dna_align_feature.hit_start"),
            resultSet.getInt("dna_align_feature.hit_end"),
            resultSet.getInt("dna_align_feature.hit_strand")
          )
        );

        align.setPercentageIdentity(resultSet.getDouble("dna_align_feature.perc_ident"));
        align.setEvalue(resultSet.getDouble("dna_align_feature.evalue"));
        align.setScore(resultSet.getDouble("dna_align_feature.score"));
        align.setCigarString(resultSet.getString("dna_align_feature.cigar_line"));
        align.setDriver(getDriver());
      }

    } catch (InvalidLocationException exception) {
      throw new AdaptorException("Error when building Location", exception);
    } catch (SQLException exception) {
      throw new AdaptorException("SQL error when building object", exception);
    }

    return align;
  }

  public List fetch(String logicalName) throws AdaptorException {
    return fetchByNonLocationConstraint(getAnalysisIDCondition(logicalName));
  }
  
  public DnaDnaAlignment fetch(long internalID) throws AdaptorException {
    return (DnaDnaAlignment)super.fetchByInternalID(internalID);
  }

  public long store(DnaDnaAlignment feature) throws AdaptorException {
    throw new NotImplementedYetException("Not yet implemented for new API");
  }

  public void delete(DnaDnaAlignment feature) throws AdaptorException {
    throw new NotImplementedYetException("Not yet implemented for new API");
  }

  public void delete(long internalID) throws AdaptorException {
    throw new NotImplementedYetException("Not yet implemented for new API");
  }


}
