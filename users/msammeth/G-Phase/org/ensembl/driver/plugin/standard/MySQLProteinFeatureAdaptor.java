/*
  Copyright (C) 2003 EBI, GRL

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
package org.ensembl.driver.plugin.standard;


import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

import org.ensembl.datamodel.Analysis;
import org.ensembl.datamodel.ProteinFeature;
import org.ensembl.datamodel.Translation;
import org.ensembl.datamodel.impl.ProteinFeatureImpl;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.AnalysisAdaptor;
import org.ensembl.driver.ProteinFeatureAdaptor;
import org.ensembl.driver.TranslationAdaptor;




/**
 * Peptide feature.
 */
public class MySQLProteinFeatureAdaptor extends BaseAdaptor implements ProteinFeatureAdaptor {


  private static final Logger logger = Logger.getLogger( MySQLProteinFeatureAdaptor.class.getName() );
  public final boolean CLIP = false;
  private TranslationAdaptor translationAdaptor = null;
  private String analysisIDCondition = null;


  public MySQLProteinFeatureAdaptor( MySQLDriver driver ) {
    super( driver );
  }

  
   public String getType() {
     return TYPE;
   }


  /**
   * @return List of SimplePeptideFeature corresponding to the the specified translation.
   */
  public List fetch( Translation translation ) throws AdaptorException {
    
    List result = Collections.EMPTY_LIST;

    long translationInternalID = translation.getInternalID();


    String sql = "SELECT "
      + " f.protein_feature_id " //1
      + " ,f.translation_id " //2
      + " ,f.analysis_id " //3
      + " ,f.seq_start " //4
      + " ,f.seq_end " //5
      + " ,f.hit_start " //6
      + " ,f.hit_end " //7
      + " ,f.hit_id " //8
      + " ,f.score " //9
      + " ,f.evalue " // 10
      + " ,f.perc_ident " // 11
      + " FROM "
      + " protein_feature f "
      + " WHERE "
      + " f.translation_id=" + translationInternalID + " ";

    logger.fine( sql );

    Connection conn = null;
    try {
      conn = getConnection();
      ResultSet rs = conn.createStatement().executeQuery( sql );

      if ( rs.next() ) {

        result = new ArrayList();

        do {
          
          ProteinFeature spf = new ProteinFeatureImpl(driver);
          
          AnalysisAdaptor analysisAdaptor = driver.getAnalysisAdaptor();
          Analysis a = analysisAdaptor.fetch( rs.getLong( 3 ) );
          spf.setAnalysis( a );
          spf.setDisplayName( rs.getString(8) ); // hitID, no better name at the moment
          spf.setInternalID( rs.getLong(1) );
          spf.setPeptideStart( rs.getInt(4) );
          spf.setPeptideEnd( rs.getInt(5) );
          spf.setTranslationInternalID( rs.getLong(2) );
          spf.setTranslation( translation );
          // We ignore hit_start and hit_end
          spf.setScore( rs.getDouble(9) );
          spf.setEvalue( rs.getDouble(10) );
          spf.setPercentageIdentity( rs.getDouble(11) );

          spf.setDriver( driver );

          result.add( spf );
        } while ( rs.next() );
      }
    } catch ( Exception e ) {
      throw new AdaptorException("Failed to fetch SimplePeptides by translation" + sql, e);
    }finally {
      close( conn );
    }

    return result;
  }



  public ProteinFeature fetch( long internalID ) throws AdaptorException {

    String sql = "SELECT "
      + " f.protein_feature_id " //1
      + " ,f.translation_id " //2
      + " ,f.analysis_id " //3
      + " ,f.seq_start " //4
      + " ,f.seq_end " //5


      + " ,f.analysis_id " //6

      + " ,f.hit_start " //7
      + " ,f.hit_end " //8
      + " ,f.hit_id " //9

      + " ,f.score " //10
      + " ,f.evalue " //11
      + " ,f.perc_ident " //12

      + " FROM "
      + " protein_feature f "
      + " ,analysis ap "
      + " WHERE "
      + " f.protein_feature_id=" + internalID + " ";
    logger.fine( sql );

    Connection conn = null;
    try {
      conn = getConnection();
      ResultSet rs = conn.createStatement().executeQuery( sql );
      if ( !rs.next() ) {
        return null;
      } else {
        ProteinFeature spf = new ProteinFeatureImpl(driver);

        AnalysisAdaptor analysisAdaptor = driver.getAnalysisAdaptor();
        Analysis a = analysisAdaptor.fetch( rs.getLong( 3 ) );
        spf.setAnalysis( a );
        spf.setDisplayName( rs.getString(9) );
        spf.setInternalID( rs.getLong(1) );
        spf.setPeptideStart( rs.getInt(4) );
        spf.setPeptideEnd( rs.getInt(5) );

        final long translationID = rs.getLong(2);
        spf.setTranslationInternalID( translationID );
        // Note: this will be potentially slow unless translation adaptor
        // caches
        spf.setTranslation( driver.getTranslationAdaptor().fetch( translationID ) );
        
        spf.setScore( rs.getFloat(10) );
        spf.setEvalue( rs.getDouble(11) );
        spf.setPercentageIdentity( rs.getInt(12) );


        if ( translationAdaptor == null ) {
          translationAdaptor = ( TranslationAdaptor )getDriver().getAdaptor( "translation" );
        }
        spf.setTranslation(translationAdaptor.fetch( rs.getLong(2) ));

        spf.setDriver( driver );

        return spf;
      }
    } catch ( SQLException e ) {
      throw new AdaptorException("Rethrow + stacktrace", e);
    }finally {
      close( conn );
    }
  }


  /**
   * @return internalID assigned to feature in database.
   */
  public long store(ProteinFeature feature ) throws AdaptorException {
    // We ignore hit_start and hit_end
    String sql = "INSERT INTO protein_feature ("
      + " translation_id " // 1 
      + ",  seq_start " // 2
      + ",  seq_end " // 3
      + ",  analysis_id " // 4
      + ",  hit_id " // 5
      + ",  score " // 6
      + ",  evalue " // 7
      + ",  perc_ident " // 8
      + " ) VALUES (?, ?, ?, ?, ?, ?, ?, ? ) ";

    long internalID = 0;
    Connection conn = null;
    try {

      conn = getConnection();
      conn.setAutoCommit(false);   

      PreparedStatement ps = conn.prepareStatement(sql);
      ps.setLong( 1, feature.getTranslationInternalID() );
      ps.setInt( 2, feature.getPeptideStart() );
      ps.setInt( 3, feature.getPeptideEnd() );
      ps.setLong( 4, feature.getAnalysis().getInternalID() );
      ps.setString( 5, feature.getDisplayName() );
      ps.setDouble( 6, feature.getScore() );
      ps.setDouble( 7, feature.getEvalue() );
      ps.setDouble( 8, feature.getPercentageIdentity() );

      internalID = executeAutoInsert( ps, sql );

      conn.commit();
      feature.setDriver( driver );
      feature.setInternalID( internalID );
    } catch (Exception e) {
      rollback( conn );
      throw new AdaptorException("Failed to store SimplePeptideFeature: "
                                 + feature, e);
    }finally {
      close( conn );
    }

    addToCache( feature);

    return internalID;
  }



  /**
   * @param internalID internalID of feature to be deleted from database.
   */
  public void delete( long internalID ) throws AdaptorException {
      
    if ( internalID<1 ) return;

    deleteFromCache( internalID );

    Connection conn = null;
    try {

      conn = getConnection();
      conn.setAutoCommit(false); 

      delete( conn, internalID );
        
      conn.commit();
    } catch (Exception e) {
      rollback( conn );
      throw new AdaptorException("Failed to delete SimplePeptideFeature: "
                                 + internalID, e);
    }finally {
      close( conn );
    }

  }



  /**
   * @param feature feature to delete.
   */
  public void delete( ProteinFeature feature  ) throws AdaptorException {
    delete( feature.getInternalID() );
    feature.setInternalID(0);
  }



  /**
   * Executes sql to delete row from protein_feature table.
   */
  void delete(Connection conn, long internalID) throws AdaptorException {

    executeUpdate(conn, "delete from protein_feature where protein_feature_id = " + internalID);

  }
}
