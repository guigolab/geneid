/*
    Copyright (C) 2001 EBI, GRL

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
package org.ensembl.driver.plugin.variation;

import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.ConfigurationException;
import org.ensembl.driver.Driver;
import org.ensembl.driver.plugin.standard.BaseDriver;
import org.ensembl.driver.variation.AlleleGroupAdaptor;
import org.ensembl.driver.variation.IndividualAdaptor;
import org.ensembl.driver.variation.IndividualGenotypeAdaptor;
import org.ensembl.driver.variation.LDFeatureAdaptor;
import org.ensembl.driver.variation.PopulationAdaptor;
import org.ensembl.driver.variation.PopulationGenotypeAdaptor;
import org.ensembl.driver.variation.TranscriptVariationAdaptor;
import org.ensembl.driver.variation.VariationAdaptor;
import org.ensembl.driver.variation.VariationDriver;
import org.ensembl.driver.variation.VariationFeatureAdaptor;
import org.ensembl.driver.variation.VariationGroupAdaptor;
import org.ensembl.driver.variation.VariationGroupFeatureAdaptor;

/**
 * This driver provides access to data in 
 * Ensembl variation databases.
 * 
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp</a>
 *
 */
public class MySQLVariationDriver
  extends BaseDriver
  implements VariationDriver {

  private VariationGroupAdaptor variationGroupAdaptor;

  private AlleleGroupAdaptor alleleGroupAdaptor;

  private PopulationAdaptor populationAdaptor;

  private VariationAdaptor variationAdaptor;

  private VariationFeatureAdaptor variationFeatureAdaptor;

  private Driver coreDriver;

  private IndividualAdaptor individualAdaptor;

private IndividualGenotypeAdaptor individualGenotypeAdaptor;

private PopulationGenotypeAdaptor populationGenotypeAdaptor;

private TranscriptVariationAdaptor transcriptVariationAdaptor;

private VariationGroupFeatureAdaptor variationGroupFeatureAdaptor;

private LDFeatureAdaptor lDFeatureAdaptor;

  /**
   * @throws AdaptorException
   */
  public MySQLVariationDriver() throws AdaptorException {
    super();
  }

  /**
   * @param host
   * @param database
   * @param user
   * @throws AdaptorException
   */
  public MySQLVariationDriver(String host, String database, String user)
    throws AdaptorException {
    super(host, database, user, false);
  }

  /**
   * @param host
   * @param database
   * @param user
   * @param password
   * @throws AdaptorException
   */
  public MySQLVariationDriver(
    String host,
    String database,
    String user,
    String password)
    throws AdaptorException {
    super(host, database, user, password, false);
  }

  /**
   * @param host
   * @param database
   * @param user
   * @param password
   * @param port
   * @throws AdaptorException
   */
  public MySQLVariationDriver(
    String host,
    String database,
    String user,
    String password,
    String port)
    throws AdaptorException {
    super(host, database, user, password, port, false);
  }

  /**
   * 
   * @return true if driver has been initiailised with a database
   * parameter.
   */
  private boolean initialisedForAdaptors() {
    return getConfiguration().getProperty("connection_string") != null;
  }

  /**
   * Does nothing because the adaptors are created on demand.
   */
  protected void loadAdaptors()
    throws AdaptorException, ConfigurationException {

  }

  /**
   * @see org.ensembl.driver.variation.VariationDriver#getVariationFeatureAdaptor()
   */
  public VariationFeatureAdaptor getVariationFeatureAdaptor()
    throws AdaptorException {

    if (variationFeatureAdaptor == null && initialisedForAdaptors())
      addAdaptor(
        variationFeatureAdaptor = new MySQLVariationFeatureAdaptor(this));

    return variationFeatureAdaptor;

  }

  /**
   * @see org.ensembl.driver.variation.VariationDriver#getVariationAdaptor()
   */
  public VariationAdaptor getVariationAdaptor() throws AdaptorException {

    if (variationAdaptor == null && initialisedForAdaptors())
      addAdaptor(variationAdaptor = new MySQLVariationAdaptor(this));

    return variationAdaptor;
  }

  /**
   * @see org.ensembl.driver.variation.VariationDriver#getPopulationAdaptor()
   */
  public PopulationAdaptor getPopulationAdaptor() throws AdaptorException {

    if (populationAdaptor == null && initialisedForAdaptors())
      addAdaptor(populationAdaptor = new MySQLPopulationAdaptor(this));

    return populationAdaptor;
  }


  /**
   * @see org.ensembl.driver.variation.VariationDriver#getVariationGroupAdaptor()
   */
  public VariationGroupAdaptor getVariationGroupAdaptor() throws AdaptorException {

    if (variationGroupAdaptor == null && initialisedForAdaptors())
      addAdaptor(variationGroupAdaptor = new MySQLVariationGroupAdaptor(this));

    return variationGroupAdaptor;
  }


  /**
   * @see org.ensembl.driver.variation.VariationDriver#getAlleleGroupAdaptor()
   */
  public AlleleGroupAdaptor getAlleleGroupAdaptor() throws AdaptorException {

    if (alleleGroupAdaptor == null && initialisedForAdaptors())
      addAdaptor(alleleGroupAdaptor = new MySQLAlleleGroupAdaptor(this));

    return alleleGroupAdaptor;
  }
  

  /**
   * @see org.ensembl.driver.variation.VariationDriver#getIndividualAdaptor()
   */
  public IndividualAdaptor getIndividualAdaptor() throws AdaptorException {

    if (individualAdaptor == null && initialisedForAdaptors())
      addAdaptor(individualAdaptor = new MySQLIndividualAdaptor(this));

    return individualAdaptor;
  }
  
  /**
   * @see org.ensembl.driver.variation.VariationDriver#getCoreDriver()
   */
  public Driver getCoreDriver() {
    return coreDriver;
  }

  /**
   * @see org.ensembl.driver.variation.VariationDriver#setCoreDriver(Driver coreDriver)
   */
  public void setCoreDriver(Driver coreDriver) throws AdaptorException {
    this.coreDriver = coreDriver;
    // the coordinate system adaptor in core needs to 
    // be able to see the meta_coord table in the variation
    // database to discover which coordinate systems 
    // variation related types are stored in.
    coreDriver.getCoordinateSystemAdaptor().addDataSource(getDataSource());
  }

/**
 * @see org.ensembl.driver.variation.VariationDriver#getIndividualGenotypeAdaptor()
 */
public IndividualGenotypeAdaptor getIndividualGenotypeAdaptor() throws AdaptorException {
	if (individualGenotypeAdaptor == null && initialisedForAdaptors())
    addAdaptor(individualGenotypeAdaptor = new MySQLIndividualGenotypeAdaptor(this));

	return individualGenotypeAdaptor;
}

/**
 * @see org.ensembl.driver.variation.VariationDriver#getPopulationGenotypeAdaptor()
 */
public PopulationGenotypeAdaptor getPopulationGenotypeAdaptor() throws AdaptorException {
	if (populationGenotypeAdaptor == null && initialisedForAdaptors())
    addAdaptor(populationGenotypeAdaptor = new MySQLPopulationGenotypeAdaptor(this));

	return populationGenotypeAdaptor;
}

/**
 * @see org.ensembl.driver.variation.VariationDriver#getTranscriptVariationAdaptor()
 */
public TranscriptVariationAdaptor getTranscriptVariationAdaptor() throws AdaptorException {
	if (transcriptVariationAdaptor == null && initialisedForAdaptors())
    addAdaptor(transcriptVariationAdaptor = new MySQLTranscriptVariationAdaptor(this));

	return transcriptVariationAdaptor;
}

/* (non-Javadoc)
 * @see org.ensembl.driver.variation.VariationDriver#getVariationGroupFeatureAdaptor()
 */
public VariationGroupFeatureAdaptor getVariationGroupFeatureAdaptor() throws AdaptorException {
  if (variationGroupFeatureAdaptor == null && initialisedForAdaptors())
    addAdaptor(variationGroupFeatureAdaptor = new MySQLVariationGroupFeatureAdaptor(this));

  return variationGroupFeatureAdaptor;
}

/**
 * @throws AdaptorException
 * @see org.ensembl.driver.variation.VariationDriver#getLDFeatureAdaptor()
 */
public LDFeatureAdaptor getLDFeatureAdaptor() throws AdaptorException {
	if (lDFeatureAdaptor == null && initialisedForAdaptors())
    addAdaptor(lDFeatureAdaptor = new MySQLLDFeatureAdaptor(this));

  return lDFeatureAdaptor;
}

}
