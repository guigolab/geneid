package org.ensembl.driver.plugin.compara;

import java.sql.Connection;
import java.util.logging.Logger;

import org.ensembl.datamodel.compara.ComparaDataFactory;
import org.ensembl.datamodel.impl.compara.ComparaDataFactoryImpl;
import org.ensembl.driver.Adaptor;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.ConfigurationException;
import org.ensembl.driver.plugin.standard.MySQLQueryAdaptor;

/**
 * This driver behaves exactly like the "standard" core MySQLDriver, with the
 * exception that it deals with adapters and objects over the 
 * ensembl compara databases.
 * @see org.ensembl.driver.plugin.standard.MySQLDriver
 */
public class ComparaMySQLDriver
  extends org.ensembl.driver.plugin.standard.MySQLDriver {

  private static final Logger logger =
    Logger.getLogger(ComparaMySQLDriver.class.getName());
  protected ComparaDataFactory factory = null;

  public ComparaMySQLDriver() throws AdaptorException {
  }

  /**
   * Create adaptors and store reference to them to enable access via driver
   * in future.
   */
  protected void loadAdaptors()
    throws AdaptorException, ConfigurationException {

    factory = new ComparaDataFactoryImpl(this);

    Connection testConnection = getConnection();

    // Quit with warning if database not available. In this case the driver
    // will have be NO adaptors.
    try {
      testConnection.close();
    } catch (Exception e) {
      throw new AdaptorException("Failed to connect to database", e);
    }

    // query adaptor can be used without specifying a database
    addAdaptor(new MySQLQueryAdaptor(this));
    addAdaptor(new MySQLGenomeDBAdaptor(this));
    addAdaptor(new MySQLDnaFragmentAdaptor(this));
    addAdaptor(new MySQLGenomicAlignAdaptor(this));
    addAdaptor(new MySQLDnaDnaAlignFeatureAdaptor(this));
    addAdaptor(new MySQLHomologyAdaptor(this));
    addAdaptor(new MySQLMethodLinkAdaptor(this));

    // Configure adaptors loaded so far.
    Adaptor[] as = getAdaptors();
    for (int i = 0; i < as.length; i++) {
      Object o = as[i];
      if (o instanceof ComparaBaseAdaptor)
         ((ComparaBaseAdaptor) o).configure();

    } //end while

  } //end loadAdaptors

}
