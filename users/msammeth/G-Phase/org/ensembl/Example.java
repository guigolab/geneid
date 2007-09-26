/*
 * Created on Mar 3, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package org.ensembl;

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


import java.text.ParseException;
import java.util.Iterator;
import java.util.List;

import org.ensembl.datamodel.CoordinateSystem;
import org.ensembl.datamodel.Exon;
import org.ensembl.datamodel.ExternalDatabase;
import org.ensembl.datamodel.ExternalRef;
import org.ensembl.datamodel.Gene;
import org.ensembl.datamodel.GeneSnapShot;
import org.ensembl.datamodel.KaryotypeBand;
import org.ensembl.datamodel.Location;
import org.ensembl.datamodel.SequenceRegion;
import org.ensembl.datamodel.Transcript;
import org.ensembl.datamodel.TranscriptSnapShot;
import org.ensembl.datamodel.Translation;
import org.ensembl.datamodel.TranslationSnapShot;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.ConfigurationException;
import org.ensembl.driver.Driver;
import org.ensembl.driver.DriverManager;
import org.ensembl.driver.ExonAdaptor;
import org.ensembl.driver.GeneAdaptor;
import org.ensembl.driver.StableIDEventAdaptor;
import org.ensembl.util.SystemUtil;

/**
 * Example code illustrating how to use various parts of ensj-core. It shows
 * how to retrieve genes and exons from the public mysql database on 
 * ensembldb.ensembl.org
 *
 * <p>The source has been split into logicaly distinct sections which focus
 * on different parts ensj-core: initialising a driver, fetching data items
 * by internal id, creating different types of locations, using locations to
 * retrieve genes and accessing exons that belong to a gene.
 * 
 * <p>You can retrieve other data types such as transcripts and
 * translations using equivalant mechanisms.
 *
 * <p>
 * <dl>
 * 
 * <dt>Basic Usage: connecting to default database specified in a 
 * file in the distributojn. If this fails
 * it might be because the database specified in the file
 * has been removed from the database server. 
 * In this case try the advanced option below specifying an 
 * existing database.</dt>
 * <dd><code>java org.ensembl.Example</code></dd>
 * 
 * <dt>Advanced Usage: connecting to user specified database. See
 * createDriver(String) for details about the file format.</dt>
 * <dd><code>java org.ensembl.Example SOME_PROPERTIES_FILE</code></dd>
 * 
 * </dl>
 * @author Craig Melsopp 
 * @see <a href="Example.java">Example.java source</a>
 * @see #createDriver(String)
 *  
 */

public class Example {

  /**
   * Runs the Example program.
   * @param args is args contains one or more values the first
   * is used as the database parameters file.
   * @throws AdaptorException if a problem occurs retrieving data
   * @throws ConfigurationException if a problem occurs initialising the
   * driver 
   * @throws ParseException if a problem occurs parsing the string used to
   * construct a location
   */
  public static void main(String[] args)
    throws ConfigurationException, AdaptorException, ParseException {

  	// Dumps some key state information about the runtime environment-
  	// useful for debugging purposes
  	System.out.println(" *** RUNTIME CONFIGURATION (useful for debugging) *** ");
//  	System.out.println(SystemUtil.environmentDump());
		
		// default configuration file
		String configFilename = "resources/data/example_database.properties";
		if (args.length>0) 
    		configFilename = args[0];	
      
    Driver driver = createDriver(configFilename);

    System.out.println("\n\n\n *** ENSJ TEST OUTPUT *** ");
    
    displayDriverState(driver);

    fetchAnExonByInternalID(driver, 5639);

    Location[] locations = createLocations();

    countGenesAndExonsInEachLocation(driver, locations);

    fetchGeneByStableIDAndViewPeptide(driver);

    // create coordinate systems for later use when constructing locations    
    CoordinateSystem chromosomeCS = new CoordinateSystem("chromosome");
    CoordinateSystem cloneCS = new CoordinateSystem("clone");
    CoordinateSystem contigCS = new CoordinateSystem("contig");

    // create locations that we can use for fetching and converting later
    Location contigLoc = new Location(contigCS, "AL159978.14.1.206442");
    Location cloneLoc = new Location(cloneCS, "AB000878.1");
    Location chromosomeLoc = new Location("chromosome:22:21m-21.2m");
    // string 'convenince' constructor

    //  fetch information about "chromsomes"
    fetchSequenceRegionsSuchAsChromosomeOrContig(driver, chromosomeCS);

    //  fetch information about "contigs"
    fetchSequenceRegionsSuchAsChromosomeOrContig(driver, contigCS);

    fetchGenesByLocation(driver, cloneLoc);

    fetchGenesByLocation(driver, contigLoc);

    convertLocationToCoordinateSystemAndGetTheSeqRegionNames(
      driver,
      contigLoc,
      chromosomeCS);

    convertLocationToCoordinateSystemAndGetTheSeqRegionNames(
      driver,
      chromosomeLoc,
      contigCS);

    fetchDeletedGeneFromArchive(driver, "ENSG00000178007", 1);

    fetchKaryotypes(driver, chromosomeCS, "1");

    showExternalRefsForAGene(driver, "ENSG00000169861");

  }

  /**
   * Converts location into the target coordinate system and prints the
   * results.
   * @param driver driver to use
   * @param sourceLoc source location to be converted
   * @param targetCS target coordinate system to convert location into.
   * @throws AdaptorException if a problem occured during the conversion
   */
  public static void convertLocationToCoordinateSystemAndGetTheSeqRegionNames(
    Driver driver,
    Location sourceLoc,
    CoordinateSystem targetCS)
    throws AdaptorException {

    Location targetLoc =
      driver.getLocationConverter().convert(sourceLoc, targetCS);

    System.out.println("Source Location = " + sourceLoc);
    System.out.println("Target Location = " + targetLoc);

    // this is how we find out what seq regions our sourceLocation mapped
    // to in the target coordinate system. We need a loop because the location
    // maybe multi part.    
    for (Location node = targetLoc; node != null; node = node.next())
      System.out.println(
        "Target Location Sequence Region = " + node.getSeqRegionName());

    System.out.println();
  }

  /**
   * Loads the driver specified in the configFilename.
   * The file should contain at least these entries:
   * <pre>
   * host=ensembldb.ensembl.org
   * user=anonymous
   * database=SOME DATABASE
   * </pre>
   * @param configFilename name of file containing database parameters.
   * @return driver ready to use.
   * @throws ConfigurationException if a problem occured during the driver
   * initilisation.
   */
  public static Driver createDriver(String configFilename)
    throws ConfigurationException {

    // Load a driver for the latest human release.
    Driver driver = DriverManager.loadDriver(configFilename);

    return driver;

  }

  /**
   * Fetch an exon by it's internal id and print it.
   * @param driver driver to use to retrieve exon
   * @param exonInternalID id of exon to retrieve
   * @throws AdaptorException if a problem occured during the retrieval
   */
  public static void fetchAnExonByInternalID(Driver driver, int exonInternalID)
    throws AdaptorException {

    ExonAdaptor exonAdaptor = driver.getExonAdaptor();

    // Fetch an exon based on it's internalID. The same approach can be used
    // for all the adaptors which support fetch( internalID ). 
    Exon exon = exonAdaptor.fetch(exonInternalID);

    // All of the org.datamodel.impl classes support java's toString()
    // method. This means that can print them to find out their current
    // state.
    System.out.println(
      "exon with internal id " + exonInternalID + " = " + exon + "\n");

  }

  /**
   * Create and then return some ensj locations.
   * @return array of locations that can used in queries.
   */
  public static Location[] createLocations() {

    // Ensj provides support for using genomic locations. These are used to
    // represent the genomic of most of the biological datatypes, via the
    // getLocation() method and in addition can be used to specify database
    // queries. 

    // Every location must have a co-ordinate system; this is defined by
    // the CoordinateSystem object passed to the Location on creation.
    Location[] locations = new Location[1];

    // Create an assembly location. Assembly locations represent part of a
    // genome assembly, in this case part of chromosome 12.
    locations[0] = new Location(new CoordinateSystem("chromosome"), "12",
      // chromosome name
    1, // start
    100000, // end
  1); // strand

    return locations;
  }

  /**
   * For each of the locations print the number of genes and the number of
   * exons.
   * @param driver driver to execute queries against
   * @param locations locations to use in queries
   * @throws AdaptorException if problem occurs during retrieval
   */
  public static void countGenesAndExonsInEachLocation(
    Driver driver,
    Location[] locations)
    throws AdaptorException {

    // The easiest way to get a handle on an adaptor is if you already have
    // it's parent driver.
    GeneAdaptor geneAdaptor = driver.getGeneAdaptor();

    // Count the number of genes and exons in each of the locations.
    for (int i = 0; i < locations.length; ++i) {

      System.out.println("Location = " + locations[i]);

      List genes = geneAdaptor.fetch(locations[i]);

      if (genes == null || genes.size() == 0) {
        System.out.println("No Genes found.");
      } else {
        int geneCount = 0;
        int exonCount = 0;
        Iterator iter = genes.iterator();
        while (iter.hasNext()) {
          Gene gene = (Gene) iter.next();
          geneCount++;
          exonCount += gene.getExons().size();
        }

        System.out.println("num genes = " + geneCount);
        System.out.println("num exons = " + exonCount);

      }

      System.out.println(); // blank line to split result sections
    }
  }

  /**
   * Fetch a gene by it's stable ID and then print the peptide corresponding
   * to its first transcript.
   * @param driver driver to execute queries against
   * @throws AdaptorException if problem occurs during retrieval
   */
  public static void fetchGeneByStableIDAndViewPeptide(Driver driver)
    throws AdaptorException {

    GeneAdaptor geneAdaptor = driver.getGeneAdaptor();

    Gene gene = geneAdaptor.fetch("ENSG00000179902");

    Transcript transcript = (Transcript) gene.getTranscripts().get(0);
    Translation translation = transcript.getTranslation();
    String peptide = translation.getPeptide();

    System.out.println(
      "Peptide for " + translation.getAccessionID() + " : " + peptide);

  }

  /**
   * This method illustrates several ways to retrieve
   * genes from a specified location.
   * 
   * With locations containing few genes the memory and speed
   * differences between the methods will be small but for
   * locations containing many genes the difference
   * is potentially huge.
   * 
   * 
   * @param driver driver to execute queries against
   * @param location location to fetch genes from
   * @throws AdaptorException if problem occurs during retrieval
   */
  public static void fetchGenesByLocation(Driver driver, Location location)
    throws AdaptorException {

    int nExons = 0;
    int nGenes = 0;

    // use the adaptor directly from the driver
    // to load all of the genes into memory in one go.
    // The transcripts, translations and exons are 
    // lazy loaded on demand many lazy load requests require
    // a separate database access.
    List genes = driver.getGeneAdaptor().fetch(location);

    // load all of the genes with their child transcripts, translations
    // and exons preloaded. This will often provide faster
    // access to the child data than lazy loading it because 
    // it requires fewer database accesses.
    List genesWithChildren = driver.getGeneAdaptor().fetch(location, true);

    // iterating over the genes provides a compromise between
    // loading all of the genes with children (fastest + largest memory usage)
    // and loading the genes one at a time and lazy loading their
    // children (slowest + minumum memory requirement). In this
    // case we also preload the child data. Iterators are fairly
    // eficient in terms of both speed and memory usage and is
    // very useful for large datasets which are too big to fit in
    // memory.
    Iterator geneIterator =
      driver.getGeneAdaptor().fetchIterator(location, true);

    // make sure the exons are loaded and report the numbers
    // of genes and exons loaded. These should be the same for
    // list/iterator.

    nExons = 0;
    nGenes = genes.size();
    for (int i = 0, n = genes.size(); i < n; i++) {
      Gene g = (Gene) genes.get(i);
      nExons += g.getExons().size();
    }
    System.out.println(
      location.toString()
        + " has "
        + nGenes
        + " genes and "
        + nExons
        + " exons.");

    nExons = 0;
    nGenes = genesWithChildren.size();
    for (int i = 0, n = genes.size(); i < n; i++) {
      Gene g = (Gene) genes.get(i);
      nExons += g.getExons().size();
    }
    System.out.println(
      location.toString()
        + " has "
        + nGenes
        + " genes and "
        + nExons
        + " exons.");

    nGenes = 0;
    nExons = 0;
    while (geneIterator.hasNext()) {
      nGenes++;
      nExons += ((Gene) geneIterator.next()).getExons().size();
    }

    System.out.println(
      location.toString()
        + " has "
        + nGenes
        + " genes and "
        + nExons
        + " exons.");
    System.out.println(); // blank line to split result sections

  }

  /**
   * Fetches information about a gene from the archive.
   * @param driver driver to get get data from
   * @param geneStableID deleted stable ID
   * @param geneVersion version of gene
   * @throws AdaptorException if problem occurs during retrieval
   */
  public static void fetchDeletedGeneFromArchive(
    Driver driver,
    String geneStableID,
    int geneVersion)
    throws AdaptorException {

    StableIDEventAdaptor adaptor = driver.getStableIDEventAdaptor();

    // Find stableIDs in the current release that relate to the geneStableID
    List relatedIDs = adaptor.fetchCurrent(geneStableID);
    for (Iterator iter = relatedIDs.iterator(); iter.hasNext();) {
      String relatedID = (String) iter.next();
      System.out.println(
        geneStableID
          + " is related to "
          + relatedID
          + " in the current release.");
    }

    // This section requires schema version >= 15 which are currently in development only. 

    // Find the snapshot of the Gene's structure when it changed or was deleted.
    GeneSnapShot geneSnapshot =
      adaptor.fetchGeneSnapShot(geneStableID, geneVersion);
    // we already have the ID and version but this shows how to get them from the snapshot
    String gStableID = geneSnapshot.getArchiveStableID().getStableID();
    String gVersion = geneSnapshot.getArchiveStableID().getStableID();

    TranscriptSnapShot[] transcriptSnapShots =
      geneSnapshot.getTranscriptSnapShots();
    for (int i = 0; i < transcriptSnapShots.length; i++) {

      TranscriptSnapShot tSnapShot = transcriptSnapShots[i];
      String tStableID = tSnapShot.getArchiveStableID().getStableID();
      int tVersion = tSnapShot.getArchiveStableID().getVersion();

      TranslationSnapShot tnSnapShot = tSnapShot.getTranslationSnapShot();
      String tnStableID = tnSnapShot.getArchiveStableID().getStableID();
      int tnVersion = tnSnapShot.getArchiveStableID().getVersion();

      System.out.println(
        gStableID
          + "."
          + gVersion
          + " -> "
          + tStableID
          + "."
          + tVersion
          + " -> "
          + tnStableID
          + "."
          + tnVersion);

      // If there is a peptide associated with the translation print it too.
      String peptide = tnSnapShot.getPeptide();
      if (peptide != null)
        System.out.println("Peptide: " + peptide);
    }

    System.out.println();

  }

  /**
   * Fetches information about karyotypes (chromosome bands) for the
   * specified chromosome.
   * @param driver driver to get get data from
   * @param coordinateSystem coordinate system containing the chromosomeName
   * @param chromosomeName name of the chromosome of interest
   * @throws AdaptorException if problem occurs during retrieval
   */
  public static void fetchKaryotypes(
    Driver driver,
    CoordinateSystem coordinateSystem,
    String chromosomeName)
    throws AdaptorException {

    List l =
      driver.getKaryotypeBandAdaptor().fetch(coordinateSystem, chromosomeName);

    System.out.println(
      "Chromosome " + chromosomeName + " has " + l.size() + " karyotypes.");

    KaryotypeBand kb = (KaryotypeBand) l.get(0);
    Location loc = kb.getLocation();
    int start = loc.getStart();
    int end = loc.getEnd();

    System.out.println(
      "The first karyotype on chromosome "
        + chromosomeName
        + " is from "
        + start
        + "bp to "
        + end
        + "bp.");

    System.out.println();

  }

  /**
   * Fetches information about the sequence regions.
   * @param driver driver to get get data from
   * @throws AdaptorException if problem occurs during retrieval
   */
  public static void fetchSequenceRegionsSuchAsChromosomeOrContig(
    Driver driver,
    CoordinateSystem coordinateSystem)
    throws AdaptorException {

    // What sequence regions are are available?
    // In the case of the "chromosome" coordinate system these are chromosomes
    SequenceRegion[] seqRegions =
      driver.getSequenceRegionAdaptor().fetchAllByCoordinateSystem(
        coordinateSystem);
    System.out.println(
      "There are "
        + seqRegions.length
        + " sequence regions in the "
        + coordinateSystem.getName()
        + "."
        + coordinateSystem.getVersion()
        + " coordinate system.");

    SequenceRegion sr = seqRegions[0];
    System.out.println(
      coordinateSystem.getName()
        + " "
        + sr.getName()
        + " has length "
        + sr.getLength());

    System.out.println();
  }

  /**
   * Fetches the specified gene from the database and 
   * , if available, prints a summary of it's external refs.
   * 
   * @param driver driver to get get data from
   * @param geneAccession accession (ensembl stable id) of the gene to retrieve 
   * @throws AdaptorException
   */
  public static void showExternalRefsForAGene(Driver driver, String geneAccession)
    throws AdaptorException {
    Gene gene = driver.getGeneAdaptor().fetch(geneAccession);
    List xrefs = gene.getExternalRefs();
    if (xrefs.size() == 0) {
      System.out.println("No xrefs for gene" + geneAccession);
    } else {

      for (int i = 0, n = xrefs.size(); i < n; i++) {
        ExternalRef xref = (ExternalRef) xrefs.get(i);
        ExternalDatabase xdb = xref.getExternalDatabase();
        System.out.println(
          geneAccession
            + " has xref "
            + xref.getDisplayID()
            + " in "
            + xdb.getName()
            + "."
            + xdb.getVersion());
      }
    }
    
  }

  
  /**
   * Print details about the driver's configuration.
   * @param driver driver of interest.
   * @throws AdaptorException
   */
  public static void displayDriverState(Driver driver) throws AdaptorException{
    System.out.println("Driver connection: " + driver.toString());
    System.out.println("Driver configuration: " +driver.getConfiguration());
    System.out.println();
  }

} // Example

