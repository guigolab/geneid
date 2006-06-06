/*
 * Copyright (C) 2004 EBI, GRL
 * 
 * This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License along with this library; if not, write to the Free
 * Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

package org.ensembl.idmapping;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.ensembl.datamodel.Gene;
import org.ensembl.datamodel.Transcript;
import org.ensembl.datamodel.Translation;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.TranslationAdaptor;

/**
 * Update the gene_archive and peptide_archive tables.
 */
public class Archiver {

    private Config conf;

    public Archiver(Config conf) {

        this.conf = conf;

    }

    // -------------------------------------------------------------------------
    /**
     * Add a row to the peptide archive for each existing translation that has had a version change.
     * 
     * i.e. deleted or upversioned,  but not new one.
     * 
     * @param changedTranslationEvents A List of StableIDEventContainer objects for the changed translations.
     */
    public void writeNewPeptideArchiveToFile(Cache cache, List changedTranslationEvents, String fileName) {

        // get the peptide sequences for each of the translations
        Map peptideSequences = new HashMap();
        try {

            TranslationAdaptor ta = conf.getSourceDriver().getTranslationAdaptor();
            Iterator it = changedTranslationEvents.iterator();
            while (it.hasNext()) {
                StableIDEventContainer event = (StableIDEventContainer) it.next();
                Translation t = (Translation) (cache.getSourceTranslationsByStableID().get(event.getOldStableID()));
                peptideSequences.put(t.getAccessionID(), t.getPeptide());
            }

        } catch (AdaptorException ae) {
            ae.printStackTrace();
        }

        // write to file
        try {

            OutputStreamWriter writer = new OutputStreamWriter(new FileOutputStream(fileName));
            Iterator it = changedTranslationEvents.iterator();
            while (it.hasNext()) {
                StableIDEventContainer event = (StableIDEventContainer) it.next();
                writer.write(event.getOldStableID() + "\t" + event.getOldVersion() + "\t"
                        + (String) peptideSequences.get(event.getOldStableID()) + "\n");
            }
            writer.close();

        } catch (IOException e) {

            e.printStackTrace();
        }

    }

    // -------------------------------------------------------------------------
    /**
     * Update the gene archive; for each gene where the version has changed, write a row for each translation.
     */
    public void writeNewGeneArchiveToFile(Cache cache, List changedGeneEvents, long mappingSessionID, String fileName) {

        try {

            OutputStreamWriter writer = new OutputStreamWriter(new FileOutputStream(fileName));

            Iterator it = changedGeneEvents.iterator();
            while (it.hasNext()) {

                StableIDEventContainer event = (StableIDEventContainer) it.next();
                Gene gene = (Gene) cache.getSourceGenesByStableID().get(event.getOldStableID());
                List transcripts = gene.getTranscripts();
                Iterator tit = transcripts.iterator();
                while (tit.hasNext()) {

                    Transcript transcript = (Transcript) tit.next();
                    Translation translation = transcript.getTranslation();
                    if (translation != null) {
                        writer.write(gene.getAccessionID() + "\t" + gene.getVersion() + "\t" + transcript.getAccessionID() + "\t"
                                + transcript.getVersion() + "\t" + translation.getAccessionID() + "\t" + translation.getVersion()
                                + "\t" + mappingSessionID + "\n");
                    }
                }

            }
            writer.close();

        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    // -------------------------------------------------------------------------
}