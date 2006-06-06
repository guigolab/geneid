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

import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.ensembl.datamodel.Gene;
import org.ensembl.datamodel.Transcript;

/**
 * Build gene scores based on transcript scores.
 */
public class GeneScoreBuilder extends ScoreBuilder {

    public GeneScoreBuilder(Cache cache) {

        super(cache);

    }

    public GeneScoreBuilder() {

        super();

    }

    //  -------------------------------------------------------------------------
    /**
     * Build a matrix of source/target gene scores based on transcript scores.
     */
    public ScoredMappingMatrix buildScoringMatrix(ScoredMappingMatrix transcriptScores) {

        // find which gene pairs actually score, store result in a "flag" matrix
        // this avoids having to do all-vs-all comparison
        ScoredMappingMatrix geneFlagMatrix = findScoringGenes(cache.getSourceGenesByInternalID().values(), cache
                .getTargetGenesByInternalID().values(), transcriptScores);

        ScoredMappingMatrix geneScores = buildMatrixFromGenes(cache.getSourceGenesByInternalID(), cache
                .getTargetGenesByInternalID(), geneFlagMatrix, transcriptScores);

        return geneScores;

    }

    // -------------------------------------------------------------------------
    /**
     * Score source genes with target genes.
     *  
     */
    private ScoredMappingMatrix buildMatrixFromGenes(Map sourceGenes, Map targetGenes, ScoredMappingMatrix geneFlagMatrix,
            ScoredMappingMatrix transcriptScores) {

        //System.out.println("# Building scoring matrix from genes and 'flag' matrix");
        ScoredMappingMatrix geneScores = new ScoredMappingMatrix();

        // iterate over ALL source genes, but only score against target genes that are in the flag matrix
        Iterator sit = sourceGenes.keySet().iterator();
        int i = 0;
        while (sit.hasNext()) {

            long sourceGeneID = ((Long) sit.next()).longValue();
            Gene sourceGene = (Gene) sourceGenes.get(new Long(sourceGeneID));

            List sourceTranscripts = sourceGene.getTranscripts();

            // find which entries in the flag matrix have this source
            // gene as the source
            List scoringTargetGeneEntries = geneFlagMatrix.sourceEntries(sourceGeneID);
            Iterator sgit = scoringTargetGeneEntries.iterator();

            // compare this source gene with each scoring target gene
            while (sgit.hasNext()) {

                Entry e = (Entry) sgit.next();
                Gene targetGene = (Gene) targetGenes.get(new Long(e.getTarget()));

                long targetGeneID = targetGene.getInternalID();
                float sourceGeneScore = 0.0f;
                float targetGeneScore = 0.0f;
                long sourceGeneLength = 0; // will be cumulative length of all transcripts
                long targetGeneLength = 0; // will be cumulative length of all transcripts

                // We are only interested in scoring with transcripts that are in the target gene.
                // The sourceTranscripts scored mapping matrix may contain scores for transcripts
                // that aren't in this gene so store a map of the target gene's transcripts

                Map targetGeneTranscriptMap = new HashMap();
                Iterator tgtit = targetGene.getTranscripts().iterator();
                while (tgtit.hasNext()) {
                    Transcript t = (Transcript) tgtit.next();
                    targetGeneTranscriptMap.put(new Long(t.getInternalID()), t);
                }

                // now find the highest scoring entry that has a target transcript belonging to the
                // target gene (hence second argument to findHighestScoreEntry())
                Iterator seit = sourceTranscripts.iterator();
                while (seit.hasNext()) {

                    Transcript sourceTranscript = (Transcript) seit.next();
                    List targetEntries = transcriptScores.sourceEntries(sourceTranscript.getInternalID());
                    // not all source exons will have corresponding targets
                    if (targetEntries.size() > 0) {
                        Entry highestScoreEntry = findHighestScoreTargetEntry(targetEntries, targetGeneTranscriptMap, false);
                        if (highestScoreEntry != null) { // may be none
                            // note transcript length here is combined length of all transcript's exons
                            sourceGeneScore += highestScoreEntry.getScore() * sourceTranscript.getLength();
                            sourceGeneLength += sourceTranscript.getLength();
                        }
                    }

                }

                // calculate target gene score similarly for target transcripts
                Map sourceGeneTranscriptMap = new HashMap();
                Iterator steit = sourceGene.getTranscripts().iterator();
                while (steit.hasNext()) {
                    Transcript t = (Transcript) steit.next();
                    sourceGeneTranscriptMap.put(new Long(t.getInternalID()), t);
                }
                Iterator teit = targetGene.getTranscripts().iterator();
                while (teit.hasNext()) {

                    Transcript targetTranscript = (Transcript) teit.next();
                    List sourceEntries = transcriptScores.targetEntries(targetTranscript.getInternalID());
                    // not all target exons will have corresponding targets
                    if (sourceEntries.size() > 0) {
                        Entry highestScoreEntry = findHighestScoreSourceEntry(sourceEntries, sourceGeneTranscriptMap, false);
                        if (highestScoreEntry != null) { // may be none
                            // note transcript length here is combined length of all transcript's exons
                            targetGeneScore += highestScoreEntry.getScore() * targetTranscript.getLength();
                            targetGeneLength += targetTranscript.getLength();
                        }
                    }

                }

                // calculate gene score and store in geneScores
                float geneScore = 0.0f;
                if (sourceGeneLength + targetGeneLength > 0) {

                    geneScore = (sourceGeneScore + targetGeneScore) / (sourceGeneLength + targetGeneLength);

                } else {

                    System.err.println("Error: Combined lengths of source gene " + sourceGeneID + " and target gene "
                            + targetGeneID + " are zero!");
                }

                geneScores.addScore(sourceGeneID, targetGeneID, geneScore);

            } // while target gene

            i++;

        } // while source gene

        return geneScores;

    }

    //---------------------------------------------------------------------
    /**
     * Find which source and target genes actually have scoring transcripts.
     * 
     * @return A ScoredMappingMatrix with entries only for scoring source/target gene pairs.
     */
    public ScoredMappingMatrix findScoringGenes(Collection allSourceGenes, Collection allTargetGenes,
            ScoredMappingMatrix transcriptScores) {

        // build target hashmap of genes keyed on transcript internal ID (may be
        // >1 gene per transcript, so this is a map of lists)
        Map targetGeneByTranscript = cache.getTargetGeneByTranscriptInternalID();

        // build "flag" score matrix
        ScoredMappingMatrix result = new ScoredMappingMatrix();

        Iterator astIt = allSourceGenes.iterator();
        while (astIt.hasNext()) {
            Gene sourceGene = (Gene) astIt.next();

            Iterator tIt = sourceGene.getTranscripts().iterator();
            while (tIt.hasNext()) {

                Transcript t = (Transcript) tIt.next();
                // find all scoring target transcripts
                List scoringEntries = transcriptScores.sourceEntries(t.getInternalID());
                Iterator scoringEntriesIterator = scoringEntries.iterator();
                while (scoringEntriesIterator.hasNext()) {

                    Entry scoringEntry = (Entry) scoringEntriesIterator.next();
                    // use target hashmap to get transcripts from target exons
                    Gene targetGene = (Gene) targetGeneByTranscript.get(new Long(scoringEntry.getTarget()));

                    if (targetGene != null) {

                        result.addScore(sourceGene.getInternalID(), targetGene.getInternalID(), 1.0f);

                    }

                }

            }
        }

        return result;

    }

    // -------------------------------------------------------------------------

}