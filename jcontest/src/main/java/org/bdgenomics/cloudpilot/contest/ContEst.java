package org.bdgenomics.cloudpilot.contest;

import org.apache.log4j.Logger;

import java.util.*;

public class ContEst {

    public static void main(String[] args) {

        int readLength = 100;
        double meanError = 0.001;
        int readDepth = 10;
        int steps = 100;
        int N = 100;

        Genome g = Genome.randomGenome(3, 10000);

        ArrayList<VariantSite> vs = g.randomVariantSites(200);

        Random rand = new Random();

        for(int i = 0; i < N; i++) {
            double c = rand.nextDouble() / 2.0;

            ArrayList<Read> reads = new ArrayList<Read>();
            for(VariantSite v : vs) {
                reads.addAll(v.sampleOverlappingReads(g, c, readLength, meanError, readDepth));
            }
            ContEst contest = new ContEst(vs, reads);

            double[] lls = contest.likelihoodGrid(steps);

            int argMax = argMax(lls);
            double spacing = 1.0 / (steps);
            double estimatedC = spacing * argMax;

            System.out.println(String.format("%.05f <- %.05f\t%.05f", c, estimatedC, (estimatedC - c)));
        }

    }

    private static Logger logger = Logger.getLogger(ContEst.class);

    private PositionTree<VariantSite> variantTree;
    private Collection<VariantSite> variants;
    private Collection<Read> reads;

    public ContEst(Collection<VariantSite> variants, Collection<Read> reads) {
        this.variants = variants;
        this.reads = reads;
        variantTree = PositionTree.createTree(variants);
    }

    /**
     * Given a collection of Reads which overlap a VariantSite, this function displays the
     * aligned reads so that the corresponding bases of each read line up.
     *
     * The display also highlights the portion of each read which overlaps the variant.
     *
     * @param vs The variant site that all the reads must overlap
     * @param rs The collection of reads overlapping the variant site.
     */
    public void printSiteAndReads(VariantSite vs, Collection<Read> rs) {
        int maxLeft = 0, maxRight = 0;
        for(Read r : rs) {
            maxLeft = Math.max(maxLeft, vs.offset() - r.alignment.offset());
            maxRight = Math.max(maxRight, r.alignment.offset() + r.sequence.length() - vs.offset() + 1);
        }

        for(Read r : rs) {
            int roffset = vs.offset() - r.alignment.offset();
            int leftPadding = maxLeft - roffset;
            String rseq = StringUtils.leftPad("", leftPadding, ' ') + StringUtils.upperCaseOffset(r.sequence, roffset);
            System.out.println(rseq);
        }

        System.out.println(String.format("%s* <- %s/%s (%.2f)",
                StringUtils.leftPad("", maxLeft + 1, ' '),
                vs.referenceAllele, vs.alternateAllele, vs.contaminatingPopulationFrequency));
    }

    /**
     * Finds the arg max_i (values[i]) -- choosing arbitrarily between ties.
     *
     * @param values The array of values in which to find the maximum
     * @return The index i whose value achieves the maximum in the array
     */
    public static int argMax(double[] values) {
        int maxId = -1;
        for(int i = 0; i < values.length; i++) {
            if(i == 0 || values[i] >= values[maxId]) {
                maxId = i;
            }
        }
        return maxId;
    }

    /**
     * Given an array of log-likelihoods, produces a corresponding array of log-likelihoods
     * that have been normalized (in NON-log space, i.e. ret[i] = ll[i] - log(sum(exp(ll[i]))) )
     *
     * @param logValues The array of log-likelihoods to be normalized
     * @return The normalized array of log-likelihoods
     */
    public static double[] logNormalize(double[] logValues) {
        double[] ns = new double[logValues.length];
        double logSum = 0.0;
        for(int i = 0; i < logValues.length; i++) {
            logSum = i == 0 ? logValues[i] :
                    Math.log(Math.exp(logSum) + Math.exp(logValues[i]));
        }

        for(int i = 0; i < logValues.length; i++) {
            ns[i] = logValues[i] - logSum;
        }
        return ns;
    }

    /**
     * Calculate a grid (array) of likelihoods, at equally-spaced points between [0.0, 1.0].
     *
     * @param steps The number of equally-spaced points in the likelihood grid.
     * @return The grid of likelihoods.
     */
    public double[] likelihoodGrid(int steps) {

        double[] lls = new double[steps];
        double step = 1.0 / (lls.length);

        int i = 0;

        for(double pc = 0.0; i < lls.length && pc <= 1.0; i++, pc += step) {
            double ll = logLikelihoodC(pc);
            lls[i] = ll;
        }

        return lls;
    }

    /**
     * Calculates, given the data (reads and variants) the likelihood of a given 'c' value -- that is,
     * the probability of generating the data given that the rate of contamination is 'c'.
     *
     * @param c The given rate of contamination
     * @return The likelihood of this rate of contamination
     */
    public double logLikelihoodC(double c) {
        Map<VariantSite, Double> siteLikelihoods = new HashMap<>();
        Map<VariantSite, Integer> siteDepths = new HashMap<>();
        Map<VariantSite, ArrayList<Read>> siteReads = new HashMap<>();

        for(VariantSite vs : variants) {
            siteReads.put(vs, new ArrayList<Read>());
            siteDepths.put(vs, 0);
        }

        /*
        For each read,
          1. find all the variant sites that overlap it
          2. calculate the likelihood P( c | ... ) given that read at that site
          3. multiply each of those likelihoods with the stored likelihoods for each site
         */
        for(Read r : reads) {
            for(VariantSite site : variantTree.findOverlapping(r.region(), new Accumulator.List<VariantSite>()).list()) {
                double likelihood = site.readLogLikelihood(c, r);
                if(!siteLikelihoods.containsKey(site)) {
                    siteLikelihoods.put(site, likelihood);
                } else {
                    siteLikelihoods.put(site, siteLikelihoods.get(site) + likelihood);
                }

                // siteDepths are mostly used for debugging and display
                siteDepths.put(site, siteDepths.get(site) + 1);
                siteReads.get(site).add(r);
            }
        }

        double[] vsLikelihoodArray = new double[siteLikelihoods.size()];
        int i = 0;
        for(VariantSite vs : siteLikelihoods.keySet()) {
            vsLikelihoodArray[i++] = siteLikelihoods.get(vs);
        }

        double sum = 0.0;
        for(i = 0; i < vsLikelihoodArray.length; i++) {
            sum += vsLikelihoodArray[i];
        }

        return sum;
    }

}
