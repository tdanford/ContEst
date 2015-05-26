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

    public static int argMax(double[] values) {
        int maxId = -1;
        for(int i = 0; i < values.length; i++) {
            if(i == 0 || values[i] >= values[maxId]) {
                maxId = i;
            }
        }
        return maxId;
    }

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

    public double logLikelihoodC(double c) {
        Map<VariantSite, Double> siteLikelihoods = new HashMap<>();
        Map<VariantSite, Integer> siteDepths = new HashMap<>();
        Map<VariantSite, ArrayList<Read>> siteReads = new HashMap<>();

        for(VariantSite vs : variants) {
            siteReads.put(vs, new ArrayList<Read>());
            siteDepths.put(vs, 0);
        }

        for(Read r : reads) {
            for(VariantSite site : variantTree.findOverlapping(r.region(), new Accumulator.List<VariantSite>()).list()) {
                double likelihood = site.readLogLikelihood(c, r);
                if(!siteLikelihoods.containsKey(site)) {
                    siteLikelihoods.put(site, likelihood);
                } else {
                    siteLikelihoods.put(site, siteLikelihoods.get(site) + likelihood);
                }

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

        /*
        for(VariantSite vs : variants) {
            logger.info(String.format("Site %s: depth %d", vs, siteDepths.get(vs)));
            printSiteAndReads(vs, siteReads.get(vs));
            System.out.println();
        }
        */

        return sum;
    }

}
