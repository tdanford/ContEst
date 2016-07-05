package org.bdgenomics.cloudpilot.contest;

import java.util.*;

public class Genome {

    public Map<String,Chromosome> chroms;

    public Genome() {}

    public Genome(Chromosome... cs) {
        this.chroms = new TreeMap<String,Chromosome>();
        for(Chromosome c : cs) {
            chroms.put(c.name, c);
        }
    }

    public ArrayList<VariantSite> randomVariantSites(int spacing) {
        ArrayList<VariantSite> vs = new ArrayList<VariantSite>();
        for(Chromosome c : chroms.values()) {
            vs.addAll(Chromosome.randomVariantSites(c, spacing));
        }
        return vs;
    }

    public static Genome randomGenome(int chromCount, int chromLength) {
        Chromosome[] cs = new Chromosome[chromCount];
        for(int i =0 ; i < cs.length; i++) {
            cs[i] = Chromosome.randomChromosome("chrom" + i, chromLength);
        }
        return new Genome(cs);
    }
}
