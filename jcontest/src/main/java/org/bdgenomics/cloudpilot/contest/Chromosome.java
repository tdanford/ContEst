package org.bdgenomics.cloudpilot.contest;

import java.util.*;

public class Chromosome {

    public String name;
    public String sequence;

    public Chromosome(String n,String seq) {
        this.name = n;
        this.sequence = seq;
    }

    public int length() { return sequence.length(); }

    public static Chromosome randomChromosome(String name, int length) {
        return new Chromosome(name, SeqUtils.randomString(length));
    }

    public static ArrayList<VariantSite> randomVariantSites(Chromosome c, int spacing) {
        ArrayList<VariantSite> vs = new ArrayList<VariantSite>();

        for(int offset = spacing; offset < c.length() ; offset += spacing) {
            String ref = c.sequence.substring(offset, offset+1);
            String alt = SeqUtils.randomAlternate(ref, 0);

            vs.add(new VariantSite(c.name, offset, ref, alt, 0.25));
        }

        return vs;
    }
}
