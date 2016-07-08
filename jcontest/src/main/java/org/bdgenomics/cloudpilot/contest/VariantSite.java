package org.bdgenomics.cloudpilot.contest;

import org.apache.log4j.Logger;

import java.util.*;
import java.io.*;

public class VariantSite extends Position {

    private static Logger logger = Logger.getLogger(VariantSite.class);

    public String referenceAllele, alternateAllele;

    /*
    The probability of seeing the allele A_i [ the reference allele ] at site i in the
    contaminating population -- labeled f_i in the paper.
     */
    public double contaminatingPopulationFrequency;

    public VariantSite(String c, int loc, String ref, String alt, double f) {
        super(c, loc);
        this.referenceAllele = ref;
        this.alternateAllele = alt;
        this.contaminatingPopulationFrequency = f;
    }

    private static final List<String> chars = Arrays.asList("A", "T", "G", "C");

    public List<Read> sampleOverlappingReads(Genome g, double c, int readLength, double meanError, int readCount) {
        ArrayList<Read> reads = new ArrayList<>();
        for(int i = 0; i < readCount; i++) {
            reads.add(sampleOverlappingRead(g, c, readLength, SeqUtils.rand, meanError));
        }
        return reads;
    }

    public Read sampleOverlappingRead(Genome g, double c, int readLength, Random rand, double meanError) {
        int variantOffset = rand.nextInt(readLength);
        int readPosition = offset() - variantOffset;
        double[] errorProbs = Read.sampleErrors(rand, readLength, meanError);

        /*
        StringBuilder sb = new StringBuilder();
        for(double e : errorProbs) {
            sb.append(String.format(" %.04f", e));
        }
        logger.info(String.format("Errors: %s", sb.toString()));
        */

        String allele = referenceAllele;

        boolean isContaminant = rand.nextDouble() <= c;
        if(isContaminant) {
            boolean isAlternate = rand.nextDouble() <= contaminatingPopulationFrequency;
            if(isAlternate) {
                allele = alternateAllele;
            }
        }

        boolean isError = rand.nextDouble() <= errorProbs[variantOffset];
        if(isError) {
            int alleleIndex = chars.indexOf(allele);
            allele = chars.get((alleleIndex + 1 + rand.nextInt(chars.size()-1)) % chars.size());
        }

        StringBuilder seq = new StringBuilder(g.chroms.get(chrom()).sequence.substring(readPosition, readPosition+readLength));
        seq.setCharAt(variantOffset, allele.charAt(0));

        boolean strand = true;
        Alignment align = new Alignment(chrom(), readPosition, strand);
        return  new Read(seq.toString(), errorProbs, align);
    }


    public double readLogLikelihood(double c, Read r) { return Math.log(readLikelihood(c, r)); }

    public double readLikelihood( double c, Read r ) {
        int offset = offset() - r.alignment.offset();
        String readBase = r.sequence.substring(offset, offset+1);
        double error = r.errorProb(offset);
        double notError = 1.0 - error;
        double error3 = error / 3.0;
        double f = contaminatingPopulationFrequency;

        if(readBase.equals(referenceAllele)) {
            double not_c = (1.0-c) * notError;
            double _c = c * (f * notError + (1.0 - f) * error3);
            return not_c + _c;

        } else if (readBase.equals(alternateAllele)) {
            double not_c = (1.0 - c) * error3;
            double _c = c * (f * error3 + (1.0 - f) * notError);
            return not_c + _c;

        } else {
            return error3;
        }
    }

    public int hashCode() {
        int code = 17;
        code += super.hashCode(); code *= 37;
        code += referenceAllele.hashCode(); code *= 37;
        return code;
    }

    public boolean equals(Object o) {
        if(!(o instanceof VariantSite)) { return false; }
        VariantSite v = (VariantSite)o;
        return super.equals(v) &&
                v.referenceAllele.equals(referenceAllele);
    }

    public static Iterator<VariantSite> vcfVariants(File f) throws IOException {
        return new VariantSiteIterator(new FileInputStream(f));
    }

    /**
     * Parses VariantSites out of a VCF file.
     */
    public static class VariantSiteIterator implements Iterator<VariantSite> {

        private BufferedReader br;
        private LinkedList<VariantSite> nextSites;

        public VariantSiteIterator(InputStream is) throws IOException {
            br = new BufferedReader(new InputStreamReader(is, "UTF-8"));
            nextSites = new LinkedList<VariantSite>();
            findNextSite();
        }

        private void findNextSite() {
            if(!nextSites.isEmpty()) { return; }
            String line = null;
            try {
                while (br != null && (line = br.readLine()) != null && line.startsWith("#")) {
                    line = null;
                }

                if(line != null) {
                    String[] array = line.split("\t");
                    String chrom = array[0];
                    int pos = Integer.parseInt(array[1]);
                    String id = array[2];
                    String ref = array[3];
                    String[] alts = array[4].split(",");

                    for(int i = 0; i < alts.length; i++) {
                        nextSites.addLast(new VariantSite(chrom, pos, ref, alts[i], 0.0));
                    }
                } else {
                    if(br != null) {
                        br.close();
                    }
                }

            } catch(IOException e) {
                e.printStackTrace(System.err);
            }
        }

        @Override
        public boolean hasNext() { return !nextSites.isEmpty(); }

        @Override
        public VariantSite next() {
            VariantSite vs = nextSites.removeFirst();
            findNextSite();
            return vs;
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException();
        }
    }

}
