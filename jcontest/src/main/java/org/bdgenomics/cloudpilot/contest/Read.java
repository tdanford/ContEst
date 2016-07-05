package org.bdgenomics.cloudpilot.contest;

import java.util.Random;
import java.util.UUID;

public class Read {

    public String id;
    public String sequence, qualities;
    public Alignment alignment;

    private static final double log10 = Math.log(10.0);

    public Read(String seq, double[] errors, Alignment align) {
        this(seq, encodeErrors(errors), align);
    }

    public Read(String seq, String qual, Alignment align) {
        this.id = UUID.randomUUID().toString();
        this.sequence = seq;
        this.qualities = qual;
        this.alignment = align;
    }

    public int hashCode() {
        return id.hashCode();
    }

    public static String encodeErrors(double[] e) {
        StringBuilder sb = new StringBuilder();
        for(int i = 0; i < e.length; i++) {
            sb.append(encodeError(e[i]));
        }

        return sb.toString();
    }

    public static double decodeError(char c) {
        int offset = c - '!';
        return Math.exp(log10 * offset / -10.0);
    }

    public static char encodeError(double e) {
        int offset = (int)Math.floor(-10.0 * Math.log(e) / log10); /// TODO fix me
        return (char)('!' + offset);
    }

    public static double[] sampleErrors(Random rand, int length, double meanError) {
        double[] es = new double[length];
        for(int i = 0; i < es.length; i++) {
            es[i] = sampleError(rand, meanError);
        }
        return es;
    }

    public static double sampleError(Random rand, double meanError) {
        double error = 1.0;
        do {
            error = meanError * -Math.log(rand.nextDouble());
        } while(error > 1.0);
        return error;
    }

    public Range region() { return new Range(alignment.chrom(), alignment.offset(), alignment.offset() + sequence.length()); }

    public double errorProb(int offset) {
        return decodeError(qualities.charAt(offset));
    }

    public boolean equals(Object o) {
        if(!(o instanceof Read)) { return false; }
        Read r = (Read)o;
        return r.id.equals(id);
    }

    public String toString() {
        return String.format("%s:(%s // %s)@%s", id, sequence, qualities, alignment);
    }
}
