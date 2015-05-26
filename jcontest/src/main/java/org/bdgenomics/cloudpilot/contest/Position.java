package org.bdgenomics.cloudpilot.contest;

public class Position implements Comparable<Position> {

    private String chrom;
    private int offset;

    public Position(String c, int o) {
        this.chrom = c;
        this.offset = o;
    }

    public String chrom() { return chrom; }
    public int offset() { return offset; }

    public int hashCode() {
        return offset + 37 * (chrom.hashCode() + 17);
    }

    public boolean equals(Object o) {
        if(!(o instanceof Position)) { return false; }
        Position p = (Position)o;
        return chrom.equals(p.chrom) && offset == p.offset;
    }

    public String toString() { return String.format("%s:%d", chrom, offset); }

    public int compareTo(Position p) {
        if(!chrom.equals(p.chrom)) { return chrom.compareTo(p.chrom); }
        if(offset < p.offset) { return -1; }
        if(offset > p.offset) { return 1; }
        return 0;
    }

}
