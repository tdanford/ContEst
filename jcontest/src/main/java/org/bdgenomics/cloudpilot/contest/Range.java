package org.bdgenomics.cloudpilot.contest;

public class Range implements Comparable<Range> {

    private String chrom;
    private int start, end;

    public Range(String c, int s, int e) {

        this.chrom = c;
        this.start = s;
        this.end = e;

        if(s < 0) { throw new IllegalArgumentException("start was negative"); }
        if(s > e) { throw new IllegalArgumentException("start was greater than end"); }
        if(c == null) { throw new IllegalArgumentException("'chrom' was null"); }
    }

    public int start() { return start; }
    public int end() { return end; }
    public String chrom() { return chrom; }
    public int length() { return end - start; }

    public Position startPosition() { return new Position(chrom, start); }
    public Position endPosition() { return new Position(chrom, end); }

    public boolean contains(Position p) {
        return chrom.equals(p.chrom()) && start <= p.offset() && end > p.offset();
    }

    public boolean overlaps(Range r) {
        return chrom.equals(r.chrom) && start < r.end() && end > r.start();
    }

    public Range splitLeft(Position p) {
        if(!contains(p)) {
            return this;
        } else {
            return new Range(chrom, start, p.offset());
        }
    }

    public Range splitRight(Position p) {
        if(!contains(p)) {
            return this;
        } else {
            return new Range(chrom, p.offset(), end);
        }
    }

    public int hashCode() {
        int code = 17;
        code += chrom.hashCode(); code *= 37;
        code += start; code *= 37;
        code += end; code *= 37;
        return code;
    }

    public boolean equals(Object o) {
        if(!(o instanceof Range)) { return false; }
        Range r = (Range)o;
        return r.chrom.equals(chrom) && start == r.start && end == r.end;
    }

    public int compareTo(Range r) {
        if(!chrom.equals(r.chrom)) { return chrom.compareTo(r.chrom); }
        if(start < r.start) { return -1; }
        if(start > r.start) { return 1; }
        if(length() > r.length()) { return -1; }
        if(length() < r.length()) { return 1; }
        return 0;
    }
}
