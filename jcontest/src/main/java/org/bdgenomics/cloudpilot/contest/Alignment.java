package org.bdgenomics.cloudpilot.contest;

public class Alignment extends Position {

    public boolean strand;

    public Alignment(String c, int off, boolean str) {
        super(c, off);
        this.strand = str;
    }

    public int hashCode() {
        int code = 17;
        code += super.hashCode(); code *= 37;
        code += (strand ? 1 : 0) ; code *= 37;
        return code;
    }

    public boolean equals(Object o) {
        if(!(o instanceof Alignment)) { return false; }
        Alignment a = (Alignment)o;
        return super.equals(o) && strand == a.strand;
    }

    public String toString() { return String.format("%s:%s", super.toString(), strand? "+" : "-"); }
}
