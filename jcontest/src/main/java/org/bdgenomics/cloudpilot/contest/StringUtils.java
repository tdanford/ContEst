package org.bdgenomics.cloudpilot.contest;

public class StringUtils {


    public static String leftPad(String base, int length, char pad) {
        StringBuilder sb = new StringBuilder();
        for(int i = 0; i < length - base.length(); i++) {
            sb.append(pad);
        }
        sb.append(base);
        return sb.toString();
    }

    public static String rightPad(String base, int length, char pad) {
        StringBuilder sb = new StringBuilder();
        sb.append(base);
        for(int i = 0; i < length - base.length(); i++) {
            sb.append(pad);
        }
        return sb.toString();
    }

    public static String upperCaseOffset(String base, int offset) {
        StringBuilder sb = new StringBuilder(base.toLowerCase());
        sb.setCharAt(offset, Character.toUpperCase(sb.charAt(offset)));
        sb.insert(offset+1, ' ');
        sb.insert(offset, ' ');
        return sb.toString();
    }

}
