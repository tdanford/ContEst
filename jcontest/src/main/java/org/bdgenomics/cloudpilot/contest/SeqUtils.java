package org.bdgenomics.cloudpilot.contest;

import java.util.*;



public class SeqUtils {

    public static final char[] letters = { 'A', 'T', 'G', 'C' };
    public static final List<Character> letterList = Arrays.asList('A', 'T', 'G', 'C');

    public static Random rand = new Random();

    public static char randomLetter() { return letters[rand.nextInt(letters.length)]; }

    public static char randomAlternate(char c) {
        int idx = letterList.indexOf(c);
        return letters[(idx + rand.nextInt(letters.length-1) + 1) % letters.length];
    }

    public static String randomString(int len) {
        StringBuilder sb = new StringBuilder();
        for(int i = 0 ; i < len; i++) {
            sb.append(randomLetter());
        }
        return sb.toString();
    }

    public static String randomAlternate(String base, int alternateOffset) {
        StringBuilder sb = new StringBuilder(base);
        sb.setCharAt(alternateOffset, randomAlternate(sb.charAt(alternateOffset)));
        return sb.toString();
    }

}
