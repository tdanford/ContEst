package org.bdgenomics.cloudpilot.contest;

import org.junit.Test;
import static org.junit.Assert.*;

public class ReadQualityScoreTest {

    @Test
    public void testQualityCharRoundTrip() {
        char phredChar = (char)('!' + 25);
        double error = Read.decodeError(phredChar);
        char reencoded = Read.encodeError(error);

        assertEquals(phredChar, reencoded);
    }

    @Test
    public void testQualityScoreRoundTrip() {
        double error = 0.015;
        char encoded = Read.encodeError(error);
        double decoded = Read.decodeError(encoded);

        assertTrue(Math.abs(error - decoded) <= 0.001);
    }



}
