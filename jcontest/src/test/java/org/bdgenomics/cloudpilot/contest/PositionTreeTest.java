package org.bdgenomics.cloudpilot.contest;

import org.apache.log4j.Logger;
import org.junit.Test;

import java.util.*;

import static org.junit.Assert.*;


public class PositionTreeTest {

    private static Logger logger = Logger.getLogger(PositionTreeTest.class);

    @Test
    public void testPositionTreeCreation() {
        ArrayList<Position> pos = new ArrayList<>();
        for(int i = 0; i < 10; i++) {
            pos.add(new Position("chrom", i));
        }
        PositionTree<Position> ptree = PositionTree.createTree(pos);
        assertNotNull(ptree);
    }

    @Test
    public void testEmptyPositionTreeCreation() {
        assertNotNull(PositionTree.createTree(new ArrayList<Position>()));
    }

    @Test
    public void testOnePositionTreeCreation() {
        assertNotNull(PositionTree.createTree(new ArrayList<Position>(Arrays.asList(new Position("foo", 10)))));
    }

    @Test
    public void testSimpleOverlap() {
        Position[] ps = new Position[] {
                new Position("foo", 10),
                new Position("foo", 50),
                new Position("foo", 51)
        };

        PositionTree<Position> ptree = PositionTree.createTree(Arrays.asList(ps));
        assertArrayEquals(
                new Position[]{new Position("foo", 10)},
                ptree.findOverlapping(new Range("foo", 5, 15), Accumulator.Array.create(new Position[0])).array());
    }

    @Test
    public void testMultipleOverlapping() {
        ArrayList<Position> ps = new ArrayList<Position>();
        for(int i = 0; i < 10; i++) {
            ps.add(new Position("foo", i));
            ps.add(new Position("foo", i));
        }

        logger.info(String.valueOf(ps));

        Range query = new Range("foo", 0, 5);

        PositionTree<Position> ptree = PositionTree.createTree(ps);
        assertEquals(
                String.valueOf(ptree.findOverlapping(query, new Accumulator.List()).list()),
                10,
                ptree.findOverlapping(query, new Accumulator.Counter<Position>()).count());
    }

}
