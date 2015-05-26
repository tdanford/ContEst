package org.bdgenomics.cloudpilot.contest;

import org.apache.log4j.Logger;

import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.Collection;

public class PositionTree<P extends Position> {

    private static Logger logger = Logger.getLogger(PositionTree.class);

    public Position split;
    public PositionTree<P> leftTree, rightTree;

    public P[] positions;
    public int start, end;

    public static <P extends Position> PositionTree<P> createTree(Collection<P> pos) {
        if(pos.isEmpty()) {
            return new PositionTree<P>();
        }

        P first = pos.iterator().next();
        P[] array = (P[])Array.newInstance(first.getClass(), pos.size());
        array = (P[])pos.toArray(array);
        Arrays.sort(array);

        logger.debug(String.format("inserting %s into PositionTree", String.valueOf(Arrays.asList(array))));

        return new PositionTree<>(array, 0, array.length);
    }

    public PositionTree() {
        this.start = this.end = -1;
        this.positions = (P[])null;
    }

    public PositionTree(P[] positions, int start, int end) {
        logger.debug(String.format("PositionTree(%d,%d)", start, end));
        this.positions = positions;
        this.start = start;
        this.end = end;

        if(!positions[start].equals(positions[end-1])) {
            int splitIndex = findSplit(positions, start, end);

            if(splitIndex > start && splitIndex < end) {
                split = positions[splitIndex];
                logger.debug(String.format("Found split %d: %s", splitIndex, split));

                leftTree = new PositionTree<P>(positions, start, splitIndex);
                rightTree = new PositionTree<P>(positions, splitIndex, end);
            }
        } else {
            logger.debug(String.format("Values %s into PositionTree leaf", String.valueOf(Arrays.asList(positions).subList(start, end))));
        }
    }

    public <A extends Accumulator<P>> A findOverlapping(Range query, A acc) {

        if(split != null) {
            if(query.contains(split)) {
                acc = leftTree.findOverlapping(query.splitLeft(split), acc);
                return rightTree.findOverlapping(query.splitRight(split), acc);

            } else {
                if(query.startPosition().compareTo(split) < 0) {
                    return leftTree.findOverlapping(query, acc);
                } else {
                    return rightTree.findOverlapping(query, acc);
                }
            }

        } else {
            for(int i = start; i < end; i++) {
                if(query.contains(positions[i])) {
                    acc.process(positions[i]);
                }
            }
            return acc;
        }
    }

    private static <P extends Position> int findSplit(P[] positions, int start, int end) {
        int k = (start+end)/2;
        while(k < end && positions[k].equals(positions[k+1])) {
            k += 1;
        }
        return k+1;
    }
}
