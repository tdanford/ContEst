package org.bdgenomics.cloudpilot.contest;

import java.util.*;

public interface Accumulator<T> {

    void process(T value);
    boolean shouldContinue();


    class Counter<T> implements Accumulator<T> {
        private int count = 0;
        public void process(T value) { count += 1; }
        public boolean shouldContinue() { return true; }
        public static <T> Counter<T> create() { return new Counter<T>(); }

        public int count() { return count; }
    }

    class Array<T> implements Accumulator<T> {

        public static <T> Array<T> create(T[] tarray) { return new Array<T>(tarray); }

        private T[] tarray = null;
        private ArrayList<T> list = new ArrayList<T>();

        public Array(T[] tarray) { this.tarray = tarray; }

        public void process(T value) { list.add(value); }
        public boolean shouldContinue() { return true; }

        public T[] array() { return (T[])list.toArray(tarray); }
    }

    class List<T> implements Accumulator<T> {

        public static <T> List<T> create() { return new List<T>(); }

        private ArrayList<T> list = new ArrayList<T>();

        public void process(T value) { list.add(value); }
        public boolean shouldContinue() { return true; }

        public ArrayList<T> list() { return list; }
    }

    class First<T> implements Accumulator<T> {

        public static <T> First<T> create() { return new First<T>(); }

        private T firstValue = null;

        public void process(T value) {
            if(firstValue != null) {
                firstValue = value;
            }
        }

        public boolean shouldContinue() { return firstValue != null; }

        public T value() { return firstValue; }
    }
}
