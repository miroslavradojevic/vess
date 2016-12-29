package com.braincadet.vess;

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

public class Test implements PlugInFilter {
    public int setup(String s, ImagePlus imagePlus) {
        long nr_bytes = imagePlus.getWidth()*imagePlus.getHeight()*imagePlus.getStackSize(); // byte images used
        IJ.log(imagePlus.getShortTitle() + "," + bytesToMegabytes(nr_bytes) + "," + ((float)nr_bytes/Integer.MAX_VALUE));
        return DOES_8G;
    }
    public void run(ImageProcessor imageProcessor) {}

    private static final long MEGABYTE = 1024L * 1024L;

    public static long bytesToMegabytes(long bytes) {
        return bytes / MEGABYTE;
    }
}
