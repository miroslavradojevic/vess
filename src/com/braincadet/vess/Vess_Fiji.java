package com.braincadet.vess;


import features.TubenessProcessor;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.io.OpenDialog;
import ij.plugin.ImageCalculator;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Experimental plugin that runs Fiji's tubularity (dependent on VIB-lib-2.1.0.jar library from Fiji distribution)
 *
 */
public class Vess_Fiji implements PlugIn {

    int N, M, P; // width, height, stack size
    String imdir, imnameshort;
    String sigs = "2,4"; // Gaussian filtering scales
    ArrayList<Float> sig = new ArrayList<Float>(); // scales as list of floating point numbers

    /**
     * Default run method for the ImageJ plugin.
     * @param s IJ Plugin argument
     */
    public void run(String s) {

        // read input image, store the most recent path in Prefs
        String in_folder = Prefs.get("com.braincadet.vess.dir", System.getProperty("user.home"));
        OpenDialog.setDefaultDirectory(in_folder);
        OpenDialog dc = new OpenDialog("Select image");
        in_folder = dc.getDirectory();
        Prefs.set("com.braincadet.vess.dir", in_folder);
        String image_path = dc.getPath();
        if (image_path == null) return;

        ImagePlus ip_load = new ImagePlus(image_path);

        if (ip_load == null) {
            IJ.log(image_path + " was null");
            return;
        }
        if (ip_load.getType() != ImagePlus.GRAY8) {
            IJ.log("Image needs to be GRAY8.");
            return;
        }

        N = ip_load.getWidth();
        M = ip_load.getHeight();
        P = ip_load.getStack().getSize();

        ip_load.setCalibration(null);

        imnameshort = ip_load.getShortTitle();
        imdir = ip_load.getOriginalFileInfo().directory;

        GenericDialog gd = new GenericDialog("Vesselness");
        gd.addStringField("sigmas",      Prefs.get("com.braincadet.vess.sigmas", sigs), 10);
        gd.showDialog();
        if (gd.wasCanceled()) return;
        sigs = gd.getNextString();
        Prefs.set("com.braincadet.vess.sigmas", sigs);

        IJ.log("source="+image_path);
        IJ.log("sigmas="+sigs);

        String[] readLn = sigs.trim().split(","); // read vesselness scales into list
        sig.clear();
        for (int i = 0; i < readLn.length; i++) {
            sig.add(   Float.valueOf(readLn[i].trim()).floatValue()   );
        }

        //------------------------------------------------------------------------
        long t1, t2;
        t1 = System.currentTimeMillis();

        ImageStack is_tness = new ImageStack(N, M); // prepare stack, fill with zeros
        for (int i = 0; i < P; i++) {
            float[] tt = new float[N * M];
            Arrays.fill(tt, 0f);
            is_tness.addSlice(new FloatProcessor(N, M, tt));
        }
        ImagePlus ip_tness = new ImagePlus(imnameshort+"_Vess_Fiji", is_tness); // float image plus object
//        ip_load = new ImagePlus(image_path); // reassign!!
        for (int i = 0; i < sig.size(); i++) {
            TubenessProcessor tp = new TubenessProcessor((double)sig.get(i), false);
            ImagePlus result = tp.generateImage(ip_load);
            ImageCalculator ic = new ImageCalculator();
            ic.run("Max 32-bit stack", ip_tness, result);
        }
        ip_tness.setCalibration(null);
        t2 = System.currentTimeMillis();
        IJ.log("t = " + IJ.d2s((t2 - t1) / 1000f,2) + " [sec]");

        ip_tness.show();

        // tubeness min-max normalize and store in an array for later and extract locations in a separate array
//        float[] tness = floatimage2array(ip_tness);
//        float tnessmin = Float.POSITIVE_INFINITY;
//        float tnessmax = Float.NEGATIVE_INFINITY;
//        for (int i = 0; i < tness.length; i++) {
//            if (tness[i] < tnessmin) tnessmin = tness[i];
//            if (tness[i] > tnessmax) tnessmax = tness[i];
//        }
//        byte[] J8_2 = new byte[N*M*P];
//        for (int i = 0; i < J8_2.length; i++) {
//            J8_2[i] = (tnessmax - tnessmin > Float.MIN_VALUE) ?   (byte) Math.round(((tness[i] - tnessmin) / (tnessmax - tnessmin))*255)    : (byte) 0;
//        }
//        new ImagePlus(imnameshort+"_Vess_Fiji",  array2imagestack(J8_2, N, M, P)).show();

    }
}
