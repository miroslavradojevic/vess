package com.braincadet.vess;

import features.TubenessProcessor;
import ij.*;
import ij.gui.GenericDialog;
import ij.io.OpenDialog;
import ij.plugin.ImageCalculator;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

public class Vess implements PlugIn {

    int N, M, P, SZ;
    String imdir, imnameshort;

    String sigs = "2,4,6";
    ArrayList<Float> sig = new ArrayList<Float>();
    float alpha = 0.5f;
    float beta = 0.5f;
    float BetaOne = 0.5f;
    float BetaTwo = 15;
    float C = 500;
    float zdist = 2.0f;
    boolean blackwhite = false;

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
        SZ = N * M * P;

        ip_load.setCalibration(null);

        imnameshort = ip_load.getShortTitle();
        imdir = ip_load.getOriginalFileInfo().directory;

        GenericDialog gd = new GenericDialog("Vesselness");
        gd.addStringField("sigmas", Prefs.get("com.braincadet.vess.sigmas", sigs), 10);
        gd.addNumericField("alpha", Prefs.get("com.braincadet.vess.alpha", alpha), 2, 5, " (3D)");
        gd.addNumericField("beta", Prefs.get( "com.braincadet.vess.beta", beta), 2, 5, " (3D)");
        gd.addNumericField("betaone", Prefs.get("com.braincadet.vess.betaone", BetaOne), 2, 5, " (2D)");


        gd.showDialog();
        if (gd.wasCanceled()) return;

        sigs = gd.getNextString();
        Prefs.set("com.braincadet.vess.sigmas", sigs);

        alpha = (float) gd.getNextNumber();
        Prefs.set("com.braincadet.vess.alpha", alpha);

        beta = (float) gd.getNextNumber();
        Prefs.set("com.braincadet.vess.beta", beta);

        BetaOne = (float) gd.getNextNumber();
        Prefs.set("com.braincadet.vess.betaone", beta);

        IJ.log("ip_load="+image_path);
        IJ.log("sigmas="+sigs);
        IJ.log("alpha="+alpha);
        IJ.log("beta="+beta);
        IJ.log("BetaOne="+BetaOne);


        IJ.log(" -- prefiltering...");

        long t1prep = System.currentTimeMillis();

        ImageStack is_tness = new ImageStack(N, M);
        for (int i = 0; i < P; i++) {
            float[] tt = new float[N * M];
            Arrays.fill(tt, 0f);
            is_tness.addSlice(new FloatProcessor(N, M, tt));
        }

        ImagePlus ip_tness = new ImagePlus("tness", is_tness);
        String[] readLn = sigs.trim().split(",");

        for (int i = 0; i < readLn.length; i++) {
            float sig = Float.valueOf(readLn[i].trim()).floatValue();
            TubenessProcessor tp = new TubenessProcessor(sig, false);
            ImagePlus result = tp.generateImage(ip_load);
            ImageCalculator ic = new ImageCalculator();
            // average, multipy
//          IJ.run(result, "Multiply...", "value=" + IJ.d2s(1f/readLn.length,3) + " stack");
//          ic.run("Add 32-bit stack", ip_tness, result); // result of the addition is placed in ip_tness
            // max
            ic.run("Max 32-bit stack", ip_tness, result);
        }

        ip_tness.setCalibration(null);

        if (true) {
            ImagePlus temp = ip_tness.duplicate();
            IJ.run(temp, "8-bit", ""); // convert to 8 bit before saving
            IJ.saveAs(temp, "Zip", imdir + File.separator + "tness," + sigs + ".zip");
        }

        // tubeness min-max normalize and store in an array for later and extract locations in a separate array
        float tnessmin = Float.POSITIVE_INFINITY;
        float tnessmax = Float.NEGATIVE_INFINITY;
        float[] tness = new float[SZ];

        for (int z = 1; z <= P; z++) { // layer count, zcoord is layer-1

            float[] slc_float = (float[]) ip_tness.getStack().getPixels(z);

            for (int x = 0; x < N; x++) {
                for (int y = 0; y < M; y++) {
                    int ii = (z - 1) * (N * M) + y * N + x;
                    tness[ii] = slc_float[y * N + x];
                    if (tness[ii] < tnessmin) tnessmin = tness[ii];
                    if (tness[ii] > tnessmax) tnessmax = tness[ii];
                }
            }
        }

        for (int i = 0; i < SZ; i++) {
            tness[i] = (tnessmax - tnessmin > Float.MIN_VALUE) ? ((tness[i] - tnessmin) / (tnessmax - tnessmin)) : 0;
        }

        long t2prep = System.currentTimeMillis();
        IJ.log("t_prep = " + IJ.d2s((t2prep - t1prep) / 1000f,2) + " [sec]");

        if (true) return;

        Frangi f = new Frangi(sig, zdist, alpha, beta, C, BetaOne, BetaTwo);

        byte[] ip_array = f.ip2array(ip_load);

        // delete ip_load (free up some memory)
        IJ.log("removing slices... " + ip_load.getStack().getSize());
        while (ip_load.getStack().getSize()>1) {
            IJ.log("before " + ip_load.getStack().getSize());
            ip_load.getStack().deleteLastSlice();
            IJ.log("after " + ip_load.getStack().getSize());
        }

        ip_load = null;

        IJ.log("done");


        float[] fgauss = new float[ip_array.length];
        f.imgaussian(ip_array, N, M, 4f, fgauss);

        new ImagePlus("orig", f.array2imagestack(ip_array, N,M,P)).show();
        new ImagePlus("gauss", f.array2imagestack(fgauss, N, M, P)).show();

//        byte[] fgauss_byte8 = new byte[fgauss.length];
//        f.float2byte(fgauss, fgauss_byte8);
//        new ImagePlus("", f.array2imagestack(fgauss_byte8, N, M, P)).show();

    }

}
