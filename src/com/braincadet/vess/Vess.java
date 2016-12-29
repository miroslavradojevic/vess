package com.braincadet.vess;

import ij.IJ;
import ij.ImagePlus;
import ij.Macro;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;

import java.util.ArrayList;

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

        gd.showDialog();
        if (gd.wasCanceled()) return;

        sigs = gd.getNextString();
        Prefs.set("com.braincadet.vess.sigmas", sigs);

            alpha = (float) gd.getNextNumber();
            Prefs.set("com.braincadet.vess.alpha", alpha);

            beta = (float) gd.getNextNumber();
            Prefs.set("com.braincadet.vess.beta", beta);

        IJ.log("ip_load="+image_path);
        IJ.log("sigmas="+sigs);
        IJ.log("alpha="+alpha);
        IJ.log("beta="+beta);


        Frangi f = new Frangi(sig, zdist, alpha, beta, C, BetaOne, BetaTwo);

        byte[] ip_array = f.ip2array(ip_load);

        // delete ip_load (free up some memory)
//        while (ip_load.getStack().getSize()>0) {
//            ip_load.getStack().deleteLastSlice();
//        }

        new ImagePlus("", f.array2imagestack(ip_array, N, M, P)).show();
        ip_load.show();

        IJ.log((Integer.MAX_VALUE)+"");
        IJ.log((N*M*P)+"");

        IJ.log(""+ (float)(Integer.MAX_VALUE/(N*M*P)));


    }

}
