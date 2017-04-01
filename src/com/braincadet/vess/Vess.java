package com.braincadet.vess;

import ij.*;
import ij.gui.GenericDialog;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageConverter;

import java.io.*;
import java.util.ArrayList;

/**
 * <h1>ImageJ plugin that computes Frangi's vesselness.</h1>
 *
 * Uses Frangi class
 * <p>
 * The Vess runs Frangi's vesselness for selected byte8 2d/3d image.
 * Gives the vesselness image and the eigenvectors.
 * Threaded implementation available.
 * </p>
 * <p>
 *     <i>Usage:</i> <b>Plugins>BrainCadet, "Vess"</b>
 * </p>
 *
 * @author  Miroslav Radojevic
 * @version 1.0
 * @since   2017-03-09
 */
public class Vess implements PlugIn {

    int N, M, P; // width, height, stack size

    String imdir, imnameshort, outdir;

    // parameters wit the default initial values
    String sigs = "2,4"; // Gaussian filtering scales
    ArrayList<Float> sig = new ArrayList<Float>(); // scales as list of floating point numbers
    float zdist = 2.0f;
    float alpha = 0.5f;
    float beta = 0.5f;
    float C = 500f;
    float BetaOne = 0.5f;
    float BetaTwo = 15f;
    boolean darkforeground = false;
    boolean directions = true;
    boolean saveoutput = true;
    boolean threaded3d = true;

    // simple struct containing image information
    private class Image8 {
        public int Width;
        public int Height;
        public int Length;
        String short_title;
        String image_dir;
        byte[] data; // byte8 data array
    }

    private static final long MEGABYTE = 1024L * 1024L;

    private static long BytesToMegabytes(long bytes) {
        return bytes / MEGABYTE;
    }

//    IJ.log("TOTAL MEMORY = "+ BytesToMegabytes( Runtime.getRuntime().totalMemory()) + " Mb" ); // memory logging

    private Image8 load_image(String image_path){

        ImagePlus input_image = new ImagePlus(image_path);

        if (input_image == null) {
            IJ.log(image_path+" == null");
            return null;
        }

        if (input_image.getType() != ImagePlus.GRAY8) {
            IJ.log(image_path + " is not GRAY8.");
            return null;
        }

        Image8 img = new Image8();
        img.Width   = input_image.getWidth();
        img.Height  = input_image.getHeight();
        img.Length  = input_image.getStack().getSize();
        img.short_title = input_image.getShortTitle();
        img.image_dir = input_image.getOriginalFileInfo().directory;

        img.data = new byte[img.Width*img.Height*img.Length]; // read image into byte[]
        for (int z = 1; z <= img.Length; z++) { // layer count, zcoord is layer-1
            byte[] slc = (byte[]) input_image.getStack().getPixels(z);
            for (int x = 0; x < img.Width; x++) {
                for (int y = 0; y < img.Height; y++) {
                    img.data[(z - 1) * (img.Width * img.Height) + y * img.Width + x] = slc[y * img.Width + x];// & 0xff;
                }
            }
        }

        return img;

    }

    /**
     * PlugIn main method that runs the ImageJ plugin.
     * @param s String PlugIn run argument.
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

        Image8 ip_load8 = load_image(image_path); // handy wrapper not to keep ImagePlus object in memory during plugin runtime
        if (ip_load8 == null) {IJ.log("could not load " + image_path); return;}

        N = ip_load8.Width;
        M = ip_load8.Height;
        P = ip_load8.Length;

        imnameshort = ip_load8.short_title;
        imdir       = ip_load8.image_dir;

        GenericDialog gd = new GenericDialog("Vesselness");
        gd.addStringField("sigmas",         Prefs.get("com.braincadet.vess.sigmas", sigs), 10);
        gd.addCheckbox("darkforeground",    Prefs.get("com.braincadet.vess.darkforeground", darkforeground));
        gd.addCheckbox("directions",        Prefs.get("com.braincadet.vess.directions", directions));
        gd.addCheckbox("saveoutput",        Prefs.get("com.braincadet.vess.saveoutput", saveoutput));

        if (P>1) { // 3d
            gd.addNumericField("alpha",     Prefs.get("com.braincadet.vess.alpha", alpha), 2, 10, "");
            gd.addNumericField("beta",      Prefs.get("com.braincadet.vess.beta", beta),   2, 10, "");
            gd.addNumericField("c",         Prefs.get("com.braincadet.vess.c", C),         2, 10, "");
            gd.addNumericField("zdist",     Prefs.get("com.braincadet.vess.zdist", zdist), 2, 10, "");
            gd.addCheckbox("threaded3d",    Prefs.get("com.braincadet.vess.threaded3d", threaded3d));
        }
        else { // 2d
            gd.addNumericField("betaone",   Prefs.get("com.braincadet.vess.betaone", BetaOne), 2, 10, "");
            gd.addNumericField("betatwo",   Prefs.get("com.braincadet.vess.betatwo", BetaTwo), 2, 10, "");
        }

        gd.showDialog();
        if (gd.wasCanceled()) return;

        sigs = gd.getNextString();
        Prefs.set("com.braincadet.vess.sigmas", sigs);

        darkforeground = gd.getNextBoolean();
        Prefs.set("com.braincadet.vess.darkforeground", darkforeground);

        directions = gd.getNextBoolean();
        Prefs.set("com.braincadet.vess.directions", directions);

        saveoutput = gd.getNextBoolean();
        Prefs.set("com.braincadet.vess.saveoutput", saveoutput);

        if (P>1) {
            alpha = (float) gd.getNextNumber();
            Prefs.set("com.braincadet.vess.alpha", alpha);

            beta = (float) gd.getNextNumber();
            Prefs.set("com.braincadet.vess.beta", beta);

            C = (float) gd.getNextNumber();
            Prefs.set("com.braincadet.vess.c", C);

            zdist = (float) gd.getNextNumber();
            Prefs.set("com.braincadet.vess.zdist", zdist);

            threaded3d = gd.getNextBoolean();
            Prefs.set("com.braincadet.vess.threaded3d", threaded3d);

        }
        else {
            BetaOne = (float) gd.getNextNumber();
            Prefs.set("com.braincadet.vess.betaone", BetaOne);

            BetaTwo = (float) gd.getNextNumber();
            Prefs.set("com.braincadet.vess.betatwo", BetaTwo);

        }

        IJ.log("---");
        IJ.log("source="+image_path);
        IJ.log("sigmas="+sigs);

        IJ.log("darkforeground="+ darkforeground);
        IJ.log("directions="+directions);
        IJ.log("saveoutput="+saveoutput);

        IJ.log("zdist="+zdist);
        IJ.log("alpha="+alpha);
        IJ.log("beta="+beta);
        IJ.log("C="+C);
        IJ.log("threaded3d="+threaded3d);

        IJ.log("BetaOne="+BetaOne);
        IJ.log("BetaTwo="+BetaTwo);
        IJ.log("--");

        String delindir;

        if (P>1) {
            delindir = "VESS.sig.dar.dir.sav.alp.bet.c.zdi.thr_" +
                    sigs                        + "_" +
                    ((darkforeground)?"1":"0")  + "_" +
                    ((directions)?"1":"0")      + "_" +
                    ((saveoutput)?"1":"0")      + "_" +
                    IJ.d2s(alpha,2)             + "_" +
                    IJ.d2s(beta, 2)             + "_" +
                    IJ.d2s(C, 2)                + "_" +
                    IJ.d2s(zdist, 2)            + "_" +
                    ((threaded3d)?"1":"0")      + File.separator;
        }
        else {
            delindir = "VESS.sig.dar.dir.sav.betaone.betatwo_" +
                    sigs                        + "_" +
                    ((darkforeground)?"1":"0")  + "_" +
                    ((directions)?"1":"0")      + "_" +
                    ((saveoutput)?"1":"0")      + "_" +
                    IJ.d2s(BetaOne,2)           + "_" +
                    IJ.d2s(BetaTwo, 2)          + File.separator;
        }

        outdir      = imdir + delindir; // imnameshort + "_vess" + File.separator;
//        createAndCleanDir(outdir);
        createDir(outdir);

        String[] readLn = sigs.trim().split(","); // read vesselness scales into list
        sig.clear();
        for (int i = 0; i < readLn.length; i++) sig.add(Float.valueOf(readLn[i].trim()).floatValue());

        byte[] I = ip_load8.data; // byteimage2array(ip_load); // store it in byte array for processing

        float[] J; // output vesselness in 1d array

        byte[]  Vx   = (directions)? new byte[I.length] : null; // set Vx,Vy,Vz=null if not necessary to calculate Vx, Vy, Vz
        byte[]  Vy   = (directions)? new byte[I.length] : null;
        byte[]  Vz   = (directions)? new byte[I.length] : null;

        Frangi f = new Frangi(sig, zdist, alpha, beta, C, BetaOne, BetaTwo, darkforeground);

        long t1,t2; // time variables

        if (P>1){
            // 3d
            t1 = System.currentTimeMillis();
            J = f.run3d_float(I, N, M, P, Vx, Vy, Vz, threaded3d);
            t2 = System.currentTimeMillis();
            IJ.log("t = " + IJ.d2s((t2 - t1) / 1000f,2) + " [sec]");

        }
        else {
            // 2d
            t1 = System.currentTimeMillis();
            J = f.run2d_float(I, N, M, P, Vx, Vy, Vz);
            t2 = System.currentTimeMillis();
            IJ.log("t = " + IJ.d2s((t2 - t1) / 1000f,2) + " [sec]");
        }

        ImagePlus ip_vess   = new ImagePlus(imnameshort+"_Vess",  array2imagestack(J, N, M, P));

        if (saveoutput) {
            ImageConverter ic = new ImageConverter(ip_vess);
            ic.convertToGray8();
            IJ.saveAs(ip_vess, "TIFF", outdir + imnameshort + "_Vess.tif");
            IJ.log("export:\t\t" + outdir + imnameshort + "_Vess.tif");
        }
        else {
            ip_vess.show();
        }

        if (Vx!=null && Vy!=null && Vz!=null) { // directions

            //---------------------
            ImagePlus ip_vx     = new ImagePlus(imnameshort+"_Vx",  array2imagestack(Vx, N, M, P));
            if (saveoutput) {
                IJ.saveAs(ip_vx, "TIFF", outdir + imnameshort + "_Vx.tif");
                IJ.log("export:\t\t" + outdir + imnameshort + "_Vx.tif");
            }
            else { // if not saved, show only
                ip_vx.show();
            }

            //---------------------
            ImagePlus ip_vy     = new ImagePlus(imnameshort+"_Vy",  array2imagestack(Vx, N, M, P));
            if (saveoutput) {
                IJ.saveAs(ip_vy, "TIFF", outdir + imnameshort + "_Vy.tif");
                IJ.log("export:\t\t" + outdir + imnameshort + "_Vy.tif");
            }
            else {
                ip_vy.show();
            }

            //---------------------
            ImagePlus ip_vz     = new ImagePlus(imnameshort+"_Vz",  array2imagestack(Vz, N, M, P));
            if (saveoutput) {
                IJ.saveAs(ip_vz, "TIFF", outdir + imnameshort + "_Vz.tif");
                IJ.log("export:\t\t" + outdir + imnameshort + "_Vz.tif");
            }
            else {
                ip_vz.show();
            }

            //---------------------
            if (saveoutput) {

                String swcfilepath = outdir + imnameshort + "_Vxyz.swc";
                try {

                    // clean file before exporting
                    PrintWriter pw = new PrintWriter(swcfilepath);
                    pw.close();

                    PrintWriter logWriter = new PrintWriter(new BufferedWriter(new FileWriter(swcfilepath, true)));

                    logWriter.println("# input:\t" + imnameshort);
                    logWriter.println("# save local vessel direction vector V=(vx, vy, vz) each 2 consicutive lines of the swc format");
                    logWriter.println("# odd  line id: start voxel location (column, row, layer) = (x, y, z)");
                    logWriter.println("# even line id: direction (x+vx, y+vy, z+vz), vector saved in the followup comment");
                    logWriter.println("# https://bitbucket.org/miroslavradojevic/vess");
                    logWriter.println("# miro@braincadet.com");

                    int si = 0;

                    // export .swc file with the local eigen vectors
                    for (int z = 0; z < P; z++) { // go through the voxels
                        for (int x = 0; x < N; x++) {
                            for (int y = 0; y < M; y++) {

                                // extract float[3] vector from the Vx, Vy and Vz from corresponding byte images
                                // convert byte [-128,+127] to float [-1,+1]
                                float vx = Frangi.to_vec(Vx[z*N*M+y*N+x]);
                                float vy = Frangi.to_vec(Vy[z*N*M+y*N+x]);
                                float vz = Frangi.to_vec(Vz[z*N*M+y*N+x]);
                                float v = (float) Math.sqrt(Math.pow(vx,2)+Math.pow(vy,2)+Math.pow(vz,2));

                                if (v>Float.MIN_VALUE)  {
                                    vx /= v;
                                    vy /= v;
                                    vz /= v;
                                }
                                else {
                                    vx = 0;
                                    vy = 0;
                                    vz = 0;
                                }

                                if (J[z*N*M+y*N+x]>0.1) { // 0.05 Float.MIN_VALUE

                                    logWriter.println(
                                            (++si) + " " +
                                                    IJ.d2s(   6, 0) + " " + // color 6=,
                                                    IJ.d2s(x+.5, 1) + " " +
                                                    IJ.d2s(y+.5, 1) + " " +
                                                    IJ.d2s(z   , 1) + " " +
                                                    IJ.d2s(   1, 0) + " " + // radius=1
                                                    IJ.d2s(si+1,0) //
                                    ); // root

                                    logWriter.println(
                                            (++si) + " " +
                                                    IJ.d2s(   6, 0) + " " + // color 6=,
                                                    IJ.d2s(x+.5+  2 * vx, 5) + " " + // J[z*N*M+y*N+x]
                                                    IJ.d2s(y+.5+  2 * vy, 5) + " " +
                                                    IJ.d2s(z   +  2 * vz, 5) + " " +
                                                    IJ.d2s(1, 0) + " " + // radius=1
                                                    IJ.d2s(si-1,0)
                                                    +
                                                    "# vx=" + IJ.d2s(vx,5) +
                                                    "  vy=" + IJ.d2s(vy,5) +
                                                    "  vz=" + IJ.d2s(vz,5)
                                    );
                                }
                            }
                        }
                    }

                    logWriter.close();

                } catch (IOException e) {}

                IJ.log("export:\t\t"+ swcfilepath);

            }
        }
    } // run()

    private byte[] byteimage2array(ImagePlus ip) {

        if (ip.getType() != ImagePlus.GRAY8) {
            IJ.log("byteimage2array() : input image needs to be byte8.");
            return null;
        }

        int W = ip.getWidth();
        int H = ip.getHeight();
        int L = ip.getStackSize();

        byte[] img = new byte[W*H*L];

        for (int z = 0; z < L; z++) {
            byte[] ai = (byte[]) ip.getStack().getProcessor(z+1).getPixels();
            for (int j = 0; j < W*H; j++) {
                int x = j % W;
                int y = j / W;
                img[z*W*H+y*W+x] = ai[j];
            }
        }

        return img;

    }

    private float[] floatimage2array(ImagePlus ip) {

        if (ip.getType() != ImagePlus.GRAY32) { // ImagePlus.GRAY16
            IJ.log("floatimage2array() : input image needs to be float.");
            return null;
        }

        int W = ip.getWidth();
        int H = ip.getHeight();
        int L = ip.getStackSize();

        float[] img = new float[W*H*L];

        for (int z = 0; z < L; z++) {
            float[] ai = (float[]) ip.getStack().getProcessor(z+1).getPixels();
            for (int j = 0; j < W*H; j++) {
                int x = j % W;
                int y = j / W;
                img[z*W*H+y*W+x] = ai[j];
            }
        }

        return img;

    }

    private ImageStack array2imagestack(byte[] a, int W, int H, int L) {

        ImageStack is = new ImageStack(W, H);

        for (int z = 0; z < L; z++) {
            byte[] ai = new byte[W*H];
            for (int j = 0; j < W * H; j++) {
                int x = j % W;
                int y = j / W;
                ai[j] = a[z*W*H+y*W+x];
            }

            is.addSlice(new ByteProcessor(W, H, ai));

        }

        return is;

    }

    private ImageStack array2imagestack(float[] a, int W, int H, int L) {

        ImageStack is = new ImageStack(W, H);

        for (int z = 0; z < L; z++) {
            float[] ai = new float[W*H];
            for (int j = 0; j < W * H; j++) {
                int x = j % W;
                int y = j / W;
                ai[j] = a[z*W*H+y*W+x];
            }

            is.addSlice(new FloatProcessor(W, H, ai));

        }

        return is;

    }

    private static void createAndCleanDir(String dirpath) {
        File f1 = new File(dirpath);
        if (!f1.exists()) {f1.mkdirs();}
        else { // clean it, reset with each exec
            File[] files = f1.listFiles();
            for (File file : files)
                if (!file.delete()) System.out.println("Failed to delete " + file);
        }
    }

    private static void createDir(String dirpath) {
        // create directory without cleaning it up
        File f1 = new File(dirpath);
        if (!f1.exists()) {f1.mkdirs();}
    }

}
