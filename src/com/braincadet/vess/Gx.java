package com.braincadet.vess;

import ij.IJ;

/**
 * <h1>Threaded implementation of image Gaussian filtering along x-dim (image width).</h1>
 * Filtering along single dimension is used as a stage in the steered filtering implementation filtering sequence.
 * Gx class runs Gaussian filtering along x-dim for selected byte8 2d/3d image.
 * Implemented as thread class which splits computation between the available CPUs, column index (x coordinate) is used as the threading index
 * <p>Gx, Gy, Gz usage:</p>
 * * <pre>
 * {@code
 * int CPU_NR = Runtime.getRuntime().availableProcessors();
 * Gx.load(w, h, l, I, F, sig);
 * Gx gx_jobs[] = new Gx[CPU_NR];
 * for (int i = 0; i < gx_jobs.length; i++) {
 * gx_jobs[i] = new Gx(i*w/CPU_NR, (i+1)*w/CPU_NR);
 * gx_jobs[i].start();
 * }
 * for (int i = 0; i < gx_jobs.length; i++) {
 * try {
 * gx_jobs[i].join();
 * } catch (InterruptedException e) {
 * e.printStackTrace();
 * }
 * }
 * // result should be in float[] F by now Gx.F and F are the same...
 * float[] K = new float[w*h*l];
 * Gy.load(w, h, l, F, K, sig);
 * Gy gy_jobs[] = new Gy[CPU_NR];
 * for (int i = 0; i < gy_jobs.length; i++) {
 * gy_jobs[i] = new Gy(i*h/CPU_NR, (i+1)*h/CPU_NR);
 * gy_jobs[i].start();
 * }
 * for (int i = 0; i < gy_jobs.length; i++) {
 * try {
 * gy_jobs[i].join();
 * } catch (InterruptedException e) {
 * e.printStackTrace();
 * }
 * }
 * Gz.load(w, h, l, K, F, sig, zdist);
 * Gz gz_jobs[] = new Gz[CPU_NR];
 * for (int i = 0; i < gz_jobs.length; i++) {
 * gz_jobs[i] = new Gz(i*l/CPU_NR, (i+1)*l/CPU_NR);
 * gz_jobs[i].start();
 * }
 * for (int i = 0; i < gz_jobs.length; i++) {
 * try {
 * gz_jobs[i].join();
 * } catch (InterruptedException e) {
 * e.printStackTrace();
 * }
 * }
 * }
 * </pre>
 *
 * @author  Miroslav Radojevic
 * @version 1.0
 * @since   2017-03-11
 */
public class Gx extends Thread {

    private int n0, n1; // range for CPU threading

    private static int w;
    private static int h;
    private static int l;

    private static byte[]  I; // input
    private static float[] F; // output

    private static float sig;
    private static int Lxy;
    private static float[] Gxy; // filtering kernel

    /**
     * Thread constructor paralelization interval.
     * @param n0 start width
     * @param n1 end width
     */
    public Gx(int n0, int n1) { // constructor defines threading range
        this.n0 = n0;
        this.n1 = n1;
    }

    /**
     * Initialize Gx static components before starting the parallel jobs.
     *
     * @param w1 Image stack width
     * @param h1 Image stack height
     * @param l1 Image stack size
     * @param I1 Input image stack as 1-d byte array
     * @param F1 Gaussian filtered image stack along x-dimension, stored as 1-d float array
     * @param sig1 Gaussian standard deviation used for the filter (corresponds to the scale)
     */
    public static void load(int w1, int h1, int l1, byte[] I1, float[] F1, float sig1) {

        w = w1;
        h = h1;
        l = l1;

        I = I1;
        F = F1;

        sig = sig1;
        Lxy = (int) Math.ceil(3 * sig);
        Gxy = new float[2*Lxy+1];

        float Gnorm = 0;
        for (int i = -Lxy; i <= Lxy; ++i) {
            Gxy[i+Lxy] = (float) Math.exp(-(i*i)/(2*sig*sig));
            Gnorm += Gxy[i+Lxy];
        }

        for (int i = 0; i < (2*Lxy+1); ++i) {
            Gxy[i] /= Gnorm;
        }

    }

    /**
     * run Thread method
     */
    public void run() {

        int i0, i1;

        for (int x = n0; x < n1; x++) { // n0, n1 are within [0, w)
            for (int y = 0; y < h; y++) {
                for (int z = 0; z < l; z++) {

                    i0 = z*w*h+y*w+x;
                    F[i0] = 0;

                    for (int x1 = x-Lxy; x1 <= x+Lxy; ++x1) {

                        if (x>=0 && x<Math.min(Lxy,w))          i1 = clamp(x1, 0, w - 1);
                        else if (x>=Lxy && x<(w-Lxy))           i1 = x1;
                        else if (x>=Math.max(w-Lxy,Lxy) && x<w) i1 = clamp(x1, 0, w - 1);
                        else                                   { i1 = Integer.MIN_VALUE;
                            IJ.log("illegal case happens!");} // this is illegal case

                        i1 = i1 + (z*w*h+y*w);

                        F[i0] += (I[i1] & 0xff) * Gxy[x1-x+Lxy];
                    }
                }
            }
        }
    }

    private static int clamp(int x, int x1, int x2) {
        int xC = (x<x1)?x1:x;
        return (xC>x2)?x2:xC;
    }

}
