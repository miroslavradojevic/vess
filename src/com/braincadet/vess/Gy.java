package com.braincadet.vess;

import ij.IJ;

/**
 * <h1>Threaded implementation of image Gaussian filtering along y-dim (image height).</h1>
 * Filtering along single dimension is used as a stage in the steered filtering implementation filtering sequence.
 * Gy class runs Gaussian filtering along y-dim for selected byte8 2d/3d image.
 * Implemented as thread class which splits computation between the available CPUs, column index (y coordinate) is used as the threading index
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
public class Gy extends Thread {

    private int n0, n1; // range for CPU threading

    private static int w;
    private static int h;
    private static int l;

    private static float[] F; // input
    private static float[] K; // output

    private static float sig;
    private static int Lxy;
    private static float[] Gxy;

    /**
     * Thread constructor paralelization interval.
     * @param n0 start height
     * @param n1 end height
     */
    public Gy(int n0, int n1) { // constructor defines threading range
        this.n0 = n0;
        this.n1 = n1;
    }

    /**
     * Initialize Gy static components before starting the parallel jobs.
     *
     * @param w1 Image stack width
     * @param h1 Image stack height
     * @param l1 Image stack size
     * @param F1 Input image stack (Gx filtered) as 1-d float array
     * @param K1 Gaussian filtered image stack along y-dimension, stored as 1-d float array
     * @param sig1 Gaussian standard deviation used for the filter (corresponds to the scale)
     */
    public static void load(int w1, int h1, int l1, float[] F1, float[] K1, float sig1) {

        w = w1;
        h = h1;
        l = l1;

        F = F1;
        K = K1;

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

        for (int y = n0; y < n1; y++) { // n0, n1 are within [0, w)
            for (int x = 0; x < w; x++) {
                for (int z = 0; z < l; z++) {

                    i0 = z*w*h+y*w+x;
                    K[i0] = 0;

                    for (int y1 = y-Lxy; y1 <= y+Lxy; ++y1) {

                        if (y>=0 && y<Math.min(Lxy,h))          i1 = clamp(y1,0,h-1)*w;
                        else if (y>=Lxy && y<(h-Lxy))           i1 = y1             *w;
                        else if (y>=Math.max(h-Lxy,Lxy) && y<h) i1 = clamp(y1,0,h-1)*w;
                        else                                    {
                            IJ.log("SHOULD NOT happen. y="+y+" [n0="+n0+",n1="+n1+"], Lxy="+Lxy+", h=" +h+" | " + Math.min(Lxy,h)+"  |  " + Math.max(h-Lxy,Lxy));
                            i1 = Integer.MIN_VALUE; // this is illegal case
                        }

                        i1 = i1 + (z*w*h+x);

                        K[i0] += F[i1] * Gxy[y1-y+Lxy];

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
