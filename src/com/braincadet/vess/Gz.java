package com.braincadet.vess;

/**
 * <h1>Threaded implementation of image Gaussian filtering along z-dim (image stack size).</h1>
 * Filtering along single dimension is used as a stage in the steered filtering implementation filtering sequence.
 * Gz class runs Gaussian filtering along z-dim for selected byte8 2d/3d image.
 * Implemented as thread class which splits computation between the available CPUs, stack layer index (z coordinate) is used as the threading index
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
public class Gz extends Thread {

    private int n0, n1; // range for CPU threading

    private static int w;
    private static int h;
    private static int l;

    private static float[] K; // input
    private static float[] F; // output

    private static float sigz;
    private static int Lz;
    private static float[] Gz;

    /**
     * Thread constructor paralelization interval.
     * @param n0 start layer
     * @param n1 end layer
     */
    public Gz(int n0, int n1) { // constructor defines threading range
        this.n0 = n0;
        this.n1 = n1;
    }

    /**
     * Initialize Gz static components before starting the parallel jobs.
     *
     * @param w1 Image stack width
     * @param h1 Image stack height
     * @param l1 Image stack size
     * @param K1 Input image stack (Gx,Gy filtered) as 1-d float array
     * @param F1 Gaussian filtered image stack stored as 1-d float array
     * @param sig1 Gaussian standard deviation used for the filter (corresponds to the scale)
     * @param zdist z-distance, z resolution in pixels (xy)
     */
    public static void load(int w1, int h1, int l1, float[] K1, float[] F1, float sig1, float zdist) {

        w = w1;
        h = h1;
        l = l1;

        K = K1;
        F = F1;

        sigz = sig1/zdist;
        Lz  = (int) Math.ceil(3 * sigz);
        Gz = new float[2*Lz+1];

        float Gnorm = 0;
        for (int i = -Lz; i <= Lz; ++i) {
            Gz[i+Lz] = (float) Math.exp(-(i*i)/(2*sigz*sigz));
            Gnorm += Gz[i+Lz];
        }

        for (int i = 0; i < (2*Lz+1); ++i) {
            Gz[i] /= Gnorm;
        }

    }

    /**
     * run Thread method
     */
    public void run() {

        int i0, i1;

        for (int z = n0; z < n1; z++) { // n0, n1 are within [0, l)
            for (int x = 0; x < w; ++x) {
                for (int y = 0; y < h; ++y) {

                    i0 = z*w*h+y*w+x;
                    F[i0] = 0;

                    for (int z1 = (z-Lz); z1 <= (z+Lz); ++z1) {

                        if (z>=0 && z<Math.min(Lz,l))           i1 = clamp(z1,0,l-1)*w*h;
                        else if (z>=Lz && z<(l-Lz))             i1 = z1             *w*h;
                        else if (z>=Math.max(l-Lz,Lz) && z<l)   i1 = clamp(z1,0,l-1)*w*h;
                        else                                    i1 = Integer.MIN_VALUE;  // this is illegal case

                        i1 = i1 + (y*w+x);

                        F[i0] += K[i1] * Gz[z1-z+Lz];

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
