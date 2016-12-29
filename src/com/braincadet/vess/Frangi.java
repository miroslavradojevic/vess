package com.braincadet.vess;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;

import java.util.ArrayList;

public class Frangi {

    ArrayList<Float> sig = new ArrayList<Float>();
    float zdist;
    float alpha, beta; // 3d
    float BetaOne, BetaTwo;
    float C;
    boolean blackwhite;

    public Frangi(ArrayList<Float> _sigs, float _zdist, float _alpha, float _beta, float _C, float _beta_one, float _beta_two) {

        alpha = _alpha;
        beta = _beta;
        C = _C;
        BetaOne = _beta_one;
        BetaTwo = _beta_two;
        zdist = _zdist;

        sig.clear();
        for (int i = 0; i < _sigs.size(); i++) {
            sig.add(_sigs.get(i));
        }

        blackwhite = false;
    }

    public void test(int[] a) {
        for (int i = 0; i < a.length; i++) {
            a[i] = i;
        }
    }

    public void pr(int[] a) {
        for (int i = 0; i < a.length; i++) {
            IJ.log("a["+i+"] = " + a[i]);
        }
    }

    public byte[] ip2array(ImagePlus ip) {

        IJ.log("ip2array...");

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

        IJ.log("done");

        return img;

    }

    public void float2byte(float[] I, byte[] J) {

        IJ.log("float2byte...");

        float Imin = Float.POSITIVE_INFINITY, Imax=Float.NEGATIVE_INFINITY;

        for (int i = 0; i < I.length; i++) {
            if (I[i]<Imin) Imin = I[i];
            if (I[i]>Imax) Imax = I[i];
        }

        for (int i = 0; i < I.length; i++) {
            int val = Math.round(((I[i]-Imin)/(Imax-Imin)) * 255f);
            J[i] = (byte) ((val<0)? 0 : ((val>255)? 255 : val));
        }

        IJ.log("done");

    }

    public ImageStack array2imagestack(byte[] a, int W, int H, int L) {

        IJ.log("array2imagestack...");

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

        IJ.log("done.");

        return is;

    }

    public ImageStack array2imagestack(float[] a, int W, int H, int L) {

        IJ.log("array2imagestack...");

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

        IJ.log("done.");

        return is;

    }

    public void imgaussian(byte[] I, int w, int h, float sig, float[] F) {

        IJ.log("imgaussian...");

//        long i0, i1;
        int i0, i1; // ok as long as w*h<Integer.MAX_VALUE

        // gaussian filter is separated into 1D Gaussian kernel Gxy[.] that will be used for filtering along x and y
        int Lxy = (int) Math.ceil(3*sig);
        float[] Gxy = new float[2*Lxy+1];
        float Gnorm = 0;

        for (int i = -Lxy; i <= Lxy; i++) {
            Gxy[i+Lxy] = (float) Math.exp(-(i*i)/(2*sig*sig));
            Gnorm += Gxy[i+Lxy];
        }

        for (int i = 0; i < (2*Lxy+1); i++) {
            Gxy[i] /= Gnorm;
        }

        float[] K = new float[w*h];

        // Gx gaussian smoothing of the byte8 image array stored in I[.], result in K[.]
        for (int y = 0; y < h; ++y) {
            // x index clamping
            for (int x = 0; x < Math.min(Lxy,w); ++x) {
                i0 = y*w+x;
                K[i0] = 0;
                for (int x1 = x-Lxy; x1 <= x+Lxy; ++x1) {
                    i1 = y*w+clamp(x1,0,w-1);
                    K[i0] += (I[i1] & 0xff) * Gxy[x1-x+Lxy];
                }
            }
            // no clamp
            for (int x = Lxy; x < (w-Lxy); ++x) {
                i0 = y*w+x;
                K[i0] = 0;
                for (int x1 = x-Lxy; x1 <= x+Lxy; ++x1) {
                    i1 = y*w+x1;
                    K[i0] += (I[i1] & 0xff) * Gxy[x1-x+Lxy];
                }
            }
            // x index clamping
            for (int x = Math.max(w-Lxy,Lxy); x < w; ++x) {
                i0 = y*w+x;
                K[i0] = 0;
                for (int x1 = x-Lxy; x1 <= x+Lxy; ++x1) {
                    i1 = y*w+clamp(x1,0,w-1);
                    K[i0] += (I[i1] & 0xff) * Gxy[x1-x+Lxy];
                }
            }
        }

        // Gy gaussian smoothing of the float image array stored in K[.], result in F[.]
        for (int x = 0; x < w; ++x) {
            // y index clamping
            for (int y = 0; y < Math.min(Lxy,h); ++y) {
                i0 = y*w+x;
                F[i0] = 0;
                for (int y1 = y-Lxy; y1 <= y+Lxy; ++y1) {
                    i1 = clamp(y1,0,h-1)*w+x;
                    F[i0] += K[i1] * Gxy[y1-y+Lxy];
                }
            }
            // no clamp
            for (int y = Lxy; y < (h-Lxy); ++y) {
                i0 = y*w+x;
                F[i0] = 0;
                for (int y1 = (y-Lxy); y1 <= (y+Lxy); ++y1) {
                    i1 = y1*w+x;
                    F[i0] += K[i1] * Gxy[y1-y+Lxy];
                }
            }
            // y index clamping
            for (int y = Math.max(h-Lxy,Lxy); y < h; ++y) {
                i0 = y*w+x;
                F[i0] = 0;
                for (int y1 = y-Lxy; y1 <= y+Lxy; ++y1) {
                    i1 = clamp(y1,0,h-1)*w+x;
                    F[i0] += K[i1] * Gxy[y1-y+Lxy];
                }
            }
        }
//        delete [] K; K = 0;

        IJ.log("done");

    }

    public void imgaussian(byte[] I, int w, int h, int l, float sig, float zdist, float[] F) {

        int i0, i1;

        float sigz = sig/zdist;

        // gaussian filter is separated into 1D Gaussian kernels Gxy[.] and Gz[.] if z coordinate is scaled
        int Lxy = (int) Math.ceil(3*sig);
        float[] Gxy = new float[2*Lxy+1];

        float Gnorm = 0;
        for (int i = -Lxy; i <= Lxy; ++i) {
            Gxy[i+Lxy] = (float) Math.exp(-(i*i)/(2*sig*sig));
            Gnorm += Gxy[i+Lxy];
        }

        for (int i = 0; i < (2*Lxy+1); ++i) {
            Gxy[i] /= Gnorm;
        }

        Gnorm = 0;
        int Lz  = (int) Math.ceil(3*sigz);
        float[] Gz = new float[2*Lz+1];
        for (int i = -Lz; i <= Lz; ++i) {
            Gz[i+Lz] = (float) Math.exp(-(i*i)/(2*sigz*sigz));
            Gnorm += Gz[i+Lz];
        }

        for (int i = 0; i < (2*Lz+1); ++i) {
            Gz[i] /= Gnorm;
        }

        // Gx gaussian smoothing of the byte8 image array stored in I[.], result in F[.]


    }

//    public void frangi3d(byte[] I, int w, int h, int l) {
//    }

    static int clamp(int x, int x1, int x2) {
        int xC = (x<x1)?x1:x;
        return (xC>x2)?x2:xC;
    }

}
