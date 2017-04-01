package com.braincadet.vess;

import ij.IJ;

import java.util.ArrayList;

/**
 * <h1>Frangi's vesselness implementation.</h1>
 * <p> <i>Frangi, Alejandro F., et al. "Multiscale vessel enhancement filtering."
 * International Conference on Medical Image Computing and Computer-Assisted Intervention.
 * Springer Berlin Heidelberg, 1998.</i></p>
 * <p>The Frangi class runs Frangi's vesselness for selected byte8 2d/3d image.</p>
 * <p>Calculates the vesselness measure and the eigenvectors.</p>
 * <p>The implementation was inspired and uses MEX and Matlab samples ported from</p>
 * <a href="https://uk.mathworks.com/matlabcentral/fileexchange/24409-hessian-based-frangi-vesselness-filter?requestedDomain=www.mathworks.com">
 *     Hessian based Frangi vesselness filter</a>
 * @author  Miroslav Radojevic
 * @version 1.0
 * @since   2017-03-09
 */
public class Frangi {

    private static float Dstep = 2f/255f; // discretization step used for converting vector component float [-1,1] into 256 discretization steps stored in byte [-128,127]

    /**
     * {@value #sig} List of Gaussian sigma values in pixels.
     */
    public ArrayList<Float> sig = new ArrayList<Float>();
    /**
     * {@value #zdist} Anisotropy z-distance between layers, expressid in xy layer pixels
     */
    public float zdist;
    /**
     * {@value #alpha} Vesselness parameter (3d image stack vesselness)
     */
    public float alpha;
    /**
     * {@value #beta} Vesselness parameter (3d image stack vesselness)
     */
    public float beta; // 3d

    /**
     * {@value #C} Frangi vesselness computation constant (see doi: 10.1007/BFb0056195 and <a href="https://uk.mathworks.com/matlabcentral/fileexchange/24409-hessian-based-frangi-vesselness-filter?requestedDomain=www.mathworks.com">Hessian based Frangi vesselness filter</a>)
     */
    public float C;

    /**
     * {@value #BetaOne} Frangi vesselness computation constant used specially in 2d case (see doi and )
     */
    public float BetaOne;

    /**
     * {@value #BetaTwo} BetaTwo Frangi vesselness computation constant used specially in 2d case (see doi and )
     */
    public float BetaTwo;

    /**
     * {@value #blackwhite} Binary flag value. darkforeground=true (dark foreground, white background), darkforeground=false (white foreground, dark background).
     */
    public boolean blackwhite;

    /**
     * Class constructor: Frangi's vesselness computation.
     *
     * @param #sigs List of Gaussian sigma values expressing the scales for the vesselness computation.
     * @param zdist z-distance between image stack layers in pixels (3d, image stack parameter)
     * @param alpha Vesselness parameter (3d)
     * @param beta Vesselness parameter (3d)
     * @param C Vesselness parameter (3d)
     * @param beta_one Vesselness parameter (2d)
     * @param beta_two Vesselness parameter (2d)
     * @param blackwhite Vesselness parameter
     */
    public Frangi(ArrayList<Float> sigs, float zdist, float alpha, float beta, float C, float beta_one, float beta_two, boolean blackwhite) {

        this.zdist = zdist; // 3d
        this.alpha = alpha; // 3d
        this.beta = beta; // 3d
        this.C = C; // 3d
        this.BetaOne = beta_one; // 2d
        this.BetaTwo = beta_two; // 2d
        this.blackwhite = blackwhite;

        sig.clear();
        for (int i = 0; i < sigs.size(); i++) sig.add(sigs.get(i));

    }

    /**
     * Frangi's vesselness computation class constructor with the default parameter values.
     * <p>zdist = 2.0f (3d) z-dim anisotropy</p>
     * <p>alpha = 0.5f (3d)</p>
     * <p>beta = 0.5f (3d)</p>
     * <p>C = 500f (3d)</p>
     * <p>BetaOne = 0.5f (2d)</p>
     * <p>BetaTwo = 15f (2d)</p>
     * <p>darkforeground = false (black foreground, white background)</>
     *
     * @param sigs ArrayList<Float> with Gaussian cross section standard deviations corresponding to the scales
     */
    public Frangi(ArrayList<Float> sigs) {
        // initialize with the default values for the remainder of the parameters

        this.zdist = 2.0f; // 3d
        this.alpha = 0.5f; // 3d
        this.beta = 0.5f; // 3d
        this.C = 500f; // 3d
        this.BetaOne = 0.5f; // 2d
        this.BetaTwo = 15f; // 2d
        this.blackwhite = false;

        sig.clear();
        for (int i = 0; i < sigs.size(); i++) sig.add(sigs.get(i));
    }

    /**
     * Vesselness computation for image stack (3d).
     *
     * @param I Input image stack in 1-d byte array
     * @param w Image stack width
     * @param h Image stack height
     * @param l Image stack size
     * @param Vx Output image stack with eigen vector x-component. Byte values [0,255] range.
     * @param Vy
     * @param Vz
     * @param is_threaded Use threading flag.
     * @return Vesselness image stack saved in 1-d byte array [0,255]
     */
    public byte[] run3d_byte(
            byte[] I,
            int w,
            int h,
            int l,
            byte[] Vx,
            byte[] Vy,
            byte[] Vz,
            boolean is_threaded
    ) {

        if (is_threaded) {
            float[] J = new float[I.length];
            frangi3d(I, w, h, l, Vx, Vy, Vz, J); // vesselness saved in float[] J
            return get_vesselness(J);
        }
        else {
            float[] Jmin = new float[1];
            float[] Jmax = new float[1];
            float[] J = frangi3d(I, w, h, l, Vx, Vy, Vz, Jmin, Jmax);
            return get_vesselness(J, Jmin[0], Jmax[0]);
        }
    }

    /**
     * Vesselness computation for image stack (3d).
     *
     *
     * @param I
     * @param w
     * @param h
     * @param l
     * @param Vx
     * @param Vy
     * @param Vz
     * @param is_threaded
     * @return
     */
    public float[] run3d_float(
            byte[] I,
            int w,
            int h,
            int l,
            byte[] Vx,
            byte[] Vy,
            byte[] Vz,
            boolean is_threaded
    ){

        float[] J;

        if (is_threaded) {
            J = new float[I.length];
            frangi3d(I, w, h, l, Vx, Vy, Vz, J); // vesselness saved in float[] J
            min_max_norm(J, 0f, 1f); // [0,1] min-max normalized flat array representing the image stack
        }
        else {
            float[] Jmin = new float[1];
            float[] Jmax = new float[1];
            J = frangi3d(I, w, h, l, Vx, Vy, Vz, Jmin, Jmax);
            min_max_norm(J, Jmin[0], Jmax[0], 0f, 1f);
        }

        return J;

    }

    /**
     * Vesselness computation for image (2d).
     *
     * @param I Input image in 1-d byte array (accepts only Byte8 images)
     * @param w Image width
     * @param h Image height
     * @param l Image size
     * @param Vx Eigen vector x-component (skipped for Vx==null)
     * @param Vy Eigen vector y-component (skipped for Vy==null)
     * @param Vz Eigen vector z-component (skipped for Vz==null)
     * @return Vesselness image saved in 1-d byte array
     */
    public byte[] run2d_byte(
            byte[] I,
            int w,
            int h,
            int l,
            byte[] Vx,
            byte[] Vy,
            byte[] Vz
    ) {

        float[] Jmin = new float[1];
        float[] Jmax = new float[1];

        float[] J = frangi2d(I, w, h, l, Vx, Vy, Vz, Jmin, Jmax);

        return get_vesselness(J, Jmin[0], Jmax[0]);

    }

    /**
     * Vesselness computation for 2d image.
     *
     *
     * @param I
     * @param w
     * @param h
     * @param l
     * @param Vx
     * @param Vy
     * @param Vz
     * @return
     */
    public float[] run2d_float(
            byte[] I,
            int w,
            int h,
            int l,
            byte[] Vx,
            byte[] Vy,
            byte[] Vz
    ){

        float[] Jmin = new float[1];
        float[] Jmax = new float[1];

        float[] J = frangi2d(I, w, h, l, Vx, Vy, Vz, Jmin, Jmax);

        min_max_norm(J, Jmin[0], Jmax[0], 0f, 1f); // // min-max normalize J[] into [0,1]

        return J;

    }

    private float[] frangi2d(
            byte[] I,
            int w,
            int h,
            int l,
            byte[] Vx,
            byte[] Vy,
            byte[] Vz,
            float[] Jmin,
            float[] Jmax) {

        float[] J     = new float[w*h*l];
        float[] Dxx   = new float[w*h*l];
        float[] Dxy   = new float[w*h*l];
        float[] Dyy   = new float[w*h*l];

        float tmp, mag;
        float v2x, v2y, v1x, v1y;
        float mu1, mu2;
        float Lambda1, Lambda2;
        float Vecx, Vecy, Vecn;
        float Rb, S2, Ifiltered;
        int val;

        float beta  = (float) (2*Math.pow(BetaOne,2));
        float c     = (float) (2*Math.pow(BetaTwo,2));

        Jmin[0] = Float.POSITIVE_INFINITY;
        Jmax[0] = Float.NEGATIVE_INFINITY;

        for (int si = 0; si < sig.size(); ++si) {

            hessian2d(I, w, h, sig.get(si), Dyy, Dxy, Dxx);

            for (int i = 0; i < (w*h*l); ++i) {

                // compute the eigenvectors of J, v1 and v2
                tmp = (float) (Math.sqrt(Math.pow(Dxx[i]-Dyy[i],2) + 4*Math.pow(Dxy[i],2)));
                v2x = 2*Dxy[i];
                v2y = Dyy[i] - Dxx[i] + tmp;

                // normalize
                mag = (float) Math.sqrt(Math.pow(v2x,2) + Math.pow(v2y,2));
                if (mag>0) {
                    v2x /= mag;
                    v2y /= mag;
                }

                v1x = -v2y; // eigenvectors
                v1y = v2x;

                // eigenvalues
                mu1 = 0.5f * (Dxx[i] + Dyy[i] + tmp);
                mu2 = 0.5f * (Dxx[i] + Dyy[i] - tmp);

                //sort eigenvalues abs(Lambda1)>abs(Lambda2)
                boolean check = Math.abs(mu1)<Math.abs(mu2); // abs(mu1)>abs(mu2) was the original, it was swithched as the x-y axes were switched
                Lambda1 = (check)?mu2:mu1;
                Lambda2 = (check)?mu1:mu2;
                Vecx    = (check)?v2x:v1x;
                Vecy    = (check)?v2y:v1y;

                // extract vesselness measure using Lambda1 and Lambda2
                Lambda1 = (Lambda1==0)?Float.MIN_VALUE:Lambda1;
                Rb = (float)  Math.pow(Lambda2/Lambda1,2);
                S2 = (float) (Math.pow(Lambda1,2) + Math.pow(Lambda2,2));
                Ifiltered = (float) (Math.exp(-Rb/beta) * (1-Math.exp(-S2/c)));

                if (blackwhite)
                    Ifiltered = (Lambda1<0)?0:Ifiltered;
                else
                    Ifiltered = (Lambda1>0)?0:Ifiltered;

                // add result of this scale to output
                if (si==0) {

                    J[i] = Ifiltered;

                    if (J[i]<Jmin[0]) Jmin[0] = J[i];
                    if (J[i]>Jmax[0]) Jmax[0] = J[i];

                    Vecn = (float) Math.sqrt(Vecx*Vecx+Vecy*Vecy);

                    if (Vx != null && Vy!=null && Vz!=null) {

                        // convert [-1,+1] --> byte [-128,+127] (java's byte data type)

                        Vx[i] = to_byte(Vecx/Vecn); // Vecx/Vecn     --> Vx[i]
                        Vy[i] = to_byte(Vecy/Vecn); // Vecy/Vecn     --> Vy[i]
                        Vz[i] = to_byte(0.0);       // 0             --> Vz[i]

//                        val = (int) Math.round((((Vecx/Vecn)+1f)/2f)*255f);
//                        val = (val<0)?0:(val>255)?255:val; // int [0, 255]
//                        Vx[i] = (byte) val;

//                        val = (int) Math.round((((Vecy/Vecn)+1f)/2f)*255f);
//                        val = (val<0)?0:(val>255)?255:val; // int [0, 255]
//                        Vy[i] = (byte) val;

//                        val = (int) Math.round(((0f+1f)/2f)*255f);
////                        val = (val<0)?0:(val>255)?255:val;
//                        Vz[i] = (byte) val;

                    }

                }
                else {
                    if (Ifiltered>J[i]) {

                        J[i] = Ifiltered; // keep maximum filter response

                        if (J[i]<Jmin[0]) Jmin[0] = J[i];
                        if (J[i]>Jmax[0]) Jmax[0] = J[i];

                        Vecn = (float) Math.sqrt(Vecx*Vecx+Vecy*Vecy);

                        if (Vx != null && Vy != null && Vz != null) {

                            // convert [-1,+1] --> byte [-128,+127] (java's byte data type)

                            Vx[i] = to_byte(Vecx/Vecn); // Vecx/Vecn     --> Vx[i]
                            Vy[i] = to_byte(Vecy/Vecn); // Vecy/Vecn     --> Vy[i]
                            Vz[i] = to_byte(0.0);       // 0             --> Vz[i]

//                            val = (int) Math.round((((Vecx/Vecn)+1)/2)*255.0);
//                            val = (val<0)?0:(val>255)?255:val;
//                            Vx[i] = (byte) val;

//                            val = (int) Math.round((((Vecy/Vecn)+1)/2)*255.0);
//                            val = (val<0)?0:(val>255)?255:val;
//                            Vy[i] = (byte) val;

//                            Vz[i] = (byte) 0;

                        }

                    }
                }

            }

        }

        return J;
    }

    public static byte to_byte(double a) {
        // convert expected double value of the vector component [Vx, Vy, Vz] from [-1,+1] range and save into byte value [-128,+127]
        // precision loss Dstep/2 ~ 0.004
        return (byte) (clamp(Math.round(((float)a+1f)/Dstep), 0, 255) - 128);
    }

    public static float to_vec(byte a){
        // convert existing byte value [-128,127] used to save the vector component into float value [-1,+1]
        return ((int)a+128)*Dstep-1f;
    }

    private void float2byte(float[] I, byte[] J) {

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

    // 2d
    private void imgaussian(
            byte[] I,
            int w,
            int h,
            float sig,
            float[] F) {

        int i0, i1; // ok as long as w*h<Integer.MAX_VALUE, or w*h*l<Integer.MAX_VALUE (fulfilled on big neuron train set)

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
    }

    private void imgaussian_threaded(
            byte[] I,
            int w,
            int h,
            int l,
            float sig,
            float zdist,
            float[] F) {

        // use Gx, Gy, and Gz threaded3d classes for single axis convolutions
        // Gx gaussian smoothing of the byte8 image array stored in I[.], result in F[.]
        int CPU_NR = Runtime.getRuntime().availableProcessors();

        Gx.load(w, h, l, I, F, sig);
        Gx gx_jobs[] = new Gx[CPU_NR];
        for (int i = 0; i < gx_jobs.length; i++) {
            gx_jobs[i] = new Gx(i*(w/CPU_NR), (i+1)*(w/CPU_NR));
            gx_jobs[i].start();
        }
        for (int i = 0; i < gx_jobs.length; i++) {
            try {
                gx_jobs[i].join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

        // result should be in float[] F by now Gx.F and F are the same...
        float[] K = new float[w*h*l];
        Gy.load(w, h, l, F, K, sig);
        Gy gy_jobs[] = new Gy[CPU_NR];
        for (int i = 0; i < gy_jobs.length; i++) {
            gy_jobs[i] = new Gy(i*(h/CPU_NR), (i+1)*(h/CPU_NR));
            gy_jobs[i].start();
        }
        for (int i = 0; i < gy_jobs.length; i++) {
            try {
                gy_jobs[i].join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }


        Gz.load(w, h, l, K, F, sig, zdist);
        Gz gz_jobs[] = new Gz[CPU_NR];
        for (int i = 0; i < gz_jobs.length; i++) {
            gz_jobs[i] = new Gz(i*(l/CPU_NR), (i+1)*(l/CPU_NR));
            gz_jobs[i].start();
        }
        for (int i = 0; i < gz_jobs.length; i++) {
            try {
                gz_jobs[i].join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

    }

    private void imgaussian(
            byte[] I,
            int w,
            int h,
            int l,
            float sig,
            float zdist,
            float[] F) {

        int i0, i1;

        float sigz = sig/zdist;

        // gaussian filter is separated into 1D Gaussian kernels Gxy[.] and Gz[.] while z coordinate is scaled down
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
        for (int y = 0; y < h; ++y) {
            for (int z = 0; z < l; ++z) {
                // x index clamping
                for (int x = 0; x < Math.min(Lxy,w); ++x) {
                    i0 = z*w*h+y*w+x;
                    F[i0] = 0;
                    for (int x1 = x-Lxy; x1 <= x+Lxy; ++x1) {
                        i1 = z*w*h+y*w+clamp(x1,0,w-1);
                        F[i0] += (I[i1] & 0xff) * Gxy[x1-x+Lxy];
                    }
                }
                // no clamp
                for (int x = Lxy; x < (w-Lxy); ++x) {
                    i0 = z*w*h+y*w+x;
                    F[i0] = 0;
                    for (int x1 = x-Lxy; x1 <= x+Lxy; ++x1) {
                        i1 = z*w*h+y*w+x1;
                        F[i0] += (I[i1] & 0xff) * Gxy[x1-x+Lxy];
                    }
                }
                // x index clamping
                for (int x = Math.max(w-Lxy,Lxy); x < w; ++x) {
                    i0 = z*w*h+y*w+x;
                    F[i0] = 0;
                    for (int x1 = x-Lxy; x1 <= x+Lxy; ++x1) {
                        i1 = z*w*h+y*w+clamp(x1,0,w-1);
                        F[i0] += (I[i1] & 0xff) * Gxy[x1-x+Lxy];
                    }
                }

            }
        }

        // Gy gaussian smoothing of the float image array stored in F[.], result in K[.]
        float[] K = new float[w*h*l];
        for (int x = 0; x < w; ++x) {
            for (int z = 0; z < l; ++z) {
                // y index clamping
                for (int y = 0; y < Math.min(Lxy,h); ++y) {
                    i0 = z*w*h+y*w+x;
                    K[i0] = 0;
                    for (int y1 = y-Lxy; y1 <= y+Lxy; ++y1) {
                        i1 = z*w*h+clamp(y1,0,h-1)*w+x;
                        K[i0] += F[i1] * Gxy[y1-y+Lxy];
                    }
                }
                // no clamp
                for (int y = Lxy; y < (h-Lxy); ++y) {
                    i0 = z*w*h+y*w+x;
                    K[i0] = 0;
                    for (int y1 = (y-Lxy); y1 <= (y+Lxy); ++y1) {
                        i1 = z*w*h+y1*w+x;
                        K[i0] += F[i1] * Gxy[y1-y+Lxy];
                    }
                }
                // y index clamping
                for (int y = Math.max(h-Lxy,Lxy); y < h; ++y) {
                    i0 = z*w*h+y*w+x;
                    K[i0] = 0;
                    for (int y1 = y-Lxy; y1 <= y+Lxy; ++y1) {
                        i1 = z*w*h+clamp(y1,0,h-1)*w+x;
                        K[i0] += F[i1] * Gxy[y1-y+Lxy];
                    }
                }
            }
        }

        // Gz gaussian smoothing with sigz, smoothing of the float image array stored in K[.], result in F[.]
        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {
                // z index clamping
                for (int z = 0; z < Math.min(Lz,l); ++z) {
                    i0 = z*w*h+y*w+x;
                    F[i0] = 0;
                    for (int z1 = (z-Lz); z1 <= (z+Lz); ++z1) {
                        i1 = clamp(z1,0,l-1)*w*h+y*w+x;
                        F[i0] += K[i1] * Gz[z1-z+Lz];
                    }
                }
                // no clamp
                for (int z = Lz; z < (l-Lz); ++z) {
                    i0 = z*w*h+y*w+x;
                    F[i0] = 0;
                    for (int z1 = (z-Lz); z1 <= (z+Lz); ++z1) {
                        i1 = z1*w*h+y*w+x;
                        F[i0] += K[i1] * Gz[z1-z+Lz];
                    }
                }
                // z index clamping
                for (int z = Math.max(l-Lz,Lz); z < l; ++z) {
                    i0 = z*w*h+y*w+x;
                    F[i0] = 0;
                    for (int z1 = (z-Lz); z1 <= (z+Lz); ++z1) {
                        i1 = clamp(z1,0,l-1)*w*h+y*w+x;
                        F[i0] += K[i1] * Gz[z1-z+Lz];
                    }
                }

            }
        }

//        String F_out = "";
//        for (int i = 0; i < Math.min(F.length,5000); i+=100) { // every 100th
//            F_out += "[" + IJ.d2s(F[i], 3) + ", " + IJ.d2s(K[i], 3) + "] ";
//        }
//        IJ.log(F_out);

    }

    private void hessian3d(
            byte[] I,
            int w,
            int h,
            int l,
            float sig,
            float zdist,
            float[] Dzz,
            float[] Dyy,
            float[] Dyz,
            float[] Dxx,
            float[] Dxy,
            float[] Dxz,
            boolean is_threaded) {

        float[] F = new float[w*h*l];

        // faster gaussian filter convolutions for threaded3d filtering
        if (is_threaded) {
            imgaussian_threaded(I, w, h, l, sig, zdist, F);
        }
        else {
            imgaussian(I, w, h, l, sig, zdist, F);
        }

        int x,y,z;

        // todo: the derivations can be computed in paralel as well

        // Dz
        float[] DD = new float[w*h*l];
        for (int i = 0; i < (w*h*l); ++i) {

            x = i%w; z = i/(w*h); y = i/w-z*h;

            if (z==0)           DD[i] =      F[(z+1)*w*h+y*w+x]  - F[i];
            else if (z<(l-1))   DD[i] = .5f*(F[(z+1)*w*h+y*w+x]  - F[(z-1)*w*h+y*w+x]);
            else if (z==(l-1))  DD[i] =      F[i]                - F[(z-1)*w*h+y*w+x];
        }
        // Dzz
        for (int i = 0; i < (w*h*l); ++i) {
            x = i%w; z = i/(w*h); y = i/w-z*h;
            if (z==0)           Dzz[i] =      DD[(z+1)*w*h+y*w+x]  - DD[i];
            else if (z<(l-1))   Dzz[i] = .5f*(DD[(z+1)*w*h+y*w+x]  - DD[(z-1)*w*h+y*w+x]);
            else if (z==(l-1))  Dzz[i] =      DD[i]                - DD[(z-1)*w*h+y*w+x];

            Dzz[i] *= (sig*sig); // correct for scaling
        }
//        IJ.log("Dzz["+sig+"],");

        // Dy
        for (int i = 0; i < (w*h*l); ++i) {
            x = i%w; z = i/(w*h); y = i/w-z*h;
            if (y==0)            DD[i] =      F[z*w*h+(y+1)*w+x] - F[i];
            else if (y<(h-1))    DD[i] = .5f*(F[z*w*h+(y+1)*w+x]  - F[z*w*h+(y-1)*w+x]);
            else if (y==(h-1))   DD[i] =      F[i]                - F[z*w*h+(y-1)*w+x];
        }

        for (int i = 0; i < (w*h*l); ++i) {
            x = i%w; z = i/(w*h); y = i/w-z*h;
            // Dyy
            if (y==0)           Dyy[i] =      DD[z*w*h+(y+1)*w+x] - DD[i];
            else if (y<(h-1))   Dyy[i] = .5f*(DD[z*w*h+(y+1)*w+x] - DD[z*w*h+(y-1)*w+x]);
            else if (y==(h-1))  Dyy[i] =      DD[i]               - DD[z*w*h+(y-1)*w+x];

            Dyy[i] *= (sig*sig); // correct for scaling
            // Dyz
            if (z==0)           Dyz[i] =      DD[(z+1)*w*h+y*w+x] - DD[i];
            else if (z<(l-1))   Dyz[i] = .5f*(DD[(z+1)*w*h+y*w+x] - DD[(z-1)*w*h+y*w+x]);
            else if (z==(l-1))  Dyz[i] =      DD[i]               - DD[(z-1)*w*h+y*w+x];

            Dyz[i] *= (sig*sig); // correct for scaling
        }

//        IJ.log("Dyy[" + sig + "],");
//        IJ.log("Dyz[" + sig + "],");

        // Dx
        for (int i = 0; i < (w*h*l); ++i) {
            x = i%w; z = i/(w*h); y = i/w-z*h;
            if (x==0)           DD[i] =       F[z*w*h+y*w+(x+1)]  - F[i];
            else if (x<(w-1))   DD[i] =  .5f*(F[z*w*h+y*w+(x+1)]  - F[z*w*h+y*w+(x-1)]);
            else if (x==(w-1))  DD[i] =       F[i]                - F[z*w*h+y*w+(x-1)];
        }

        // remove F (garbage collector)

        for (int i = 0; i < (w*h*l); ++i) {
            x = i%w; z = i/(w*h); y = i/w-z*h;
            // Dxx
            if (x==0)           Dxx[i] =       DD[z*w*h+y*w+(x+1)] - DD[i];
            else if (x<(w-1))   Dxx[i] =  .5f*(DD[z*w*h+y*w+(x+1)] - DD[z*w*h+y*w+(x-1)]);
            else if (x==(w-1))  Dxx[i] =       DD[i]               - DD[z*w*h+y*w+(x-1)];

            Dxx[i] *= (sig*sig); // correct for scaling
            // Dxy
            if (y==0)           Dxy[i] =      DD[z*w*h+(y+1)*w+x] - DD[i];
            else if (y<(h-1))   Dxy[i] = .5f*(DD[z*w*h+(y+1)*w+x] - DD[z*w*h+(y-1)*w+x]);
            else if (y==(h-1))  Dxy[i] =      DD[i]               - DD[z*w*h+(y-1)*w+x];

            Dxy[i] *= (sig*sig); // correct for scaling
            // Dxz
            if (z==0)           Dxz[i] =      DD[(z+1)*w*h+y*w+x] - DD[i];
            else if (z<(l-1))   Dxz[i] = .5f*(DD[(z+1)*w*h+y*w+x] - DD[(z-1)*w*h+y*w+x]);
            else if (z==(l-1))  Dxz[i] =      DD[i]               - DD[(z-1)*w*h+y*w+x];

            Dxz[i] *= (sig*sig); // correct for scaling
        }

        // delete DD
//        IJ.log("Dxx["+sig+"],");
//        IJ.log("Dxy["+sig+"],");
//        IJ.log("Dxz["+sig+"]... ");// + ((t2-t1)/(double)CLOCKS_PER_SEC) + " sec.");

    }

    // calculate Dxx, Dxy, and D
    private void hessian2d(
            byte[] I,
            int w,
            int h,
            float sig,
            float[] Dyy,
            float[] Dxy,
            float[] Dxx
    ) {

        float[] F = new float[w*h];
        imgaussian(I, w, h, sig, F);

        int x,y;

        // Dy
        float[] DD = new float[w*h];
        for (int i = 0; i < (w*h); ++i) {
            x = i%w; y = i/w;
            if (y==0)            DD[i] =      F[(y+1)*w+x] - F[i];
            else if (y<(h-1))    DD[i] = .5f*(F[(y+1)*w+x] - F[(y-1)*w+x]);
            else if (y==(h-1))   DD[i] =      F[i]         - F[(y-1)*w+x];
        }

        for (int i = 0; i < (w*h); ++i) {
            x = i%w; y = i/w;
            // Dyy
            if (y==0)           Dyy[i] =      DD[(y+1)*w+x] - DD[i];
            else if (y<(h-1))   Dyy[i] = .5f*(DD[(y+1)*w+x] - DD[(y-1)*w+x]);
            else if (y==(h-1))  Dyy[i] =      DD[i]         - DD[(y-1)*w+x];

            Dyy[i] *= (sig*sig); // correct for scaling
        }

        // Dx
        for (int i = 0; i < (w*h); ++i) {
            x = i%w; y = i/w;
            if (x==0)           DD[i] =       F[y*w+(x+1)]  - F[i];
            else if (x<(w-1))   DD[i] =  .5f*(F[y*w+(x+1)]  - F[y*w+(x-1)]);
            else if (x==(w-1))  DD[i] =       F[i]          - F[y*w+(x-1)];
        }

//        delete [] F;   F = 0; // not necessary after the last axis derivative is done

        for (int i = 0; i < (w*h); ++i) {
            x = i%w; y = i/w;
            // Dxx
            if (x==0)           Dxx[i] =       DD[y*w+(x+1)] - DD[i];
            else if (x<(w-1))   Dxx[i] =  .5f*(DD[y*w+(x+1)] - DD[y*w+(x-1)]);
            else if (x==(w-1))  Dxx[i] =       DD[i]         - DD[y*w+(x-1)];

            Dxx[i] *= (sig*sig); // correct for scaling
            // Dxy
            if (y==0)           Dxy[i] =      DD[(y+1)*w+x] - DD[i];
            else if (y<(h-1))   Dxy[i] = .5f*(DD[(y+1)*w+x] - DD[(y-1)*w+x]);
            else if (y==(h-1))  Dxy[i] =      DD[i]         - DD[(y-1)*w+x];

            Dxy[i] *= (sig*sig); // correct for scaling
        }

//        delete [] DD;   DD = 0;

    }

    // frangi3d, non-threaded3d version
    private float[] frangi3d(
            byte[] I,
            int w,
            int h,
            int l,
            byte[] Vx,
            byte[] Vy,
            byte[] Vz,
            float[] Jmin,
            float[] Jmax
    ) {

        float[] J     = new float[w*h*l];
        float[] Dzz   = new float[w*h*l];
        float[] Dyy   = new float[w*h*l];
        float[] Dyz   = new float[w*h*l];
        float[] Dxx   = new float[w*h*l];
        float[] Dxy   = new float[w*h*l];
        float[] Dxz   = new float[w*h*l];

        double[][] Ma = new double[3][3];
        double[][] Davec = new double[3][3];
        double[]   Daeig = new double[3];

        double Lambda1, Lambda2, Lambda3;
        double LambdaAbs1, LambdaAbs2, LambdaAbs3;
        double Ra, Rb, S, expRa, expRb, expS;
        double Voxel_data;
        int val;

        Jmin[0] = Float.POSITIVE_INFINITY;
        Jmax[0] = Float.NEGATIVE_INFINITY;

        for (int si = 0; si < sig.size(); ++si) {

//            IJ.log("hessian3D[" + sig.get(si) + "]");
            hessian3d(I, w, h, l, sig.get(si), zdist, Dzz, Dyy, Dyz, Dxx, Dxy, Dxz, false); // not-threaded3d

//            IJ.log("eigen3D[" + sig.get(si) + "]");

//            int PLOT_INTERVAL = (w*h*l)/10;

            for (int i = 0; i < (w*h*l); ++i) {

//                if (i%PLOT_INTERVAL==0) IJ.log(((i/PLOT_INTERVAL)*10) + "%");

                Ma[0][0]=Dxx[i]; Ma[0][1]=Dxy[i]; Ma[0][2]=Dxz[i];
                Ma[1][0]=Dxy[i]; Ma[1][1]=Dyy[i]; Ma[1][2]=Dyz[i];
                Ma[2][0]=Dxz[i]; Ma[2][1]=Dyz[i]; Ma[2][2]=Dzz[i];

                eigen_decomposition(Ma, Davec, Daeig);

                Lambda1 = Daeig[0];
                Lambda2 = Daeig[1];
                Lambda3 = Daeig[2];
//            Vecx = Davec[0][0]; // vx
//            Vecy = Davec[1][0]; // vy
//            Vecz = Davec[2][0]; // vz
                LambdaAbs1 = Math.abs(Lambda1);
                LambdaAbs2 = Math.abs(Lambda2);
                LambdaAbs3 = Math.abs(Lambda3);

                Ra = LambdaAbs2/LambdaAbs3;
                Rb = LambdaAbs1/Math.sqrt(LambdaAbs2*LambdaAbs3);

                S = Math.sqrt(LambdaAbs1*LambdaAbs1+LambdaAbs2*LambdaAbs2+LambdaAbs3*LambdaAbs3);

                expRa = (1 - Math.exp(-((Ra*Ra)/(2*alpha*alpha))));
                expRb =      Math.exp(-((Rb*Rb)/(2*beta*beta)));
                expS  = (1 - Math.exp(-(S*S)/(2*C*C)));

                Voxel_data = expRa * expRb * expS;

                if (blackwhite) {
                    Voxel_data = (Lambda2<0)? 0 : Voxel_data;
                    Voxel_data = (Lambda3<0)? 0 : Voxel_data;
                }
                else {
                    Voxel_data = (Lambda2>0)? 0 : Voxel_data;
                    Voxel_data = (Lambda3>0)? 0 : Voxel_data;
                }

                // remove NaN
                Voxel_data = Double.isNaN(Voxel_data)? 0 : Voxel_data ;

                // add result of this scale to output
                if (si==0) {
                    J[i] = (float) Voxel_data;

                    if (J[i]<Jmin[0]) Jmin[0] = J[i]; // reference to Jmin in c++ is float[1] in java here
                    if (J[i]>Jmax[0]) Jmax[0] = J[i];

                    if (Vx != null && Vy != null && Vz != null) {

                        // convert [-1,+1] --> byte [-128,+127] (java's byte data type)

                        Vx[i] = to_byte(Davec[0][0]); // Davec[0][0]     --> Vx[i]
                        Vy[i] = to_byte(Davec[1][0]); // Davec[1][0]     --> Vy[i]
                        Vz[i] = to_byte(Davec[2][0]); // Davec[2][0]     --> Vz[i]

//                        val = (int) Math.round(  ((Davec[0][0]+1)/2) * 255 ); // vx
//                        val = (val<0)?0:(val>255)?255:val; // int [0, 255]
//                        Vx[i] = (byte) val;

//                        val = (int) Math.round(  ((Davec[1][0]+1)/2) * 255 ); // vy
//                        val = (val<0)?0:(val>255)?255:val; // int [0, 255]
//                        Vy[i] = (byte) val;

//                        val = (int) Math.round(  ((Davec[2][0]+1)/2) * 255 ); // vz
//                        val = (val<0)?0:(val>255)?255:val; // int [0, 255]
//                        Vz[i] = (byte) val;
                    }

                }
                else {
                    if (Voxel_data>J[i]) {
                        J[i] = (float) Voxel_data; // keep maximum filter response

                        if (J[i]<Jmin[0]) Jmin[0] = J[i];
                        if (J[i]>Jmax[0]) Jmax[0] = J[i];

                        if (Vx != null && Vy != null && Vz != null) {

                            // convert [-1,+1]  --> byte [-128,+127] (java's byte data type)

                            Vx[i] = to_byte(Davec[0][0]); // Davec[0][0]     --> Vx[i]
                            Vy[i] = to_byte(Davec[1][0]); // Davec[1][0]     --> Vy[i]
                            Vz[i] = to_byte(Davec[2][0]); // Davec[2][0]     --> Vz[i]

//                            val = (int) Math.round(  ((Davec[0][0]+1)/2) * 255 ); // vx
//                            val = (val<0)?0:(val>255)?255:val; // int [0, 255]
//                            Vx[i] = (byte) val;

//                            val = (int) Math.round(  ((Davec[1][0]+1)/2) * 255 ); // vy
//                            val = (val<0)?0:(val>255)?255:val; // int [0, 255]
//                            Vy[i] = (byte) val;

//                            val = (int) Math.round(  ((Davec[2][0]+1)/2) * 255 ); // vz
//                            val = (val<0)?0:(val>255)?255:val; // int [0, 255]
//                            Vz[i] = (byte) val;

                        }
                    }
                }
            } // pix
        } // sig

        return J;
    }

    // frangi3d, threaded3d version
    private void frangi3d(
            byte[] I,
            int w,
            int h,
            int l,
            byte[] Vx,
            byte[] Vy,
            byte[] Vz,
            float[] J
    ) {

        float[] Dzz   = new float[w*h*l];
        float[] Dyy   = new float[w*h*l];
        float[] Dyz   = new float[w*h*l];
        float[] Dxx   = new float[w*h*l];
        float[] Dxy   = new float[w*h*l];
        float[] Dxz   = new float[w*h*l];

        for (int si = 0; si < sig.size(); ++si) {

            hessian3d(I, w, h, l, sig.get(si), zdist, Dzz, Dyy, Dyz, Dxx, Dxy, Dxz, true); // threaded3d

            // threaded3d eigen value/vector computation

            int CPU_NR = Runtime.getRuntime().availableProcessors();

            Eigen3D eig_jobs[] = new Eigen3D[CPU_NR];

            for (int i = 0; i < eig_jobs.length; i++) {

//                IJ.log("[" + (i*(w*h*l)/CPU_NR)   + " -- " + ((i+1)*(w*h*l)/CPU_NR)   + "] " + (i*(w*h*l)>Integer.MAX_VALUE) + "," + ((i+1)*(w*h*l)>Integer.MAX_VALUE));
//                IJ.log("[" + (i*((w*h*l)/CPU_NR)) + " -- " + ((i+1)*((w*h*l)/CPU_NR)) + "] "+(w*h*l));

                eig_jobs[i] = new Eigen3D(
                        i    *((w*h*l)/CPU_NR),
                        (i+1)*((w*h*l)/CPU_NR),
                        Dxx, Dxy, Dxz, Dyy, Dyz, Dzz,
                        alpha,
                        beta,
                        C,
//                        BetaOne,
//                        BetaTwo,
                        blackwhite,
                        J, Vx, Vy, Vz,
                        si == 0);

                eig_jobs[i].start();
            }
            for (int i = 0; i < eig_jobs.length; i++) {
                try {
                    eig_jobs[i].join();
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }

        } // sigs

    }

    private static void min_max_norm(float[] J, float in_min, float in_max, float out_min, float out_max) {

        // min/max of the input array are precomputed
        for (int i = 0; i < J.length; i++) {
            J[i] = out_min + ( (J[i]-in_min)/(in_max-in_min) ) * (out_max - out_min); // [out_min, out_max]
        }

    }

    private static void min_max_norm(float[] J, float out_min, float out_max) {

        float Jmin = Float.POSITIVE_INFINITY;
        float Jmax = Float.NEGATIVE_INFINITY;

        for (int i = 0; i < J.length; i++) {
            if (J[i]<Jmin) Jmin = J[i];
            if (J[i]>Jmax) Jmax = J[i];
        }

        min_max_norm(J, Jmin, Jmax, out_min, out_max);

    }

    private static byte[] get_vesselness(float[] J) {

        byte[] J8 = new byte[J.length];

        float Jmin = Float.POSITIVE_INFINITY;
        float Jmax = Float.NEGATIVE_INFINITY;

        for (int i = 0; i < J.length; i++) {
            if (J[i]<Jmin) Jmin = J[i];
            if (J[i]>Jmax) Jmax = J[i];
        }

        // convert min-max normalized [0,255] J to J8 (byte8 image)
        if (Math.abs(Jmax-Jmin)<=Float.MIN_VALUE) {
            for (int i = 0; i < J.length; ++i) {
                J8[i] = (byte) 0;
            }
        }
        else {
            for (int i = 0; i < J.length; ++i) {
                int Jnorm = Math.round(((J[i]-Jmin)/(Jmax-Jmin)) * 255);
                Jnorm = (Jnorm<0)? 0 : (Jnorm > 255)? 255 : Jnorm; // wrap into [0,255] range
                J8[i] = (byte) Jnorm; // convert min-max normalized J to J8 (byte8 image)
            }
        }

        return J8;

    }

    private static byte[] get_vesselness(
            float[] J,
            float Jmin,
            float Jmax) {

        byte[] J8 = new byte[J.length]; // float[] J ---> byte[] J8

        // convert min-max normalized [0,255] J to J8 (byte8 image)
        if (Math.abs(Jmax-Jmin)<=Float.MIN_VALUE) {
            for (int i = 0; i < J.length; ++i) {
                J8[i] = (byte) 0;
            }
        }
        else {
            for (int i = 0; i < J.length; ++i) {
                int Jnorm = Math.round(((J[i]-Jmin)/(Jmax-Jmin)) * 255);
                Jnorm = (Jnorm<0)? 0 : (Jnorm > 255)? 255 : Jnorm; // int [0,255]
                J8[i] = (byte) Jnorm; // convert min-max normalized J to J8 (byte8 image)
            }
        }

        return J8;
    }

    private static void eigen_decomposition(double[][] A, double[][] V, double[] d) { // double[3][3] A, double[3][3] V, double[3] d

        int n = 3; // hessian dimension

//        double e[n];
        double[] e = new double[3];
//        double da[3];
        double[] da = new double[3];
        double dt, dat;
//        double vet[3];
        double[] vet = new double[3];
        int i, j;

        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                V[i][j] = A[i][j];
            }
        }

        tred2(V, d, e);
        tql2(V, d, e);

        /* Sort the eigen values and vectors by abs eigen value */
        da[0]=absd(d[0]); da[1]=absd(d[1]); da[2]=absd(d[2]);

        if((da[0]>=da[1])&&(da[0]>da[2])) {
            dt=d[2];   dat=da[2];    vet[0]=V[0][2];    vet[1]=V[1][2];    vet[2]=V[2][2];
            d[2]=d[0]; da[2]=da[0];  V[0][2] = V[0][0]; V[1][2] = V[1][0]; V[2][2] = V[2][0];
            d[0]=dt;   da[0]=dat;    V[0][0] = vet[0];  V[1][0] = vet[1];  V[2][0] = vet[2];
        }
        else if((da[1]>=da[0])&&(da[1]>da[2])) {
            dt=d[2];   dat=da[2];    vet[0]=V[0][2];    vet[1]=V[1][2];    vet[2]=V[2][2];
            d[2]=d[1]; da[2]=da[1];  V[0][2] = V[0][1]; V[1][2] = V[1][1]; V[2][2] = V[2][1];
            d[1]=dt;   da[1]=dat;    V[0][1] = vet[0];  V[1][1] = vet[1];  V[2][1] = vet[2];
        }

        if(da[0]>da[1]) {
            dt=d[1];   dat=da[1];    vet[0]=V[0][1];    vet[1]=V[1][1];    vet[2]=V[2][1];
            d[1]=d[0]; da[1]=da[0];  V[0][1] = V[0][0]; V[1][1] = V[1][0]; V[2][1] = V[2][0];
            d[0]=dt;   da[0]=dat;    V[0][0] = vet[0];  V[1][0] = vet[1];  V[2][0] = vet[2];
        }

    }

    // Symmetric Householder reduction to tridiagonal form.
    private static void tred2(double[][] V, double[] d, double[] e) {

        int n = 3;

        /*  This is derived from the Algol procedures tred2 by */
        /*  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for */
        /*  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding */
        /*  Fortran subroutine in EISPACK. */
        int i, j, k;
        double scale;
        double f, g, h;
        double hh;
        for (j = 0; j < n; j++) {d[j] = V[n-1][j]; }

        /* Householder reduction to tridiagonal form. */
        for (i = n-1; i > 0; i--) {
            /* Scale to avoid under/overflow. */
            scale = 0.0;
            h = 0.0;
            for (k = 0; k < i; k++) { scale = scale + Math.abs(d[k]); } // fabs() -- Math.abs()
            if (scale == 0.0) {
                e[i] = d[i-1];
                for (j = 0; j < i; j++) { d[j] = V[i-1][j]; V[i][j] = 0.0;  V[j][i] = 0.0; }
            } else {

                /* Generate Householder vector. */

                for (k = 0; k < i; k++) { d[k] /= scale; h += d[k] * d[k]; }
                f = d[i-1];
                g = Math.sqrt(h);
                if (f > 0) { g = -g; }
                e[i] = scale * g;
                h = h - f * g;
                d[i-1] = f - g;
                for (j = 0; j < i; j++) { e[j] = 0.0; }

                /* Apply similarity transformation to remaining columns. */

                for (j = 0; j < i; j++) {
                    f = d[j];
                    V[j][i] = f;
                    g = e[j] + V[j][j] * f;
                    for (k = j+1; k <= i-1; k++) { g += V[k][j] * d[k]; e[k] += V[k][j] * f; }
                    e[j] = g;
                }
                f = 0.0;
                for (j = 0; j < i; j++) { e[j] /= h; f += e[j] * d[j]; }
                hh = f / (h + h);
                for (j = 0; j < i; j++) { e[j] -= hh * d[j]; }
                for (j = 0; j < i; j++) {
                    f = d[j]; g = e[j];
                    for (k = j; k <= i-1; k++) { V[k][j] -= (f * e[k] + g * d[k]); }
                    d[j] = V[i-1][j];
                    V[i][j] = 0.0;
                }
            }
            d[i] = h;
        }

        /* Accumulate transformations. */

        for (i = 0; i < n-1; i++) {
            V[n-1][i] = V[i][i];
            V[i][i] = 1.0;
            h = d[i+1];
            if (h != 0.0) {
                for (k = 0; k <= i; k++) { d[k] = V[k][i+1] / h;}
                for (j = 0; j <= i; j++) {
                    g = 0.0;
                    for (k = 0; k <= i; k++) { g += V[k][i+1] * V[k][j]; }
                    for (k = 0; k <= i; k++) { V[k][j] -= g * d[k]; }
                }
            }
            for (k = 0; k <= i; k++) { V[k][i+1] = 0.0;}
        }
        for (j = 0; j < n; j++) { d[j] = V[n-1][j]; V[n-1][j] = 0.0; }
        V[n-1][n-1] = 1.0;
        e[0] = 0.0;

    }

    // Symmetric tridiagonal QL algorithm.
    private static void tql2(double[][] V, double[] d, double[] e) {

        int n = 3;

        /*  This is derived from the Algol procedures tql2, by */
        /*  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for */
        /*  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding */
        /*  Fortran subroutine in EISPACK. */

        int i, j, k, l, m;
        double f;
        double tst1;
        double eps;
        int iter;
        double g, p, r;
        double dl1, h, c, c2, c3, el1, s, s2;

        for (i = 1; i < n; i++) { e[i-1] = e[i]; }
        e[n-1] = 0.0;

        f = 0.0;
        tst1 = 0.0;
        eps = Math.pow(2.0, -52.0);
        for (l = 0; l < n; l++) {

        /* Find small subdiagonal element */

            tst1 = MAX(tst1, Math.abs(d[l]) + Math.abs(e[l])); // fabs() -- Math.abs()
            m = l;
            while (m < n) {
                if (Math.abs(e[m]) <= eps*tst1) { break; } // fabs() -- Math.abs()
                m++;
            }

        /* If m == l, d[l] is an eigenvalue, */
        /* otherwise, iterate. */

            if (m > l) {
                iter = 0;
                do {
                    iter = iter + 1;  /* (Could check iteration count here.) */
                /* Compute implicit shift */
                    g = d[l];
                    p = (d[l+1] - g) / (2.0 * e[l]);
                    r = hypot2(p, 1.0);
                    if (p < 0) { r = -r; }
                    d[l] = e[l] / (p + r);
                    d[l+1] = e[l] * (p + r);
                    dl1 = d[l+1];
                    h = g - d[l];
                    for (i = l+2; i < n; i++) { d[i] -= h; }
                    f = f + h;
                /* Implicit QL transformation. */
                    p = d[m]; c = 1.0; c2 = c; c3 = c;
                    el1 = e[l+1]; s = 0.0; s2 = 0.0;
                    for (i = m-1; i >= l; i--) {
                        c3 = c2;
                        c2 = c;
                        s2 = s;
                        g = c * e[i];
                        h = c * p;
                        r = hypot2(p, e[i]);
                        e[i+1] = s * r;
                        s = e[i] / r;
                        c = p / r;
                        p = c * d[i] - s * g;
                        d[i+1] = h + s * (c * g + s * d[i]);
                    /* Accumulate transformation. */
                        for (k = 0; k < n; k++) {
                            h = V[k][i+1];
                            V[k][i+1] = s * V[k][i] + c * h;
                            V[k][i] = c * V[k][i] - s * h;
                        }
                    }
                    p = -s * s2 * c3 * el1 * e[l] / dl1;
                    e[l] = s * p;
                    d[l] = c * p;

                /* Check for convergence. */
                } while (Math.abs(e[l]) > eps*tst1); // fabs() -- Math.abs()
            }
            d[l] = d[l] + f;
            e[l] = 0.0;
        }

    /* Sort eigenvalues and corresponding vectors. */
        for (i = 0; i < n-1; i++) {
            k = i;
            p = d[i];
            for (j = i+1; j < n; j++) {
                if (d[j] < p) {
                    k = j;
                    p = d[j];
                }
            }
            if (k != i) {
                d[k] = d[i];
                d[i] = p;
                for (j = 0; j < n; j++) {
                    p = V[j][i];
                    V[j][i] = V[j][k];
                    V[j][k] = p;
                }
            }
        }
    }

    private static double hypot2(double x, double y) { return Math.sqrt(x*x+y*y); }

    private static double absd(double val){ if(val>0){ return val;} else { return -val;} }

    private static int clamp(int x, int x1, int x2) {
        int xC = (x<x1)?x1:x;
        return (xC>x2)?x2:xC;
    }

    private static double MAX(double a, double b) {
        return (a>b)?a:b;
    }

}
