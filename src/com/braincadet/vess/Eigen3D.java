package com.braincadet.vess;

/**
 * <h1>Threaded implementation of the Eigen analysis.</h1>
 * Takes the Hessian derivatives of the 3d image stack as input,
 * Dxx, Dxy, Dxz; Dyy, Dyz; Dzz
 * and computes
 * eigen vectors Vx, Vy, Vz  (saved as output)
 * eigen values lambda 1,2,3 (not saved as output)
 * and the vesselness value at each voxel (saved as output, J)
 * using the eigen values and following the recipe from (see doi: 10.1007/BFb0056195 and <a href="https://uk.mathworks.com/matlabcentral/fileexchange/24409-hessian-based-frangi-vesselness-filter?requestedDomain=www.mathworks.com">
 *     Hessian based Frangi vesselness filter</a>).
 *
 * * <pre>
 * {@code
 *
 * int CPU_NR = Runtime.getRuntime().availableProcessors();
 * Eigen3D eig_jobs[] = new Eigen3D[CPU_NR];
 * for (int i = 0; i < eig_jobs.length; i++) {
 * eig_jobs[i] = new Eigen3D(i*(w*h*l)/CPU_NR, (i+1)*(w*h*l)/CPU_NR, Dxx, Dxy, Dxz, Dyy, Dyz, Dzz, alpha, beta, C, blackwhite, J, Vx, Vy, Vz, si == 0);
 * eig_jobs[i].start();
 * }
 * for (int i = 0; i < eig_jobs.length; i++) {
 * try {
 * eig_jobs[i].join();
 * } catch (InterruptedException e) {
 * e.printStackTrace();
 * }
 * }
 * </pre>
 */
public class Eigen3D extends Thread {

    private int n0, n1; // range for CPU threading

    private static float[] Dzz; // Hessian derivative images
    private static float[] Dyz;
    private static float[] Dyy;
    private static float[] Dxz;
    private static float[] Dxy;
    private static float[] Dxx;

    private float alpha; // vesselness parameters
    private float beta;
    private float C;
//    private float BetaOne; // these are parametrizing 2d case which is not threaded3d
//    private float BetaTwo;
    private boolean blackwhite;

    /**
     * {@value #J} 1d array image with Frangi's vesselness.
     */
    public static float[] J;  // Frangi's vesselness
    /**
     * {@value #Vx} 1d array image with output eigenvector x-component saved as byte.
     */
    public static byte[] Vx;
    /**
     * {@value #Vy} 1d array image with output eigenvector y-component saved as byte.
     */
    public static byte[] Vy;
    /**
     * {@value #Vz} 1d array image with output eigenvector z-component saved as byte.
     */
    public static byte[] Vz;

    private boolean enforce;

    /**
     * Thread class constructor.
     * Initialize with chosen paralelization interval and the Hessian derivative inputs.
     *
     * @param n0 Starting threading interval index.
     * @param n1 Ending threading interval index.
     * @param Dxx Hessian derivative Dxx component 1d array image.
     * @param Dxy Hessian derivative component 1d array image.
     * @param Dxz Hessian derivative component 1d array image.
     * @param Dyy Hessian derivative component 1d array image.
     * @param Dyz Hessian derivative component 1d array image.
     * @param Dzz Hessian derivative component 1d array image.
     */
    public Eigen3D(int n0, int n1,
                   float[] Dxx, float[] Dxy, float[] Dxz, float[] Dyy, float[] Dyz, float[] Dzz,
                   float alpha, float beta, float C, boolean blackwhite,
                   float[] J, byte[] Vx, byte[] Vy, byte[] Vz, boolean enforce) {
        this.n0 = n0;
        this.n1 = n1;

        this.Dxx = Dxx;
        this.Dxy = Dxy;
        this.Dxz = Dxz;
        this.Dyy = Dyy;
        this.Dyz = Dyz;
        this.Dzz = Dzz;

        this.alpha = alpha;
        this.beta = beta;
        this.C = C;
//        this.BetaOne = beta_one;
//        this.BetaTwo = beta_two;
        this.blackwhite = blackwhite;

        this.J = J;
        this.Vx = Vx;
        this.Vy = Vy;
        this.Vz = Vz;

        this.enforce = enforce;

    }

    /**
     * Run the class instance thread.
     */
    public void run() {

        double[][] Ma = new double[3][3];
        double[][] Davec = new double[3][3];
        double[]   Daeig = new double[3];

        double Lambda1, Lambda2, Lambda3;
        double LambdaAbs1, LambdaAbs2, LambdaAbs3;
        double Ra, Rb, S, expRa, expRb, expS;
        double Voxel_data;
        int val;

        for (int i = n0; i < n1; ++i) {

            Ma[0][0]=Dxx[i]; Ma[0][1]=Dxy[i]; Ma[0][2]=Dxz[i];
            Ma[1][0]=Dxy[i]; Ma[1][1]=Dyy[i]; Ma[1][2]=Dyz[i];
            Ma[2][0]=Dxz[i]; Ma[2][1]=Dyz[i]; Ma[2][2]=Dzz[i];

            eigen_decomposition(Ma, Davec, Daeig);

            Lambda1 = Daeig[0];
            Lambda2 = Daeig[1];
            Lambda3 = Daeig[2];

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

            Voxel_data = Double.isNaN(Voxel_data)? 0 : Voxel_data ; // remove NaN

            // add result of this scale to output
            if (enforce) {
                J[i] = (float) Voxel_data;

                if (Vx!=null && Vy!=null && Vz!=null) {

                    // convert [-1,+1] --> byte [-128,+127] (java's byte data type)

                    Vx[i] = Frangi.to_byte(Davec[0][0]); // Davec[0][0]     --> Vx[i]
                    Vy[i] = Frangi.to_byte(Davec[1][0]); // Davec[1][0]     --> Vy[i]
                    Vz[i] = Frangi.to_byte(Davec[2][0]); // Davec[2][0]     --> Vz[i]

                }

//                    val = (int) Math.round(  ((Davec[0][0]+1)/2) * 255 ); // vx
//                    val = (val<0)?0:(val>255)?255:val;
//                    Vx[i] = (byte) val; // (unsigned char) val;

//                    val = (int) Math.round(  ((Davec[1][0]+1)/2) * 255 ); // vy
//                    val = (val<0)?0:(val>255)?255:val;
//                    Vy[i] = (byte) val; // (unsigned char) val;

//                    val = (int) Math.round(  ((Davec[2][0]+1)/2) * 255 ); // vz
//                    val = (val<0)?0:(val>255)?255:val;
//                    Vz[i] = (byte) val; // (unsigned char) val;

            }
            else {
                if (Voxel_data>J[i]) { // add if higher than current value
                    J[i] = (float) Voxel_data; // keep maximum filter response

                    if (Vx!=null && Vy!=null && Vz!=null) {

                        // convert [-1,+1] --> byte [-128,+127] (java's byte data type)

                        Vx[i] = Frangi.to_byte(Davec[0][0]); // Davec[0][0]     --> Vx[i]
                        Vy[i] = Frangi.to_byte(Davec[1][0]); // Davec[1][0]     --> Vy[i]
                        Vz[i] = Frangi.to_byte(Davec[2][0]); // Davec[2][0]     --> Vz[i]

                    }

//                        val = (int) Math.round(  ((Davec[0][0]+1)/2) * 255 ); // vx
//                        val = (val<0)?0:(val>255)?255:val;
//                        Vx[i] = (byte) val; // (unsigned char) val;

//                        val = (int) Math.round(  ((Davec[1][0]+1)/2) * 255 ); // vy
//                        val = (val<0)?0:(val>255)?255:val;
//                        Vy[i] = (byte) val; // (unsigned char) val;

//                        val = (int) Math.round(  ((Davec[2][0]+1)/2) * 255 ); // vz
//                        val = (val<0)?0:(val>255)?255:val;
//                        Vz[i] = (byte) val; // (unsigned char) val;

                }
            }
        }
    }

//    public static byte[] get_vesselness() { // is typically run after the threading
//        // go through the computed J[] values and find Jmin, Jmax,
//        float Jmin = Float.POSITIVE_INFINITY;
//        float Jmax = Float.NEGATIVE_INFINITY;
//
//        for (int i = 0; i < J.length; i++) {
//            if (J[i]<Jmin) Jmin = J[i];
//            if (J[i]>Jmax) Jmax = J[i];
//        }
//
//        byte[] J8 = new byte[J.length];
//
//        // wrap them in the [0,255] and convert to byte
//        if (Math.abs(Jmax-Jmin)<=Float.MIN_VALUE) {
//            for (int i = 0; i < J.length; ++i) {
//                J8[i] = (byte) 0;
//            }
//        }
//        else {
//            for (int i = 0; i < J.length; ++i) {
//                int Jnorm = Math.round(((J[i]-Jmin)/(Jmax-Jmin)) * 255);
//                Jnorm = (Jnorm<0)? 0 : (Jnorm > 255)? 255 : Jnorm; // wrap into [0,255] range
//                J8[i] = (byte) Jnorm; // convert min-max normalized J to J8 (byte8 image)
//            }
//        }
//
//        return J8;
//    }

    private void eigen_decomposition(double[][] A, double[][] V, double[] d) {

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
    private void tred2(double[][] V, double[] d, double[] e) {

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
    private void tql2(double[][] V, double[] d, double[] e) {

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

    private double hypot2(double x, double y) { return Math.sqrt(x*x+y*y); }

    private double absd(double val){ if(val>0){ return val;} else { return -val;} }

//    private int clamp(int x, int x1, int x2) {
//        int xC = (x<x1)?x1:x;
//        return (xC>x2)?x2:xC;
//    }

    private double MAX(double a, double b) {
        return (a>b)?a:b;
    }

}
