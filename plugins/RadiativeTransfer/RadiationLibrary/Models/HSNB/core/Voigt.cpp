#include "Voigt.h"
#include <complex>

//==============================================================================

double Voigt::humlicek(double bl, double bd, double ds)
{
    double rtln2 = 0.8325546;
    double rtpi  = 1.772454;

    // Avoid divide by zero
    if (bd < 1.0e-10) bd = 1.0e-10;
    double bd2 = bd / rtln2;

    double x = ds/bd2;
    double y = bl/bd2;
    std::complex<double> t(y, -x);
    double s = std::abs(x) + y;
    std::complex<double> w4;
    if (s >= 15.0) {
        // Region I
        w4 = t*0.5641896/(0.5 + t*t);
        return w4.real()/(rtpi*bd2);
    }

    if (s >= 5.5) {
        // Region II
        std::complex<double> u = t*t;
        w4 = t*(1.410474 + u*0.5641896)/(0.75 + u*(3.0+u));
        return w4.real()/(rtpi*bd2);
    }

    if (y >= 0.195*std::abs(x) - 0.176) {
        // Region III
        w4 = (16.4955+t*(20.20933+t*(11.96482+t*(3.778987+t*.5642236))))/
            (16.4955+t*(38.82363+t*(39.27121+t*(21.69274+t*(6.699398+t)))));
        return w4.real()/(rtpi*bd2);
    }

    // Region IV
    std::complex<double> u = t*t;
    w4 = (std::exp(u)-t*(36183.31-u*(3321.9905-u*(1540.787-u*(219.0313
        -u*(35.76683-u*(1.320522-u*.56419))))))
        /(32066.6-u*(24322.84-u*(9022.228-u*(2186.181-u*(364.2191
        -u*(61.57037-u*(1.841439-u))))))));
    return w4.real()/(rtpi*bd2);
}

//==============================================================================

double drayson_helper(const double& x, const double& y)
{
    static double b[22], ri[15], hn[25];
    static double c0, d0[25], d1[25], d2[25], d3[25], d4[25];

    const double c[21] = {
         0.7093602e-7, -0.2518434e-6,  0.8566874e-6,
        -0.2787638e-5,  0.8660740e-5, -0.2565551e-4,
         0.7228775e-4, -0.1933631e-3,  0.4899520e-3,
        -0.1173267e-2,  0.2648762e-2, -0.5623190e-2,
         0.1119601e-1, -0.2084976e-1,  0.3621573e-1,
        -0.5851412e-1,  0.8770816e-1, -0.1216640,
         0.15584     , -0.184       ,  0.2
    };

    const double h = 0.201;

    const double xn[15] = {
        10.0, 9.0, 8.0, 8.0, 7.0, 6.0, 5.0, 4.0,
         3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0
    };
    const double yn[15] = {
        0.6, 0.6, 0.6, 0.5, 0.4, 0.4, 0.3, 0.3,
        0.3, 0.3, 1.0, 0.9, 0.8, 0.7, 0.7
    };
    const double hh[3] = { 0.2562121, 0.2588268e-1, 0.2820948 };
    const double xx[3] = { 0.5246476, 1.65068, 0.7071068 };
    const double nby2[19] = {
        9.5, 9.0, 8.5, 8.0, 7.5, 7.0, 6.5, 6.0, 5.5, 5.0,
        4.5, 4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.0, 0.5
    };

    double u, uu, v, vv, y2, dx;
    int i, j, n, min, max;

    // Setup data on the first call
    static bool first = true;
    if (first) {
        first = false;

        for (i = 0; i < 15; ++i)
            ri[i] = -(i+1)/2.0;

        b[0] = 0.0;
        b[1] = 0.7093602e-7;

        for (i = 0; i < 25; ++i) {
            hn[i] = h*(i+0.5);
            c0 = 4.0*hn[i]*hn[i]/25.0 - 2.0;
            for (j = 1; j < 21; ++j)
                b[j+1] = c0*b[j] - b[j-1] + c[j];
            d0[i] = hn[i]*(b[21]-b[20])/5.0;
            d1[i] = 1.0-2.0*hn[i]*d0[i];
            d2[i] = (hn[i]*d1[i] + d0[i])/ri[1];
            d3[i] = (hn[i]*d2[i] + d1[i])/ri[2];
            d4[i] = (hn[i]*d3[i] + d2[i])/ri[3];
        }
    }

    // REGION II
    if ((x < 5.0) && (y > 1.0) && x <= (1.85*(3.6-y))) {
        //r2++;
        if (y >= 1.45)
            i = (int) (y + y) - 1;
        else
            i = (int) (11.0*y) - 1;

        j = (int) (x + x + 1.85) - 1;

        max = (int) (xn[j]*yn[i] + 0.46);
        min = std::min(16, 21-2*max) - 1;

        uu = y;
        vv = x;

        for (j = min; j < 19; ++j) {
            u = nby2[j] / (uu*uu + vv*vv);
            uu = y + u*uu;
            vv = x - u*vv;
        }

        return (uu / (uu*uu + vv*vv) / 1.772454);
    }

    // REGION I
    if ((x < 5.0) && (y <= 1.0) && x+y < 5.0) {
        //r1++;
        y2 = y*y;
        n = (int)(x/h);
        dx = x - hn[n];
        u = (((d4[n]*dx+d3[n])*dx+d2[n])*dx+d1[n])*dx+d0[n];
        v = 1.0-2.0*x*u;

        // Taylor series expansion
        vv = std::exp(y2-x*x) * std::cos(2.0*x*y)/1.128379 - y*v;
        uu = -y;
        max = (int) (5.0 + (12.5-x)*0.8*y);

        for (i = 1; i < max; i += 2) {
            u = (x*v+u) / ri[i];
            v = (x*u+v) / ri[i+1];
            uu = -uu*y2;
            vv = vv + v*uu;
        }

        return (1.128379*vv);
    }

    // REGION IIIB
    if ( ((x < 5.0) && (y > 1.0) && (x >= 1.85*(3.6-y))) ||
         ((x >= 5.0) && (y >= (11.0-0.6875*x))) ) {
        //r3b++;
        y2 = y*y;
        u = x - xx[2];
        v = x + xx[2];
        return (y * (hh[2]/(y2+u*u) + hh[2]/(y2+v*v)));
    }

    // REGION IIIA
    //r3a++;
    y2 = y*y;
    u = x - xx[0];
    v = x + xx[0];
    uu = x - xx[1];
    vv = x + xx[1];
    return (y*(hh[0]/(y2+u*u)+hh[0]/(y2+v*v)+hh[1]/(y2+uu*uu)+hh[1]/(y2+vv*vv)));
}

//==============================================================================

double Voigt::drayson(double bl, double bd, double ds)
{
    double rtln2 = 0.8325546;
    double rtpi  = 1.772454;

    // Avoid divide by zero
    if (bd < 1.0e-10) bd = 1.0e-10;
    double bd2 = bd / rtln2;

    double x = std::abs(ds)/bd2;
    double y = bl/bd2;

    return drayson_helper(x, y)/(rtpi*bd2);
}

double new_drayson(const double& x, const double& y)
{
    static double b[22], ri[15], hn[25];
    static double c0, d0[25], d1[25], d2[25], d3[25], d4[25];

    const double c[21] = {
         0.7093602e-7, -0.2518434e-6,  0.8566874e-6,
        -0.2787638e-5,  0.8660740e-5, -0.2565551e-4,
         0.7228775e-4, -0.1933631e-3,  0.4899520e-3,
        -0.1173267e-2,  0.2648762e-2, -0.5623190e-2,
         0.1119601e-1, -0.2084976e-1,  0.3621573e-1,
        -0.5851412e-1,  0.8770816e-1, -0.1216640,
         0.15584     , -0.184       ,  0.2
    };

    const double h = 0.201;

    const double xn[15] = {
        10.0, 9.0, 8.0, 8.0, 7.0, 6.0, 5.0, 4.0,
         3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0
    };
    const double yn[15] = {
        0.6, 0.6, 0.6, 0.5, 0.4, 0.4, 0.3, 0.3,
        0.3, 0.3, 1.0, 0.9, 0.8, 0.7, 0.7
    };
    const double hh[3] = { 0.2562121, 0.2588268e-1, 0.2820948 };
    const double xx[3] = { 0.5246476, 1.65068, 0.7071068 };
    const double nby2[19] = {
        9.5, 9.0, 8.5, 8.0, 7.5, 7.0, 6.5, 6.0, 5.5, 5.0,
        4.5, 4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.0, 0.5
    };

    double u, uu, v, vv, y2, dx;
    int i, j, n, min, max;

    // Setup data on the first call
    static bool first = true;
    if (first) {
        first = false;

        for (i = 0; i < 15; ++i)
            ri[i] = -(i+1)/2.0;

        b[0] = 0.0;
        b[1] = 0.7093602e-7;

        for (i = 0; i < 25; ++i) {
            hn[i] = h*(i+0.5);
            c0 = 4.0*hn[i]*hn[i]/25.0 - 2.0;
            for (j = 1; j < 21; ++j)
                b[j+1] = c0*b[j] - b[j-1] + c[j];
            d0[i] = hn[i]*(b[21]-b[20])/5.0;
            d1[i] = 1.0-2.0*hn[i]*d0[i];
            d2[i] = (hn[i]*d1[i] + d0[i])/ri[1];
            d3[i] = (hn[i]*d2[i] + d1[i])/ri[2];
            d4[i] = (hn[i]*d3[i] + d2[i])/ri[3];
        }
    }

    // REGION IIIB
    if (y >= 11.0 - 0.6875*x) {
        //r3b++;
        y2 = y*y;
        u = x - xx[2];
        v = x + xx[2];
        return (y * (hh[2]/(y2+u*u) + hh[2]/(y2+v*v)));
    // REGION IIIA
    } else if (y >= 3.6 - 0.72*x) {
        //r3a++;
        y2 = y*y;
        u = x - xx[0];
        v = x + xx[0];
        uu = x - xx[1];
        vv = x + xx[1];
        return (y*(hh[0]/(y2+u*u)+hh[0]/(y2+v*v)+hh[1]/(y2+uu*uu)+hh[1]/(y2+vv*vv)));
    // REGION II
    } else if (y >= 1.0) {
        //r2++;
        if (y >= 1.45)
            i = (int) (y + y) - 1;
        else
            i = (int) (11.0*y) - 1;

        j = (int) (x + x + 1.85) - 1;

        max = (int) (xn[j]*yn[i] + 0.46);
        min = std::min(16, 21-2*max) - 1;

        uu = y;
        vv = x;

        for (j = min; j < 19; ++j) {
            u = nby2[j] / (uu*uu + vv*vv);
            uu = y + u*uu;
            vv = x - u*vv;
        }

        return (uu / (uu*uu + vv*vv) / 1.772454);
    // REGION I
    } else {
        //r1++;
        y2 = y*y;
        n = (int)(x/h);
        dx = x - hn[n];
        u = (((d4[n]*dx+d3[n])*dx+d2[n])*dx+d1[n])*dx+d0[n];
        v = 1.0-2.0*x*u;

        // Taylor series expansion
        vv = std::exp(y2-x*x) * std::cos(2.0*x*y)/1.128379 - y*v;
        uu = -y;
        max = (int) (5.0 + (12.5-x)*0.8*y);

        for (i = 1; i < max; i += 2) {
            u = (x*v+u) / ri[i];
            v = (x*u+v) / ri[i+1];
            uu = -uu*y2;
            vv = vv + v*uu;
        }

        return (1.128379*vv);
    }
}



