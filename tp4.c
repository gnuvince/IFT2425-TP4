#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define H     0.00005
#define START 0
#define END   150



static double e = 0;
static double a = 0;
static double w = 0;

/*
  Équation différentielle:

  y''(t) - y'(t)(e - y²(t)) + y(t) = a*cos(wt)

  <=>

  y''(t) = y'(t)(e - y²(t)) - y(t) + a*cos(wt)

  Avec:
    y(0)  = 0
    y'(0) = 0

  Posons: x(t) = y'(t):

    x'(t) = x(t)(e - y²(t)) - y(t) + a*cos(wt) = f1(t, y(t), x(t))
    y'(t) = x(t)                               = f2(t, y(t), x(t))
 */


double f1(double t, double yt, double xt) {
    return xt*(e - yt*yt) - yt + a*cos(w*t);
}

double f2(double t, double yt, double xt) {
    (void)t, (void)yt; /* Force evaluation of y and yt to silence compiler warnings. */
    return xt;
}



int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Utilisation: %s <a> <e> <w>\n", argv[0]);
        exit(1);
    }

    a = atof(argv[1]);
    e = atof(argv[2]);
    w = atof(argv[3]);

    double tn = START;
    double k[2][6];
    double yn = 0.0;
    double xn = 0.0;

    while (tn < END) {
        k[0][0] = H * f1(tn, yn, xn);
        k[1][0] = H * f2(tn, yn, xn);

        k[0][1] = H * f1(tn + H/4, yn + k[1][0]/4, xn + k[0][0]/4);
        k[1][1] = H * f2(tn + H/4, yn + k[1][0]/4, xn + k[0][0]/4);

        k[0][2] = H * f1(tn + H*3/8,
                         yn + k[1][0]*3/32 + k[1][1]*9/32,
                         xn + k[0][0]*3/32 + k[0][1]*9/32);
        k[1][2] = H * f2(tn + H*3/8,
                         yn + k[1][0]*3/32 + k[1][1]*9/32,
                         xn + k[0][0]*3/32 + k[0][1]*9/32);

        k[0][3] = H * f1(tn + H*12/13,
                         yn + k[1][0]*1932/2197 - k[1][1]*7200/2197 + k[1][2]*7296/7200,
                         xn + k[0][0]*1932/2197 - k[0][1]*7200/2197 + k[0][2]*7296/7200);
        k[1][3] = H * f2(tn + H*12/13,
                         yn + k[1][0]*1932/2197 - k[1][1]*7200/2197 + k[1][2]*7296/7200,
                         xn + k[0][0]*1932/2197 - k[0][1]*7200/2197 + k[0][2]*7296/7200);

        k[0][4] = H * f1(tn + H,
                         yn + k[1][0]*439/216 - k[1][1]*8 + k[1][2]*3680/513 - k[1][3]*845/4104,
                         xn + k[0][0]*439/216 - k[0][1]*8 + k[0][2]*3680/513 - k[0][3]*845/4104);
        k[1][4] = H * f2(tn + H,
                         yn + k[1][0]*439/216 - k[1][1]*8 + k[1][2]*3680/513 - k[1][3]*845/4104,
                         xn + k[0][0]*439/216 - k[0][1]*8 + k[0][2]*3680/513 - k[0][3]*845/4104);

        k[0][5] = H * f1(tn + H/2,
                         yn - k[1][0]*8/27 + k[1][1]*2 - k[1][2]*3544/2565 + k[1][3]*1859/4104 - k[1][4]*11/40,
                         xn - k[0][0]*8/27 + k[0][1]*2 - k[0][2]*3544/2565 + k[0][3]*1859/4104 - k[0][4]*11/40);
        k[1][5] = H * f2(tn + H/2,
                         yn - k[1][0]*8/27 + k[1][1]*2 - k[1][2]*3544/2565 + k[1][3]*1859/4104 - k[1][4]*11/40,
                         xn - k[0][0]*8/27 + k[0][1]*2 - k[0][2]*3544/2565 + k[0][3]*1859/4104 - k[0][4]*11/40);

        double yn1 = yn + k[1][0]*16/135 + k[1][2]*6656/12825 + k[1][3]*28561/56430 - k[1][4]*9/50 + k[1][5]*2/55;
        double xn1 = xn + k[0][0]*16/135 + k[0][2]*6656/12825 + k[0][3]*28561/56430 - k[0][4]*9/50 + k[0][5]*2/55;

        printf("%g %g\n", yn1, xn);

        yn = yn1;
        xn = xn1;
        tn += H;
    }

    return EXIT_SUCCESS;
}
