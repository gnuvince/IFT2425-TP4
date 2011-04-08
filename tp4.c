#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define H     0.00005
#define START 0
#define END   150


#define NAME_IMG_OUT "VincentEric.pgm"

#define WIDTH  512
#define LENGTH 512

#define MID_WDTH WIDTH/2
#define MID_LGTH LENGTH/2

#define WHITE 255
#define BLACK 0



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


double** fmatrix_allocate_2d(int vsize,int hsize)
{
    int i;
    double** matrix;
    double *imptr;

    matrix=(double**)malloc(sizeof(double*)*vsize);
    if (matrix==NULL) printf("probleme d'allocation memoire");

    imptr=(double*)malloc(sizeof(double)*hsize*vsize);
    if (imptr==NULL) printf("probleme d'allocation memoire");

    for(i=0;i<vsize;i++,imptr+=hsize) matrix[i]=imptr;
    return matrix;
}

//----------------------------------------------------------*/
/* Libere la memoire de la matrice 2d de double              */
/*----------------------------------------------------------*/
void free_fmatrix_2d(double** pmat)
{
    free(pmat[0]);
    free(pmat);
}

/*----------------------------------------------------------*/
/* Sauvegarde de l'image de nom <name> au format pgm        */
/*----------------------------------------------------------*/
void SaveImagePgm(char* name,double** mat,int length,int width)
{
    int i,j;
    FILE* fic;


    /*--ouverture fichier--*/
    fic=fopen(NAME_IMG_OUT, "w");
    if (fic==NULL)
    { printf(" Probleme dans la sauvegarde de %s",NAME_IMG_OUT);
        exit(-1); }
    printf("\n Sauvegarde de %s au format pgm\n",name);

    /*--sauvegarde de l'entete--*/
    fprintf(fic,"P5");
    fprintf(fic,"\n# IMG Module");
    fprintf(fic,"\n %d %d",length,width);
    fprintf(fic,"\n255\n");

    /*--enregistrement--*/
    for(i=0;i<length;i++)
        for(j=0;j<width;j++)
            fprintf(fic,"%c",(char)mat[i][j]);

    /*--fermeture fichier--*/
    fclose(fic);
}

/*-------------------------------------------------------------*/
/*--- SetPointBlack -------------------------------------------*/
/*-------------------------------------------------------------*/
void SetPointBlack(double** Mat,int lgth,int wdth,int row,int col)
{
    Mat[lgth-row][col+wdth]=0.0;
}



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
    double **matrice = fmatrix_allocate_2d(WIDTH, LENGTH);

    for (int i = 0; i < WIDTH; ++i)
        for (int j = 0; j < WIDTH; ++j)
            matrice[i][j] = WHITE;


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

        SetPointBlack(matrice, MID_LGTH, MID_WDTH, (int)(xn1*MID_WDTH) % MID_WDTH, (int)(yn1*MID_LGTH) % MID_LGTH);

        yn = yn1;
        xn = xn1;
        tn += H;
    }

    SaveImagePgm(NAME_IMG_OUT, matrice, LENGTH, WIDTH);
    free_fmatrix_2d(matrice);

    return EXIT_SUCCESS;
}
