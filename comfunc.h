#pragma once
extern double GetAveStdRMS(const double *a, int n, int opt);
/* Matrixs calculation function,
2018.1.3 
author:ShiLei   ---------- */
#pragma region Matrixs
//c[m][n]=a[m][n]+b[m][n]
void Maddn(double *a, double *b, double *c, int m, int n);
//a[m][n]=a[m][n]+b[m][n]
void Madd(double *a, double *b, int m, int n);
//c[m][n]=a[m][n]-b[m][n]
void Mminn(double *a, double *b, double *c, int m, int n);
//a[m][n]=a[m][n]-b[m][n]
void Mmin(double *a, double *b, int m, int n);
//c[m][k]=a[m][n]*b[n][k]
void Mmulnm(double *a, double *b, int m, int n, int k, double *c);
//a[m][n]=a[m][k]*b
void Mmul(double *a, int m, int n, double b);
//c[m][n]=a[m][k]*b
void Mmuln(double *a, int m, int n, double b, double *c);
//b=aT
void Mtn(double *a, int m, int n, double *b);
//a=aT
void Mt(double *a, int m, int n);
//inv(a)
double Minv(double a[], int n);
//b=inv(a)
double Minvn(double a[], int n, double *b);
//A* adjoint matrix  
double Mrem(double *a, int i, int j, int n);
//det 
double Mdet(double *a, int n);
//N[m][n]=M[m][n]
void Mequalm(double *M, int m, int n, double *N);
//M[m][n]=a
void Mequal(double *M, int m, int n, double a);
//mean of col
double Mmean(double *a, int m);
//?¨®¨º¦Ì??3????¨®¦Ì?¨¬??¡Â?¦Ì?¡ã¨¬??¡Â?¨°¨¢?¦Ì?????¡À¨¨¡¤? 
//¨¤?¨®?????¡À¨¨(Jacobi)¡¤?¡¤??¨®¨º¦Ì??3????¨®¦Ì?¨¨?2?¨¬??¡Â?¦Ì?¡ã¨¬??¡Â?¨°¨¢? 
//¡¤¦Ì???¦ÌD?¨®¨²0¡À¨ª¨º?3?1y¦Ì¨¹?¨²jt??¨¨?????¦Ì????¨¨¨°a?¨® 
//¡¤¦Ì???¦Ì?¨®¨®¨²0¡À¨ª¨º??y3¡ê¡¤¦Ì?? 
//a-3€?¨¨?an*n¦Ì?¨ºy¡Á¨¦¡ê???¡¤?¨º¦Ì??3????¨®¡ê?¡¤¦Ì??¨º¡À????????¡¤?n??¨¬??¡Â?¦Ì 
//n-???¨®¦Ì??¡Á¨ºy 
//u-3€?¨¨?an*n¦Ì?¨ºy¡Á¨¦¡ê?¡¤¦Ì??¨¬??¡Â?¨°¨¢?(¡ã?¨¢D???¡é) 
//eps-???????¨¨¨°a?¨® 
//jt-??D¨ª¡À?¨¢?¡ê?????¡Á??¨®¦Ì¨¹?¨²??¨ºy 
int Meejcb(double a[], int n, double v[], double eps, int jt);
//eye
void Munit(double* a, int n);
#pragma endregion
/* end --------------------------------------------------------- */
