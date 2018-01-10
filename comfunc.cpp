#include "comfunc.h"
#include <math.h>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <stdarg.h>
#include <iostream>
#include <algorithm> 
extern double GetAveStdRMS(const double *a, int n, int opt)
{
	int i;
	double avg = 0.0, std = 0.0, rms = 0.0, sum = 0.0;

	if (n == 0) return 99999.9;

	for (i = 0; i<n; i++) {
		sum += a[i];
		rms += a[i] * a[i];
	}
	avg = sum / n;

	if (opt == 0) return avg;

	sum = 0.0;
	for (i = 0; i<n; i++) {
		sum += (a[i] - avg)*(a[i] - avg);
	}

	std = sqrt(sum / double(n - 1));
	rms = sqrt(rms / double(n));

	if (opt == 1) return std;
	else if (opt == 2) return rms;

	return 0.0;
}
#pragma region Matrixs
void Maddn(double *a, double *b, double *c, int m, int n)
{
	for (int i = 0; i<m; i++)
		for (int j = 0; j<n; j++)
		{
			*(c + i * n + j) = *(a + i * n + j) + *(b + i * n + j);
		}
}

void Madd(double *a, double *b, int m, int n)
{
	for (int i = 0; i<m; i++)
		for (int j = 0; j<n; j++)
		{
			*(a + i * n + j) = *(a + i * n + j) + *(b + i * n + j);
		}
}

void Mminn(double *a, double *b, double *c, int m, int n)
{
	int i, j;

	for (i = 0; i<m; i++)
		for (j = 0; j<n; j++)
		{
			*(c + i * n + j) = *(a + i * n + j) - *(b + i * n + j);
		}
}

void Mmin(double *a, double *b, int m, int n)
{
	int i, j;

	for (i = 0; i<m; i++)
		for (j = 0; j<n; j++)
		{
			*(a + i * n + j) = *(a + i * n + j) - *(b + i * n + j);
		}
}

void Mmulnm(double *a, double *b, int m, int n, int k, double *c)
{
	int i, j, l, u;
	for (i = 0; i <= m - 1; i++)
		for (j = 0; j <= k - 1; j++)
		{
			u = i * k + j;
			c[u] = 0.0;
			for (l = 0; l <= n - 1; l++)
				c[u] = c[u] + a[i*n + l] * b[l*k + j];
		}
}

void Mmul(double *a, int m, int n, double b)
{
	for (int i = 0; i<m; i++)
	{
		for (int j = 0; j<n; j++)
		{
			*(a + i * n + j) = *(a + i * n + j)*b;
		}
	}
}

void Mmuln(double *a, int m, int n, double b, double *c)
{
	for (int i = 0; i<m; i++)
	{
		for (int j = 0; j<n; j++)
		{
			*(c + i * n + j) = *(a + i * n + j)*b;
		}
	}
}

void Mtn(double *a, int m, int n, double *b)
{
	for (int l = 0; l<n; l++)
		for (int k = 0; k<m; k++)
		{
			b[l*m + k] = a[k*n + l];
		}

}

void Mt(double *a, int m, int n)
{
	double *at;
	at = (double *)malloc(n*m * sizeof(double));
	for (int i = 0; i<m; i++)
		for (int j = 0; j<n; j++)
		{
			*(at + i * n + j) = *(a + i * n + j);
		}
	for (int l = 0; l<n; l++)
		for (int k = 0; k<m; k++)
		{
			a[l*m + k] = at[k*n + l];
		}
	free(at);
}

double Minv(double a[], int n)
{
	int *is, *js, i, j, k, l, u, v;
	double d, p;
	is = (int *)malloc(n * sizeof(int));
	js = (int *)malloc(n * sizeof(int));
	for (k = 0; k <= n - 1; k++)
	{
		d = 0.0;
		for (i = k; i <= n - 1; i++)
			for (j = k; j <= n - 1; j++)
			{
				l = i * n + j;
				p = fabs(a[l]);
				if (p>d)
				{
					d = p;
					is[k] = i;
					js[k] = j;
				}
			}
		if (d == 0)//d+1.0==1.0 lcc
		{
			free(is);
			free(js);
			printf("\nerror:inverse matrix is not exist\n");
			return (0);
		}
		if (is[k] != k)
			for (j = 0; j <= n - 1; j++)
			{
				u = k * n + j; v = is[k] * n + j;
				p = a[u]; a[u] = a[v]; a[v] = p;
			}
		if (js[k] != k)
			for (i = 0; i <= n - 1; i++)
			{
				u = i * n + k;
				v = i * n + js[k];
				p = a[u];
				a[u] = a[v]; a[v] = p;
			}
		l = k * n + k;
		a[l] = 1.0 / a[l];
		for (j = 0; j <= n - 1; j++)
			if (j != k)
			{
				u = k * n + j;
				a[u] = a[u] * a[l];
			}
		for (i = 0; i <= n - 1; i++)
			if (i != k)
				for (j = 0; j <= n - 1; j++)
					if (j != k)
					{
						u = i * n + j;
						a[u] = a[u] - a[i*n + k] * a[k*n + j];
					}
		for (i = 0; i <= n - 1; i++)
			if (i != k)
			{
				u = i * n + k; a[u] = -a[u] * a[l];
			}
	}
	for (k = n - 1; k >= 0; k--)
	{
		if (js[k] != k)
			for (j = 0; j <= n - 1; j++)
			{
				u = k * n + j; v = js[k] * n + j;
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}
		if (is[k] != k)
			for (i = 0; i <= n - 1; i++)
			{
				u = i * n + k;
				v = i * n + is[k];
				p = a[u];
				a[u] = a[v]; a[v] = p;
			}
	}

	free(is);
	free(js);
	return(1);
}

double Minvn(double a[], int n, double *b)
{
	for (int i = 0; i<n; i++)
		for (int j = 0; j<n; j++)
		{
			*(b + i * n + j) = *(a + i * n + j);
		}

	int *is, *js, i, j, k, l, u, v;
	double d, p;
	is = (int *)malloc(n * sizeof(int));
	js = (int *)malloc(n * sizeof(int));
	for (k = 0; k <= n - 1; k++)
	{
		d = 0.0;
		for (i = k; i <= n - 1; i++)
			for (j = k; j <= n - 1; j++)
			{
				l = i * n + j;
				p = fabs(b[l]);
				if (p>d)
				{
					d = p;
					is[k] = i;
					js[k] = j;
				}
			}
		if (d == 0)//d+1.0==1.0 lcc
		{
			free(is);
			free(js);
			printf("\nerror:inverse matrix is not exist\n");
			return (0);
		}
		if (is[k] != k)
			for (j = 0; j <= n - 1; j++)
			{
				u = k * n + j; v = is[k] * n + j;
				p = b[u]; b[u] = b[v]; b[v] = p;
			}
		if (js[k] != k)
			for (i = 0; i <= n - 1; i++)
			{
				u = i * n + k;
				v = i * n + js[k];
				p = b[u];
				b[u] = b[v]; b[v] = p;
			}
		l = k * n + k;
		b[l] = 1.0 / b[l];
		for (j = 0; j <= n - 1; j++)
			if (j != k)
			{
				u = k * n + j;
				b[u] = b[u] * b[l];
			}
		for (i = 0; i <= n - 1; i++)
			if (i != k)
				for (j = 0; j <= n - 1; j++)
					if (j != k)
					{
						u = i * n + j;
						b[u] = b[u] - b[i*n + k] * b[k*n + j];
					}
		for (i = 0; i <= n - 1; i++)
			if (i != k)
			{
				u = i * n + k; b[u] = -b[u] * b[l];
			}
	}
	for (k = n - 1; k >= 0; k--)
	{
		if (js[k] != k)
			for (j = 0; j <= n - 1; j++)
			{
				u = k * n + j; v = js[k] * n + j;
				p = b[u];
				b[u] = b[v];
				b[v] = p;
			}
		if (is[k] != k)
			for (i = 0; i <= n - 1; i++)
			{
				u = i * n + k;
				v = i * n + is[k];
				p = b[u];
				b[u] = b[v]; b[v] = p;
			}
	}

	free(is);
	free(js);
	return(1);
}

double Mrem(double *a, int i, int j, int n)
{
	int k, m;
	double  *pTemp = new double[(n - 1)*(n - 1)];

	for (k = 0; k<i; k++)
	{
		for (m = 0; m<j; m++)
		{
			pTemp[k*(n - 1) + m] = a[k*n + m];
		}
		for (m = j; m<n - 1; m++)
		{
			pTemp[k*(n - 1) + m] = a[k*n + m + 1];
		}
	}
	for (k = i; k<n - 1; k++)
	{
		for (m = 0; m<j; m++)
		{
			pTemp[k*(n - 1) + m] = a[(k + 1)*n + m];
		}
		for (m = j; m<n - 1; m++)
		{
			pTemp[k*(n - 1) + m] = a[(k + 1)*n + m + 1];
		}
	}
	double  dResult = (((i + j) % 2 == 1) ? -1 : 1)*Mdet(pTemp, n - 1);
	delete[] pTemp;
	return  dResult;
}

double Mdet(double *a, int n)   //?¨®DD¨¢D¨º?¦Ì??¦Ì    
{
	if (n == 1)
	{
		return a[0];
	}
	double sum = 0;
	for (int j = 0; j<n; j++)
	{
		sum += a[0 * n + j] * Mrem(a, 0, j, n);
	}
	return sum;
}

void Mequalm(double *M, int m, int n, double *N)
{
	for (int i = 0; i<m; i++)
	{
		for (int j = 0; j<n; j++)
		{
			*(N + i * n + j) = *(M + i * n + j);
		}
	}
}

void Mequal(double *M, int m, int n, double a)
{
	for (int i = 0; i<m; i++)
	{
		for (int j = 0; j<n; j++)
		{
			*(M + i * n + j) = a;
		}
	}
}

double Mmean(double *a, int m)
{
	double sum = 0;
	for (int i = 0; i<m; i++)
	{
		sum = sum + a[i];
	}
	return sum / m;
}

int eejcb(double a[], int n, double v[], double eps, int jt)
{
	int i, j, p, q, u, w, t, s, l;
	double fm, cn, sn, omega, x, y, d;
	l = 1;
	for (i = 0; i <= n - 1; i++)
	{
		v[i*n + i] = 1.0;
		for (j = 0; j <= n - 1; j++)
		{
			if (i != j)
			{
				v[i*n + j] = 0.0;
			}
		}
	}
	while (1 == 1)
	{
		fm = 0.0;
		for (i = 0; i <= n - 1; i++)
		{
			for (j = 0; j <= n - 1; j++)
			{
				d = fabs(a[i*n + j]);
				if ((i != j) && (d>fm))
				{
					fm = d;
					p = i;
					q = j;
				}
			}
		}
		if (fm<eps)
		{
			return(1);
		}
		if (l>jt)
		{
			return(-1);
		}
		l = l + 1;
		u = p * n + q;
		w = p * n + p;
		t = q * n + p;
		s = q * n + q;
		x = -a[u];
		y = (a[s] - a[w]) / 2.0;
		omega = x / sqrt(x*x + y * y);
		if (y<0.0)
		{
			omega = -omega;
		}
		sn = 1.0 + sqrt(1.0 - omega * omega);
		sn = omega / sqrt(2.0*sn);
		cn = sqrt(1.0 - sn * sn);
		fm = a[w];
		a[w] = fm * cn*cn + a[s] * sn*sn + a[u] * omega;
		a[s] = fm * sn*sn + a[s] * cn*cn - a[u] * omega;
		a[u] = 0.0;
		a[t] = 0.0;
		for (j = 0; j <= n - 1; j++)
		{
			if ((j != p) && (j != q))
			{
				u = p * n + j;
				w = q * n + j;
				fm = a[u];
				a[u] = fm * cn + a[w] * sn;
				a[w] = -fm * sn + a[w] * cn;
			}
		}
		for (i = 0; i <= n - 1; i++)
		{
			if ((i != p) && (i != q))
			{
				u = i * n + p;
				w = i * n + q;
				fm = a[u];
				a[u] = fm * cn + a[w] * sn;
				a[w] = -fm * sn + a[w] * cn;
			}
		}
		for (i = 0; i <= n - 1; i++)
		{
			u = i * n + p;
			w = i * n + q;
			fm = v[u];
			v[u] = fm * cn + v[w] * sn;
			v[w] = -fm * sn + v[w] * cn;
		}
	}
	return(1);
}

void Munit(double* a, int n)
{
	for (int i = 0; i<n; i++)
	{
		for (int j = 0; j<n; j++)
		{
			a[i*n + j] = 0.0;
		}
	}
	for (int i = 0; i<n; i++)
	{
		a[i*n + i] = 1.0;
	}
}
#pragma endregion
