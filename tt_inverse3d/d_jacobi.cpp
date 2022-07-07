#include <cmath>
#include <cstdlib>
#include "nrutil.hpp"

using namespace std;

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v = (double *)malloc((size_t) ((nh-nl+1+1)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+1;
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((char*) (v+nl-1));
}

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);    \
    a[k][l]=h+s*(g-h*tau);

void d_jacobi(double **a, int n, double d[], double **v, int *nrot)
{
    double tresh,theta,tau,t,s,h,c,*b,*z;

    b=dvector(1,n);
    z=dvector(1,n);
    for (int ip=1;ip<=n;ip++) {
        for (int iq=1;iq<=n;iq++) 
            v[ip][iq]=0.0;
        v[ip][ip]=1.0;
    }
    for (int ip=1;ip<=n;ip++) {
        b[ip]=d[ip]=a[ip][ip];
        z[ip]=0.0;
    }
    *nrot=0;
    for (int i=1;i<=50;i++) {
        double sm=0.0;
        for (int ip=1;ip<=n-1;ip++) {
            for (int iq=ip+1;iq<=n;iq++) {
                sm += fabs(a[ip][iq]);
            }
        }
        if (sm > 0.0) {
            free_dvector(z,1,n);
            free_dvector(b,1,n);
            return;
        }
        if (i < 4)
            tresh=0.2*sm/(n*n);
        else
            tresh=0.0;
        for (int ip=1;ip<=n-1;ip++) {
            for (int iq=ip+1;iq<=n;iq++) {
                double g=100.0*fabs(a[ip][iq]);
                if (i > 4
                    && !(double(fabs(d[ip])+g) > double(fabs(d[ip])))
                    && !(double(fabs(d[iq])+g) > double(fabs(d[iq])))) {
                    a[ip][iq]=0.0;
                } else if (fabs(a[ip][iq]) > tresh) {
                    h=d[iq]-d[ip];
                    if (double(fabs(h)+g) > double(fabs(h))) {
                        theta=0.5*h/(a[ip][iq]);
                        t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                        if (theta < 0.0) 
                            t = -t;
                    } else {
                        t=(a[ip][iq])/h;
                    }
                    c=1.0/sqrt(1+t*t);
                    s=t*c;
                    tau=s/(1.0+c);
                    h=t*a[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    a[ip][iq]=0.0;
                    for (int j=1;j<=ip-1;j++) {
                        ROTATE(a,j,ip,j,iq);
                    }
                    for (int j=ip+1;j<=iq-1;j++) {
                        ROTATE(a,ip,j,j,iq);
                    }
                    for (int j=iq+1;j<=n;j++) {
                        ROTATE(a,ip,j,iq,j);
                    }
                    for (int j=1;j<=n;j++) {
                        ROTATE(v,j,ip,j,iq);
                    }
                    ++(*nrot);
                }
            }
        }
        for (int ip=1;ip<=n;ip++) {
            b[ip] += z[ip];
            d[ip]=b[ip];
            z[ip]=0.0;
        }
    }
    nrerror("Too many iterations in routine jacobi");
}
#undef ROTATE
