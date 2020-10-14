#include<stdio.h>
#include<math.h>
#include<stdlib.h>
/*
 standard O(n^3) gaussian elimination matrix solver
*/
void d_solvec(double **a,double *b,int *iflag,int n)
{
  int i,j,k,l,flag;
  double fact;
  double temp;
  double sum;
  double eps=1e-8;

  
  for(i=0;i<n;i++)
    {
      if (fabs(a[i][i]) < eps)
	{
	  flag=1;
	  for(k=i+1;k<n && flag;k++)
	    {
	      if (a[k][i]!=0)
                {
		  flag=0;
		  for(l=0;l<n;l++)
		    {
		      temp=a[k][l];
		      a[k][l]=a[i][l];
		      a[i][l]=temp;
		    }
		  temp=b[k];
		  b[k]=b[i];
		  b[i]=temp;
                }
	    }
	  if (flag) {*iflag=0;return;}
	}
      for(k=i+1;k<n;k++)
	{
	  if (i!=k)
	    {
	      fact=-a[k][i]/a[i][i];
	      for(j=0;j<n;j++)
		{
		  a[k][j]+=fact*a[i][j];
		}
	      b[k]+=fact*b[i];
	    }
	}
    }

  for(i=n-1;i>=0;i--)
    {
      sum=0;
      for(j=i+1;j<n;j++)
	sum+=a[i][j]*b[j];
      b[i]=(b[i]-sum)/a[i][i];
    }
  *iflag=1;
  return;

}
/*
Interpolate using RBF + linear polynomial 
d=number of dimensions
np=number of points
xcloud=cloud of points xcloud[3*np]
P=target point
weights=interpolation weights
*/
void d_interprbf(double *xcloud, double *P,double *weights,int np, int d)
{
  int i,j,n,imin;
  double **M, *B;
  double dd,dmin;
  int neqns=np+d+1;
  int iflag;
  
  
  M=(double **)malloc(sizeof(double *)*neqns);
  B=(double *)malloc(sizeof(double) *neqns);
  for(i=0;i<neqns;i++) M[i]=(double *)malloc(sizeof(double)*neqns);
  
  dmin=1.0e15;
  imin=0;

  for(i=0;i<np;i++)
    {
      for(j=0;j<np;j++)
	{
	  dd=0;
	  for(n=0;n<d;n++) dd+=(xcloud[d*i+n]-xcloud[d*j+n])*(xcloud[d*i+n]-xcloud[d*j+n]);
	  M[i][j]=exp(-dd);
	}
      M[i][np]=M[np][i]=1.0;
      for(j=np+1,n=0;j<np+d+1;j++,n++) M[i][j]=M[j][i]=xcloud[d*i+n];
      dd=0;
      for(n=0;n<d;n++) dd+=(xcloud[d*i+n]-P[n])*(xcloud[d*i+n]-P[n]);
      if (dd < dmin) {
	imin=i;
	dmin=dd;
      }
      B[i]=exp(-dd);
    }
  for(i=np;i<np+d+1;i++) 
    {
      for(j=np;j<np+d+1;j++)
	M[i][j]=0.0;
      B[i]=(i==np)?1.0:P[i-(np+1)];
    }
  d_solvec(M,B,&iflag,neqns);
  
  if (iflag==0) 
    {
      for(i=0;i<np;i++) weights[i]=0.0;
      weights[imin]=1.0;
    }
  else 
    {
      for(i=0;i<np;i++)weights[i]=B[i];
    }

  for(i=0;i<neqns;i++) free(M[i]);
  free(M);
  free(B);
}

/*
 * test
void main()
{
  double xcloud[24]={0,0,0,1,0,0,1,1,0,0,1,0,0,0,1,1,0,1,1,1,1,0,1,1};
  double P[3]={0.5,0.5,0.5};
  double xcloud2[8]={0,0,1,0,1,1,0,1};
  double Q[2]={0.5,0.5};

  double weights[8];
  int i;
  
  printf("P=%f %f %f\n",P[0],P[1],P[2]);
  d_interprbf(xcloud,P,weights,8,3);
  for(i=0;i<8;i++) printf("%f ",weights[i]);
  printf("\n");
  d_interprbf(xcloud2,Q,weights,4,2);
  for(i=0;i<4;i++) printf("%f ",weights[i]);
  printf("\n");
}
*/
