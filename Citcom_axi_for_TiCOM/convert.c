/* read in velocity near the plume at processor 12, i=27, j=7 for 57h.5530 */

#include <math.h>
#include <stdio.h>
#include <fcntl.h>

#define PI	3.141593
#define PI2	PI * 2

int n1,n2,n3,n4;

main(argc,argv)
     int argc;
     char **argv;
{

  FILE *fp1,*fp2,*fp3,*fp4,*fp5;
  char input_s[200],filename1[250],filename2[250],inputf[250],outputf[250];
  int stride,nno,nx,step,i,j,i1,j1,k,jj,nz,nz1,n,m,frame,proc;
  float temp1,temp2,temp3,temp4,temp5,temp6,tmax,*t1,*t2,d1,d2;
  float scale,cosa,sina,rr,vx,vy,time;
  float rf,rt,alfa,beta,gama;
  void to_uniform();
  int nxt,nzt,ntot;
  float *x,*y,*t,*c,*xv,*yv,*tt,*xx,*yy;
  float *xx1,*yy1,*xx2,*yy2;
  int nprocx,nprocz,nx2,nz2;
  double direction,magnitude;

  if (argc < 2)   {
    fprintf(stderr,"Usage: executable PARAMETERFILE\n");
    exit(10);
  }

  frame = atoi(argv[2]);
  nx=atoi(argv[3]);
  nz=atoi(argv[4]);

  sprintf(filename1,"%s/coord.0",argv[1]);
  fp2=fopen(filename1,"r");

  sprintf(filename1,"%s/temp_comp.%d",argv[1],frame);
  fp1=fopen(filename1,"r");
  fgets(input_s,200,fp1);
  sscanf(input_s,"%d %d %g",&nno,&step,&time);

  xx1 = (float *)malloc((nx+1)*sizeof(float));
  yy1 = (float *)malloc((nz+1)*sizeof(float));
  t = (float *)malloc((nno+1)*sizeof(float));
  c = (float *)malloc((nno+1)*sizeof(float));
  xv = (float *)malloc((nno+1)*sizeof(float));
  yv = (float *)malloc((nno+1)*sizeof(float));

  for (n=1;n<=nno;n++)    {
    fgets(input_s,200,fp1);
    sscanf(input_s,"%g %g %g %g",&temp1,&temp2,&temp4,&temp5);
    t[n] = temp1;
    c[n] = temp2;
    xv[n] = temp4;
    yv[n] = temp5;
    }
  fclose(fp1);

  n1 = nx;
  n2 = nz;
  n3 = nx;
  n4 = nz;

  xx2 = (float *)malloc((n3+1)*sizeof(float));
  yy2 = (float *)malloc((n4+1)*sizeof(float));

  fgets(input_s,200,fp2);
  for (i=1;i<=n1;i++)   
  for (j=1;j<=n2;j++)   {
    fgets(input_s,200,fp2);
    sscanf(input_s,"%d %g %g",&n,&temp1,&temp2);
    xx1[i] = temp1;
    yy1[j] = temp2;  
    }
  fclose(fp2);

  for (i=1;i<=n1;i++) { 
       temp1=0;
       for (j=2;j<n2;j++)   {
	  n=j+(i-1)*n2;
	  temp1 += 0.5*(xv[n]+xv[n-1])*(yy1[j]-yy1[j-1]);
//	  c[n] = temp1;
          }
    }

  d1 = (xx1[n1])/(n1-1);
  d2 = (yy1[n2]-yy1[1])/(n2-1);
  for (i=1;i<=n3;i++)   
    xx2[i] = d1*(i-1);  
  for (i=1;i<=n4;i++)   
    yy2[i] = d2*(i-1)+yy1[1];  

  fp2=fopen("temp","w");
  fp3=fopen("strf","w");

  to_uniform(xx1,n1,yy1,n2,xx2,n3,yy2,n4,t);
  to_uniform(xx1,n1,yy1,n2,xx2,n3,yy2,n4,c);

  for (i=1;i<=n3;i++)   
  for (j=1;j<=n4;j++)   {
      n = (i-1)*n4+j;
      fprintf(fp2,"%g %g %g\n",xx2[i],yy2[j],t[n]);
      fprintf(fp3,"%g %g %g\n",xx2[i],yy2[j],c[n]);
    }
  fclose(fp2);
  fclose(fp3);

  
}

void to_uniform(xx1,n1,yy1,n2,xx2,n3,yy2,n4,t)
float *xx1,*yy1,*xx2,*yy2,*t;
int n1,n2,n3,n4;
{

FILE *fp;
float *temp1,*temp2,*temp3,dx1,dx2,area;
int i,j,k1,k2,n,io,jo,ns,ne,ms,me;
static float *temp;
static int been_here=0;

 if(been_here==0)  {
  temp = (float *)malloc((100000+1)*sizeof(float));
  been_here++;
 }

  ne = n1-1;
  me = n2-1;

  for (jo=1;jo<=n2;jo++)    {
    ns = 1;
    for (i=1;i<=n3;i++)    {

      for (io=ns;io<=ne;io++)
        if (xx2[i]<=xx1[io+1]&&xx2[i]>=xx1[io])   {
          ns = io;
          break;
          }

      dx1 = xx2[i]-xx1[ns];
      dx2 = xx1[ns+1]-xx2[i];
      area = dx1+dx2;
      n = jo + (i-1)*n2;
      k1 = jo + (ns-1)*n2;
      k2 = jo + ns*n2;
      temp[n] = (t[k1]*dx2+t[k2]*dx1)/area;
      }
    }

  for (i=1;i<=n3;i++)    {
    ms = 1;
    for (j=1;j<=n4;j++)    {
      for (jo=ms;jo<=me;jo++)
        if (yy2[j]<=yy1[jo+1]&&yy2[j]>=yy1[jo])   {
          ms = jo;
          break;
          }

      dx1 = yy2[j]-yy1[ms];
      dx2 = yy1[ms+1]-yy2[j];
      area = dx1+dx2;
      n = j + (i-1)*n4;
      k1 = ms + (i-1)*n2;
      k2 = ms + (i-1)*n2+1;
      t[n] = (temp[k1]*dx2+temp[k2]*dx1)/area;
      }
    }


 }

