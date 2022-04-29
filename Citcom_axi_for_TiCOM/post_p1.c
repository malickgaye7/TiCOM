/* program to get spectra versus depth on a regular grid */


/* read in velocity near the plume at processor 12, i=27, j=7 for 57h.5530 */

#include <math.h>
#include <stdio.h>
#include <fcntl.h>

 float dd,dens,gg,refT,diff,visc,alfa;

 main(argc,argv)
  int argc;
  char **argv;
 {

  FILE *fp,*fp0,*fp3,*fp1,*fp2,*fp4,*fp5;
  char input_s[100],filename1[250],filename2[250],inf[250],inputf[250],outputf[250];
  int sw,stride,nno,nx,step,i,j,i1,j1,k,jj,nz,nz1,n,m,frame,proc;
 int ny,nprocy,nxy,ny2,nyt,k1,machine,nmachine;
     float rate_1,rate_2,time1;
     int kpkp;
  float temp1,temp2,temp3,temp4,temp5,temp6,tmax,tmin,*t1,*tp,*gp,d1,d2;
  float scale,cosa,sina,rr,vx,vy,time;
  float rate1,rate2,rate3,rate4,xmax;
  float rf,rt,beta,gama;
  float xx1[1000],xx2[1000],yy1[1000],yy2[1000],zz1[100],zz2[100];
 int nxt,nzt,ntot;
float *topo[16],*x[16],*y[16],*z[16],*t[16],*flux_adv[16],*flux_m[16], *tt,*xx,*yy,*zz,*tflux,*tflux_adv,*tflux_m,*width,*area,*f[16],*ff;
//float **flux; // DEBUG
  int nproc,nprocx,nprocz,nx2,nz2;
  float *tp3,*xx3,*reff,*refr,scalef,d_dens,pres;
 int n1,n2,n3,n4,n5,n6;
  float botvel,topvel,botvel1,topvel1,botvel_1,topvel_1,plumeinfo[10],ztmax,xtmax,ttmax,T0,slopemelting,solidus;
  int been;
  void locate_plume();

 void to_uniform();
 void to_uniform_1d();

 const float gradient=0.1;

  
  fp0 = fopen(argv[1],"r");
  fgets(input_s,100,fp0);
  sscanf(input_s,"%s",inputf);
  fgets(input_s,100,fp0);
  sscanf(input_s,"%d %d",&nx,&nz);
  fgets(input_s,100,fp0);
  sscanf(input_s,"%d %d %d %d",&sw,&step,&frame,&kpkp);
  fclose(fp0);

  scale = 1;
  scalef = 1;

  if (argc < 2)   {
          fprintf(stderr,"Usage: executable PARAMETERFILE\n");
          exit(10);
          }

  nxy = nx*nz;

/* read in every thing */

  if (sw>0)    {
      sprintf(filename1,"%s/time_hf_top",inputf);
      fp1=fopen(filename1,"w");
      sprintf(filename1,"%s/time_hf_bot",inputf);
      fp2=fopen(filename1,"w");
      sprintf(filename1,"%s/time_vel_top",inputf);
      fp3=fopen(filename1,"w");
      sprintf(filename1,"%s/time_vel_bot",inputf);
      fp4=fopen(filename1,"w");
      sprintf(filename1,"%s/time_t",inputf);
      fp5=fopen(filename1,"w");
      
      k=frame;
      sprintf(filename1,"%s/ave.%d",inputf,k);
      printf("%s | ",filename1);
      fp0=fopen(filename1,"r");
      temp1 = 0;
      temp2 = 0;
      temp3 = 0;
      botvel1=0;
      topvel1=0;
      do   {
      fgets(input_s,100,fp0);
      sscanf(input_s,"%d %d %g %g %g",&i,&i,&time,&rate1,&rate2);
      temp6=0;
      for (i=1;i<=nz;i++)  {
          fgets(input_s,100,fp0);
          sscanf(input_s,"%g %g %g",&temp5,&rate3,&rate4);
          if (i==nz) topvel=rate4;
          if (i==1) botvel=rate4;
	  temp6 += rate3;
          }
      temp6 = temp6/nz;
      fclose(fp0);
           fprintf(fp1,"%g %g\n",time,rate1);
           fprintf(fp3,"%g %g\n",time,topvel);
           fprintf(fp2,"%g %g\n",time,rate2);
           fprintf(fp4,"%g %g\n",time,botvel);
           fprintf(fp5,"%g %g\n",time,temp6);

      if (k>=kpkp) {
           temp1 += (rate1+rate_1)*0.5*(time-time1);
           temp2 += (rate2+rate_2)*0.5*(time-time1);
           topvel1 += (topvel+topvel_1)*0.5*(time-time1);
           botvel1 += (botvel+botvel_1)*0.5*(time-time1);
           temp3 += (time-time1);
          }
      rate_1 = rate1;
      rate_2 = rate2;
      topvel_1 = topvel;
      botvel_1 = botvel;
      time1 = time;

      k = k + step;
      sprintf(filename1,"%s/ave.%d",inputf,k);
      printf("%s | ",filename1);
      } while ((fp0=fopen(filename1,"r"))!=NULL); 
      printf("\n");
      fprintf(stderr,"averaged hf at top after %d is %g for last %g time\n",kpkp,temp1/temp3,temp3);
      fprintf(stderr,"averaged hf at bottom after %d is %g for last %g time\n",kpkp,temp2/temp3,temp3);
      fprintf(fp1,"averaged hf after %d is %g for last %g time\n",kpkp,temp1/temp3,temp3);
      fprintf(fp2,"averaged hf after %d is %g for last %g time\n",kpkp,temp2/temp3,temp3);

      fprintf(stderr,"averaged vel at top after %d is %g for last %g time\n",kpkp,topvel1/temp3,temp3);
      fprintf(stderr,"averaged vel at bottom after %d is %g for last %g time\n",kpkp,botvel1/temp3,temp3);
      fprintf(fp3,"averaged vel at top after %d is %g for last %g time\n",kpkp,topvel1/temp3,temp3);
      fprintf(fp4,"averaged vel at bottom after %d is %g  for last %g time\n",kpkp,botvel1/temp3,temp3);
      fclose(fp1);
      fclose(fp2);

  }
  else if (sw==0)  {
   i=k=0;
    for (j=0;j<nprocz;j++)    {
      proc = j + i*nprocz + k*nprocz*nprocx;
      nmachine = machine + proc%nproc;
      sprintf(filename1,"/scratch_erd%02d/szhong/%s.coord.%d",nmachine,inf,proc);
      printf("%s\n",filename1);
      fp1=fopen(filename1,"r");
      fgets(input_s,100,fp1);
      sscanf(input_s,"%d",&nno);


      x[proc] = (float *)malloc((nno+1)*sizeof(float));
      y[proc] = (float *)malloc((nno+1)*sizeof(float));
      z[proc] = (float *)malloc((nno+1)*sizeof(float));


      for (n=1;n<=nno;n++)    {
        fgets(input_s,100,fp1);
        sscanf(input_s,"%g %g %g",&temp1,&temp2,&temp3);
        x[proc][n] = temp1;  
        y[proc][n] = temp3;  
        z[proc][n] = 1-temp2;  
        }
      fclose(fp1);

      } 

    for (j=0;j<nprocz;j++)    {
      proc = j;
      flux[proc] = (float *)malloc((nz+1)*sizeof(float));
      flux_adv[proc] = (float *)malloc((nz+1)*sizeof(float));
      flux_m[proc] = (float *)malloc((nz+1)*sizeof(float));
      t[proc] = (float *)malloc((nz+1)*sizeof(float));
      for (i=1;i<=nz;i++)    {
         t[proc][i] = 0;
         flux[proc][i] = 0;
         flux_adv[proc][i] = 0;
         flux_m[proc][i] = 0;
         }
      nmachine = machine + proc%nproc;
      k=frame;
      stride=0;
      sprintf(filename1,"/scratch_erd%02d/szhong/%s.ave.%d.%d",nmachine,inputf,proc,k);
      fp0=fopen(filename1,"r");
      do  {
      fgets(input_s,100,fp0);
      sscanf(input_s,"%d %d %g %g %g",&i,&i,&time,&rate1,&rate2);
      for (i=1;i<=nz;i++)    {
         fgets(input_s,100,fp0);
         sscanf(input_s,"%g %g %g %g %g %g",&temp1,&temp2,&temp2,&temp3,&temp4,&temp5);
         t[proc][i] += temp1;
         flux[proc][i] += temp3;
         flux_adv[proc][i] += temp4;
         flux_m[proc][i] += temp5;
         } 
      fclose(fp0);
      stride = stride+1;
      k = k + step;
      sprintf(filename1,"/scratch_erd%02d/szhong/%s.ave.%d.%d",nmachine,inputf,proc,k);
      printf("%s\n",filename1);
      } while ((fp0=fopen(filename1,"r"))!=NULL); 

      }

  nzt = (nz-1)*nprocz+1;
  nyt = (ny-1)*nprocy+1;
  nxt = (nx-1)*nprocx+1;
  ntot = nxt*nyt*nzt;

  n1 = nxt;
  n2 = nyt;
  n3 = nxt;
  n4 = nyt;
  n5 = nzt;

  xx = (float *)malloc((ntot+1)*sizeof(float));
  yy = (float *)malloc((ntot+1)*sizeof(float));
  zz = (float *)malloc((ntot+1)*sizeof(float));
  tt = (float *)malloc((ntot+1)*sizeof(float));
  tflux = (float *)malloc((ntot+1)*sizeof(float));
  tflux_adv = (float *)malloc((ntot+1)*sizeof(float));
  tflux_m = (float *)malloc((ntot+1)*sizeof(float));

/* combine into one block */

  for (j=0;j<nprocz;j++)    {
      proc = j;
      nz2 = nz-1;
      if(j==nprocz-1) nz2 = nz;

      for (j1=1;j1<=nz2;j1++)   {
        n = j1;
        m = j1 + j*(nz-1);

        zz[m] = z[proc][n];
        tt[m] = t[proc][n]/stride;
        tflux[m] = flux[proc][n]/stride;
        tflux_adv[m] = flux_adv[proc][n]/stride;
        tflux_m[m] = flux_m[proc][n]/stride;
        }
      }

  sprintf(filename1,"%s.aveT.%d",inputf,frame);
  fp2=fopen(filename1,"w");
  sprintf(filename1,"%s.avehf.%d",inputf,frame);
  fp3=fopen(filename1,"w");
  sprintf(filename1,"%s.avehf_adv.%d",inputf,frame);
  fp4=fopen(filename1,"w");
  sprintf(filename1,"%s.avehf_m.%d",inputf,frame);
  fp5=fopen(filename1,"w");
  for (j1=nzt;j1>=1;j1--)   {
        fprintf(fp2,"%g %g\n",zz[j1],tt[j1]);
        fprintf(fp3,"%g %g\n",zz[j1],tflux[j1]);
        fprintf(fp4,"%g %g\n",zz[j1],tflux_adv[j1]);
        fprintf(fp5,"%g %g\n",zz[j1],(tflux_m[j1]+tflux_adv[j1])/(2.0*tflux_adv[j1]));
        }
  fclose(fp2);
  fclose(fp3);
  fclose(fp4);
  fclose(fp5);

  }

 }

 void to_uniform_1d(x1,y1,n1,x2,y2,n2)
 float *x1,*y1,*x2,*y2;
 int n1,n2;
 {

 float dx1,dx2,area;
 int i,k1,k2,io,ns,ne;
 FILE *fp;

 fp=fopen("junk","w");
  ns = 1;
  ne = n1-1;
  for (i=1;i<=n2;i++)    {
     for (io=ns;io<=ne;io++)
        if (x2[i]<=x1[io+1]&&x2[i]>=x1[io])   {
          ns = io;
          break;
          }
      dx1 = x2[i]-x1[ns];
      dx2 = x1[ns+1]-x2[i];
      area = dx1+dx2;
      k1 = ns;
      k2 = ns+1;
      y2[i] = (y1[k1]*dx2+y1[k2]*dx1)/area;
fprintf(fp,"%d %g %g %g %d %d\n",i,area,dx1,dx2,k1,k2);
      }
 fclose(fp);

 return;
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

  void locate_plume(filename1,tt,xx1,zz1,n1,n5,plumeinfo)
  int n1,n5;
  float *tt,*plumeinfo,*zz1,*xx1;
  char filename1[250];
  {

  FILE *fp2;
  int been,i1,j1,n;
  float T0, slopemelting,temp4,temp3,temp2,temp1,pres,solidus;

  T0 = 1050; 
  slopemelting=100/1e9;   /* 100 C per 1e9 Pa */

  been=0;
  fp2=fopen(filename1,"w");

  for (j1=n5;j1>=1;j1--)   {
    temp1 = zz1[j1]*dd/1000;
    pres = dens*gg*temp1*1000;
    solidus = T0 + pres*slopemelting;
    
    temp3=0;
    for (i1=10;i1<=n1-10;i1++)    {
        n = j1 + (i1-1)*n5;
        temp2 = (tt[n-n5] + tt[n] + tt[n+n5])/3;
        if(temp2>temp3)  {
            temp3=temp2;
            temp4 = dd/1000*xx1[i1];
            }
      }
    if (temp3>solidus && been==0) {
       been=1;
       plumeinfo[1] = temp4;
       plumeinfo[2] = temp1;
       plumeinfo[3] = temp3;
       }
    fprintf(fp2,"%.3e %.3e %.3e %.3e %.3e\n",temp1,pres,solidus,temp3,temp4);
    }
  fclose(fp2);
 return;
 }
