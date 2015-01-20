/*************************************************************
 * Numerical Code for Detonation Simulation 
 * Numerical Scheme: CE/SE (unfinished)
 * Boundary Method: Ghost Cell Imersed Boundary Method (developing)
 * Programed by Mo Huangrui
 **************************************************************/
# include<stdio.h>
# include<stdlib.h>
# include<conio.h>
# include<malloc.h>
# include<math.h>
# include<string.h>
# include<stdarg.h>
/*************************************************************
 * data structure for model parameters
 **************************************************************/
typedef struct
{
 double ***U;//storage sapce for vector U
 double ***Un;//storage sapce for U at time n
 double ***F;//storage sapce for Vector F
 double ***G;//storage sapce for Vector G
 double ***Ru;//storage space for primitive variables
 double ***Buff;//buff storage
 double ***Eigen;//eigenvalue storage
 double ***Lmatrix;//eigenvector storage
 double ***Rmatrix;//eigenvector storage
 double ***Wxy;//storage sapce for coordinates of wedges
 FILE *File;//file pointer
 double LX;//model length in x
 double LY;//model length in y
 int MX;//mesh number in x
 int MY;//mesh number in y
 double dx;//mesh size in x
 double dy;//mesh size in y
 double Tt;//total compute time
 double Ma;//Mach number
 double ϒ;//gas constant
 double CFL;//CFL factor
 int WN;//wedge number
}MODEL;
typedef struct
{
 double X_left;
}DOMAIN_GEOMETRY;
/*************************************************************
 * function declaration 
 **************************************************************/
void PreProcess(MODEL *M);
void Solve(MODEL *M);
void PostProcess(MODEL M);
void Initialize(MODEL M);
void Boundary(MODEL M);
void TranslateU(MODEL M);
void CalculateF(MODEL M);
void CalculateG(MODEL M);
double CFL(MODEL M);
void TVDScheme(MODEL M,double dt);
void LoadWedge(MODEL M);
int InWedge(MODEL M,const char *state,int pj,int pi);
void BoundWedge(MODEL M);
void Output(MODEL M,double t);
void GeoGraph(MODEL M);
double Min(double a,double b);
int Max(int a,int b);
double Q(double x);
int sgn(double x);
double minmod(double x,double y,double z);
void Calculateλx(double λx[],double u,double a);
void Calculateλy(double λy[],double v,double a);
void CalculateLRx(double **L,double **R,double u,double v,double a);
void CalculateLRy(double **L,double **R,double u,double v,double a);
void Lx(MODEL M,double dt);
void Ly(MODEL M,double dt);
double ***AssignStorage(MODEL M,int Kd,int Jd,int Id);
void RetrieveStorage(double ***dva,int Kd,int Jd);
void Input(MODEL M,double inf,double sup,const char *format, ...);
void Error(MODEL M,const char *state);
/*************************************************************
 * main function
 **************************************************************/
void main()
{
 //odrp:Two dimension compressible flow
 MODEL TDCF={0};// define model
 PreProcess(&TDCF);// preprocessing
 Solve(&TDCF);//solve
 PostProcess(TDCF);//post processing
}
/*************************************************************
 * preprocessing 
 **************************************************************/
void PreProcess(MODEL *M)
{
 //model parameters(dimensionless)
 printf("\n            Two Dimensional Compressible Flow\n");
 printf("\n                    (By Mo Huangrui)\n");
 printf("\n\n(BOOM! A new life comes..)\n(^.^): Hello,world! Please color my life with the chaos of number~");
 printf("\n\n(o.o) Loading parameters:\n");
 printf("The Lx=");Input(*M,0,1e10,"%lf",&(M->LX));
 printf("The Ly=");Input(*M,0,1e10,"%lf",&(M->LY));
 printf("The meshx=");Input(*M,0,1e10,"%d",&(M->MX));
 printf("The meshy=");Input(*M,0,1e10,"%d",&(M->MY));
 printf("The Mach number=");Input(*M,0,50,"%lf",&(M->Ma));
 printf("The Total time=");Input(*M,0,100,"%lf",&(M->Tt));
 printf("The CFL factor=");Input(*M,0,1,"%lf",&(M->CFL));
 printf("How many wedges needed? Wedge Num=");Input(*M,0,10,"%d",&(M->WN));
 M->dx=M->LX/M->MX;M->dy=M->LY/M->MY;M->ϒ=1.4;
 M->U=AssignStorage(*M,4,M->MY+3,M->MX+3);
 M->Un=AssignStorage(*M,4,M->MY+3,M->MX+3);
 M->F=AssignStorage(*M,4,M->MY+3,M->MX+3);
 M->G=AssignStorage(*M,4,M->MY+3,M->MX+3);
 M->Ru=AssignStorage(*M,5,M->MY+3,M->MX+3);
 M->Buff=AssignStorage(*M,4,M->MY+3,M->MX+3);
 M->Eigen=AssignStorage(*M,4,Max(M->MX,M->MY)+3,4);
 M->Lmatrix=AssignStorage(*M,Max(M->MX,M->MY)+3,4,4);
 M->Rmatrix=AssignStorage(*M,Max(M->MX,M->MY)+3,4,4);
 M->Wxy=AssignStorage(*M,1,10,4);
 if((M->File=fopen("Output.dat","w"))==0){Error(*M,"Open file failed");}
 //geometry
 LoadWedge(*M);
 printf("\n(-.-) Import model succeed...\n");
 printf("\n(~.~) Ahaha,I'm ready now,press any key to let me show you what I got~\n");getch();
}
/*************************************************************
 * solver 
 **************************************************************/
void Solve(MODEL *M)
{
 double time=0,dt=0;
 int step=0,interval=(int)(M->Tt/(20*(M->dx)))+1;//control the time to output 
 Initialize(*M);//initialize
 printf("\n\n(*.*) Start computing...");
 for(step=0,time=0;time<M->Tt;step++)//time marching
 {
  printf("\n**Current time is t=%lf,step %d\n",time,step+1);
  if(step%interval==0){Output(*M,time);}
  dt=CFL(*M);//compute dt
  time+=2*dt;//time marching
  if(time>M->Tt){dt=0.5*(M->Tt+2*dt-time);}//refine "dt" to get time=totaltime
  TVDScheme(*M,dt);
 }
 printf("\n(~.~)! Compute finished~\n");
 TranslateU(*M);//calculate u,v,p and check data
 Output(*M,M->Tt);//output the result at "total time"
 GeoGraph(*M);//draw the wedges and output parameters
}
/*************************************************************
 * post processing 
 **************************************************************/
void PostProcess(MODEL M)
{
 RetrieveStorage(M.Wxy,1,10);
 RetrieveStorage(M.Lmatrix,Max(M.MX,M.MY)+3,4);
 RetrieveStorage(M.Rmatrix,Max(M.MX,M.MY)+3,4);
 RetrieveStorage(M.Eigen,4,Max(M.MX,M.MY)+3);
 RetrieveStorage(M.Buff,4,M.MY+3);
 RetrieveStorage(M.Ru,5,M.MY+3);
 RetrieveStorage(M.G,4,M.MY+3);
 RetrieveStorage(M.F,4,M.MY+3);
 RetrieveStorage(M.Un,4,M.MY+3);
 RetrieveStorage(M.U,4,M.MY+3);
 if(M.File!=0){fclose(M.File);}
 printf("\n(T.T)! My life is burned out,Maybe it's colorful,maybe full of mess,but what matters? Please don't cry for me,if do sad,sing a blue song for me...");
 getch();
}
/*************************************************************
 * load wedge
 **************************************************************/
void LoadWedge(MODEL M)
{
 int wn=0;double *x,*a,*y,*ac,pi=3.1415926;
 for(wn=0;wn<M.WN;wn++)
 {
  x=M.Wxy[0][wn]+0;a=M.Wxy[0][wn]+1;y=M.Wxy[0][wn]+2;ac=M.Wxy[0][wn]+3;
  printf("(o.o)Reading Coordinates for Wedge %d...\n",wn+1);
  printf("X value of the left endpoint=");Input(M,0,M.LX,"%lf",x);
  ac[0]=180+atan(M.LY/(-x[0]))*180/pi;//max angle of input
  printf("Angle(max allowed=%.0lfdegree)=",ac[0]);Input(M,0,ac[0],"%lf",a);
  a[0]=a[0]*pi/180;//translate to rad from degree units
  ac[0]=atan(M.LY/(M.LX-x[0]));//critical angle of wedge
  if(a[0]<ac[0]){y[0]=(M.LX-x[0])*tan(a[0]);}
  else{y[0]=M.LY+2;}
 }
}
/*************************************************************
 * initialize solver 
 **************************************************************/
void Initialize(MODEL M)
{
 double **ρ=M.U[0],**ρu=M.U[1],**ρv=M.U[2],**E=M.U[3];
 double ρ1=1.0,u1=0.0,v1=0.0,p1=0.71429;//initial condition
 int i=0,j=0;
 //initial whole domain
 for(j=0;j<M.MY+3;j++)
 {
  for(i=0;i<M.MX+3;i++)
  {
   ρ[j][i]=ρ1;ρu[j][i]=ρ1*u1;ρv[j][i]=ρ1*v1;
   E[j][i]=p1/(M.ϒ-1)+ρ1*(u1*u1+v1*v1)*0.5;
  }
 }
 Boundary(M);
 TranslateU(M);
}
/*************************************************************
 * boundary condition
 **************************************************************/
void Boundary(MODEL M)
{
 double **ρ=M.U[0],**ρu=M.U[1],**ρv=M.U[2],**E=M.U[3];
 double ρ1=1.0,u1=0,v1=0,p1=0.71429,a1=1,ρ2=0,u2=0,v2=0,p2=0,Ms2=M.Ma*M.Ma;//initial condition
 int i=0,j=0,k=0;
 u2=u1+2*a1*(M.Ma-1/M.Ma)/(M.ϒ+1);v2=0;
 ρ2=ρ1*(M.ϒ+1)*Ms2/(2+(M.ϒ-1)*Ms2);p2=p1*(2*M.ϒ*Ms2-M.ϒ+1)/(M.ϒ+1);
 //bottom wall and up side outlet 
 for(i=2;i<M.MX+1;i++)
 {
  ρ[1][i]=ρ[2][i];ρu[1][i]=ρu[2][i];ρv[1][i]=0;E[1][i]=E[2][i];//slip condition for bottom wall
  for(k=0;k<4;k++){M.U[k][M.MY+1][i]=M.U[k][M.MY][i];}//free-out condition for up side
 }
 //inlet & right side outlet
 for(j=1;j<M.MY+2;j++)
 {
  ρ[j][1]=ρ2;ρu[j][1]=ρ2*u2;ρv[j][1]=ρ2*v2;E[j][1]=p2/(M.ϒ-1)+ρ2*(u2*u2+v2*v2)*0.5;//shock fitting for inlet
  for(k=0;k<4;k++){M.U[k][j][M.MX+1]=M.U[k][j][M.MX];}//free-out condition for outlet
 }
 //virtual mesh
 //virtual bottom wall and up side outlet 
 for(i=1;i<M.MX+2;i++)
 {
  for(k=0;k<4;k++)
  {
   M.U[k][0][i]=2*M.U[k][1][i]-M.U[k][2][i];//linear extrapolate
   M.U[k][M.MY+2][i]=2*M.U[k][M.MY+1][i]-M.U[k][M.MY][i];
  }
 }
 //virtual inlet & right side outlet
 for(j=0;j<M.MY+3;j++)
 {
  for(k=0;k<4;k++)
  {
   M.U[k][j][0]=2*M.U[k][j][1]-M.U[k][j][2];//linear extrapolate
   M.U[k][j][M.MX+2]=2*M.U[k][j][M.MX+1]-M.U[k][j][M.MX];
  }
 }
 BoundWedge(M);
}
void BoundWedge(MODEL M)
{
 double **ρ=M.U[0],**ρu=M.U[1],**ρv=M.U[2],**E=M.U[3];
 int Iw=0,j=0,k=0,wn=0,top=0;
 double *x,*a,*y,yy=0,xx=0,pa=0,pb=0,a1=0,a2=0,a3=0;
 double uu=0,vv=0,pi=3.1415926;
 for(wn=0;wn<M.WN;wn++)
 {
  x=M.Wxy[0][wn]+0;a=M.Wxy[0][wn]+1;y=M.Wxy[0][wn]+2;
  top=int(y[0]/M.dy)+1;if(top>M.MY+1){top=M.MY+1;}
  for(j=1;j<=top;j++)
  {
   yy=(j-1)*M.dy;xx=yy/tan(a[0])+x[0];
   Iw=int(xx/M.dx)+1;pa=2*a[0]/pi;pb=1-pa;
   a1=cos(a[0])*cos(a[0]);
   a2=cos(a[0])*sin(a[0]);
   a3=sin(a[0])*sin(a[0]);
   uu=(pa*ρu[j][Iw-1]+pb*ρu[j+1][Iw]);
   vv=(pa*ρv[j][Iw-1]+pb*ρv[j+1][Iw]);
   ρ[j][Iw]=pa*ρ[j][Iw-1]+pb*ρ[j+1][Iw];
   E[j][Iw]=pa*E[j][Iw-1]+pb*E[j+1][Iw];
   ρu[j][Iw]=uu*a1+vv*a2;
   ρv[j][Iw]=uu*a2+vv*a3;
   //virtual mesh
   for(k=0;k<4;k++){M.U[k][j][Iw+1]=2*M.U[k][j][Iw]-M.U[k][j][Iw-1];}
  }
 }
}
/*************************************************************
 * check whether a node is in one wedge
 **************************************************************/
int InWedge(MODEL M,const char *state,int j,int i)
{
 int wn=0,Iw=0;
 double *x,*a,yy=0,xx=0;
 for(wn=0;wn<M.WN;wn++)
 {
  x=M.Wxy[0][wn]+0;a=M.Wxy[0][wn]+1;
  yy=(j-1)*M.dy;xx=yy/tan(a[0])+x[0];Iw=int(xx/M.dx)+1;
  if(!strcmp(state,"c")){if(i>=Iw){return 1;}}
  else
  {
   if(!strcmp(state,"o")){if(i>Iw+1){return 1;}}
   else
   {
	if(j==0||j==M.MY+2){return 0;}
	else{if((i-3)*M.dx>=xx){return 1;}}
   }
  }
 }
 return 0;
}
/*************************************************************
 * calculate primitive variables from conservative variables 
 **************************************************************/
void TranslateU(MODEL M)
{
 double **ρ=M.Ru[0],**u=M.Ru[1],**v=M.Ru[2],**p=M.Ru[3],**h=M.Ru[4],**ρu=M.U[1],**ρv=M.U[2],**E=M.U[3];
 double sup=1e10,KE=0;
 int i=0,j=0;
 for(j=0;j<M.MY+3;j++)
 {
  for(i=0;i<M.MX+3;i++)
  {
   if(InWedge(M,"w",j,i)){continue;}
   ρ[j][i]=M.U[0][j][i];
   u[j][i]=ρu[j][i]/ρ[j][i];v[j][i]=ρv[j][i]/ρ[j][i];
   KE=0.5*ρ[j][i]*(u[j][i]*u[j][i]+v[j][i]*v[j][i]);
   p[j][i]=(M.ϒ-1)*(E[j][i]-KE);
   h[j][i]=(M.ϒ*p[j][i]/(M.ϒ-1)+KE)/ρ[j][i];
   if((i>0)&&(i<M.MX+2)&&(j>0)&&(j<M.MY+2))
   {
    if((fabs(u[j][i])>sup)||(fabs(v[j][i])>sup))
    {
     printf("\ni=%d,j=%d,ρ=%lf,u=%lf,v=%lf,p=%lf\n",i,j,ρ[j][i],u[j][i],v[j][i],p[j][i]);
     Error(M,"Velocity overflow");
    }
    if((p[j][i]<0)||(ρ[j][i]<0))
    {
     printf("\ni=%d,j=%d,ρ=%lf,u=%lf,v=%lf,p=%lf\n",i,j,ρ[j][i],u[j][i],v[j][i],p[j][i]);
     Error(M,"Negative p or ρ");
    }
   }
  }
 }
}
/*************************************************************
 * compute dt based on CFL condition 
 **************************************************************/
double CFL(MODEL M)
{
 double maxvel=1e-100,vel=0;
 double **ρ=M.Ru[0],**u=M.Ru[1],**v=M.Ru[2],**p=M.Ru[3];
 int i=0,j=0;
 for(j=2;j<M.MY+1;j++)
 {
  for(i=2;i<M.MX+1;i++)
  {
   if(InWedge(M,"c",j,i)){continue;}
   vel=sqrt(M.ϒ*p[j][i]/ρ[j][i])+sqrt(u[j][i]*u[j][i]+v[j][i]*v[j][i]);
   if(vel>maxvel){maxvel=vel;}
  }
 }
 return M.CFL*Min(M.dy,M.dx)/maxvel;
}
/*************************************************************
 * compute flux vector F 
 **************************************************************/
void CalculateF(MODEL M)
{
 double **ρ=M.Ru[0],**u=M.Ru[1],**v=M.Ru[2],**p=M.Ru[3],**ρu=M.U[1],**ρv=M.U[2],**E=M.U[3];
 int i=0,j=0;
 for(j=0;j<M.MY+3;j++)
 {
  for(i=0;i<M.MX+3;i++)
  {
   if(InWedge(M,"w",j,i)){continue;}
   M.F[0][j][i]=ρu[j][i];M.F[1][j][i]=ρu[j][i]*u[j][i]+p[j][i];
   M.F[2][j][i]=ρu[j][i]*v[j][i];M.F[3][j][i]=(E[j][i]+p[j][i])*u[j][i];
  }
 }
}
void CalculateG(MODEL M)
{
 double **ρ=M.Ru[0],**u=M.Ru[1],**v=M.Ru[2],**p=M.Ru[3],**ρu=M.U[1],**ρv=M.U[2],**E=M.U[3];
 int i=0,j=0;
 for(j=0;j<M.MY+3;j++)
 {
  for(i=0;i<M.MX+3;i++)
  {
   if(InWedge(M,"w",j,i)){continue;}
   M.G[0][j][i]=ρv[j][i];M.G[1][j][i]=ρu[j][i]*v[j][i];
   M.G[2][j][i]=ρv[j][i]*v[j][i]+p[j][i];M.G[3][j][i]=(E[j][i]+p[j][i])*v[j][i];
  }
 }
}
/*************************************************************
 * TVD scheme 
 **************************************************************/
void TVDScheme(MODEL M,double dt)
{
 Lx(M,dt/2);Boundary(M);
 Ly(M,dt/2);Boundary(M);
 Ly(M,dt/2);Boundary(M);
 Lx(M,dt/2);Boundary(M);
}
/*************************************************************
 * some math functions 
 **************************************************************/
double Min(double a,double b)
{
 if(a<b){return a;}
 else{return b;}
}
int Max(int a,int b)
{
 if(a>b){return a;}
 else{return b;}
}
double Q(double x)
{
 double ε=0.01;
 if(fabs(x)>=ε){return fabs(x);}
 else{return (0.5*(x*x/ε+ε));}
}
int sgn(double x)
{
 if(x>0){return 1;}
 else
 {
  if(x<0){return -1;}
  else{return 0;}
 }
}
/*************************************************************
 * minmod limiter 
 **************************************************************/
double minmod(double x,double y,double z)
{
 if((x*y<=0)||(x*z)<=0){return 0;}
 else{return (sgn(x)*Min(fabs(x),Min(fabs(y),fabs(z))));}

}
/*************************************************************
 * characteristic values and vectors
 **************************************************************/
void Calculateλx(double λx[],double u,double a)
{
 λx[0]=u;    λx[1]=u;
 λx[2]=u-a;  λx[3]=u+a;
}
void Calculateλy(double λy[],double v,double a)
{
 λy[0]=v;   λy[1]=v;
 λy[2]=v-a; λy[3]=v+a;
}
void CalculateLRx(double **L,double **R,double u,double v,double a)
{
 double b1=0,b2=0,ϒ=1.4,h=0;
 b2=(ϒ-1)/(a*a);b1=0.5*b2*(u*u+v*v);
 h=a*a/(ϒ-1)+0.5*(u*u+v*v);
 L[0][0]=1-b1;           L[0][1]=b2*u;            L[0][2]=b2*v;           L[0][3]=-b2;
 L[1][0]=-b1*v;          L[1][1]=b2*u*v;          L[1][2]=1+b2*v*v;       L[1][3]=-b2*v;
 L[2][0]=0.5*(b1+u/a);   L[2][1]=-0.5*(b2*u+1/a); L[2][2]=-0.5*b2*v;      L[2][3]=0.5*b2;
 L[3][0]=0.5*(b1-u/a);   L[3][1]=-0.5*(b2*u-1/a); L[3][2]=-0.5*b2*v;      L[3][3]=0.5*b2;

 R[0][0]=1;              R[0][1]=0;               R[0][2]=1;              R[0][3]=1;
 R[1][0]=u;              R[1][1]=0;               R[1][2]=u-a;            R[1][3]=u+a;
 R[2][0]=0;              R[2][1]=1;               R[2][2]=v;              R[2][3]=v;
 R[3][0]=0.5*(u*u-v*v);  R[3][1]=v;               R[3][2]=h-a*u;          R[3][3]=h+a*u;
}
void CalculateLRy(double **L,double **R,double u,double v,double a)
{
 double b1=0,b2=0,ϒ=1.4,h=0;
 b2=(ϒ-1)/(a*a);b1=0.5*b2*(u*u+v*v);
 h=a*a/(ϒ-1)+0.5*(u*u+v*v);
 L[0][0]=-b1*u;         L[0][1]=1+b2*u*u;      L[0][2]=b2*u*v;             L[0][3]=-b2*u;
 L[1][0]=1-b1;          L[1][1]=b2*u;          L[1][2]=b2*v;               L[1][3]=-b2;
 L[2][0]=0.5*(b1+v/a);  L[2][1]=-0.5*b2*u;     L[2][2]=-0.5*(b2*v+1/a);    L[2][3]=0.5*b2;
 L[3][0]=0.5*(b1-v/a);  L[3][1]=-0.5*b2*u;     L[3][2]=-0.5*(b2*v-1/a);    L[3][3]=0.5*b2;

 R[0][0]=0;             R[0][1]=1;             R[0][2]=1;                  R[0][3]=1;
 R[1][0]=1;             R[1][1]=0;             R[1][2]=u;                  R[1][3]=u;
 R[2][0]=0;             R[2][1]=v;             R[2][2]=v-a;                R[2][3]=v+a;
 R[3][0]=u;             R[3][1]=0.5*(v*v-u*u); R[3][2]=h-a*v;              R[3][3]=h+a*v;
}
/*************************************************************
 * x direction marching
 **************************************************************/
void Lx(MODEL M,double dt)
{
 double **ρ=M.Ru[0],**u=M.Ru[1],**v=M.Ru[2],**h=M.Ru[4];
 double **λ=M.Eigen[0],**SumRΦ=M.Eigen[1],**α=M.Eigen[2],**g=M.Eigen[3];
 double ***L=M.Lmatrix,***R=M.Rmatrix;
 double Da=0,ua=0,va=0,ha=0,aa=0,r=dt/M.dx,FR=0,FL=0;
 int j=0,i=0,n=0,m=0;
 TranslateU(M);
 CalculateF(M);
 for(j=2;j<M.MY+1;j++)
 {
  for(i=0;i<M.MX+2;i++)
  {
   Da=sqrt(ρ[j][i+1]/ρ[j][i]);
   ua=(u[j][i]+Da*u[j][i+1])/(1+Da);
   va=(v[j][i]+Da*v[j][i+1])/(1+Da);
   ha=(h[j][i]+Da*h[j][i+1])/(1+Da);
   aa=sqrt((M.ϒ-1)*(ha-0.5*(ua*ua+va*va)));
   Calculateλx(λ[i],ua,aa);
   CalculateLRx(L[i],R[i],ua,va,aa);  
  }
  for(i=0;i<M.MX+2;i++)
  {
   for(n=0;n<4;n++)
   {
    α[i][n]=0;
	g[i][n]=0;
    SumRΦ[i][n]=0;
   }
  }
  for(i=0;i<M.MX+2;i++)
  {
   for(n=0;n<4;n++)
   {
    for(m=0;m<4;m++){α[i][n]=α[i][n]+L[i][n][m]*(M.U[m][j][i+1]-M.U[m][j][i]);}
   }
  }
  for(i=1;i<M.MX+1;i++)
  {
   for(n=0;n<4;n++)
   {
   g[i][n]=minmod(α[i-1][n],α[i][n],α[i+1][n]);
   }
  }
  for(i=1;i<M.MX+1;i++)
  {
   for(n=0;n<4;n++)
   {
	for(m=0;m<4;m++){SumRΦ[i][n]=SumRΦ[i][n]+R[i][n][m]*(-1/r)*(pow(r*λ[i][m],2)*g[i][m]+Q(r*λ[i][m])*(α[i][m]-g[i][m]));}
   }
  }
  for(n=0;n<4;n++)
  {
   for(i=2;i<M.MX+1;i++)
   {
	if(InWedge(M,"c",j,i)){continue;}
    FL=0.5*(M.F[n][j][i-1]+M.F[n][j][i]+SumRΦ[i-1][n]);
	FR=0.5*(M.F[n][j][i]+M.F[n][j][i+1]+SumRΦ[i][n]);
	M.U[n][j][i]=M.U[n][j][i]-r*(FR-FL);
   }
  }
 }
}
/*************************************************************
 * y direction marching
 **************************************************************/
void Ly(MODEL M,double dt)
{
 double **ρ=M.Ru[0],**u=M.Ru[1],**v=M.Ru[2],**h=M.Ru[4];
 double **λ=M.Eigen[0],**SumRΦ=M.Eigen[1],**α=M.Eigen[2],**g=M.Eigen[3];
 double ***L=M.Lmatrix,***R=M.Rmatrix;
 double Da=0,ua=0,va=0,ha=0,aa=0,r=dt/M.dy,GR=0,GL=0;
 int i=0,j=0,n=0,m=0;
 TranslateU(M);
 CalculateG(M);
 for(i=2;i<M.MX+1;i++)
 {
  for(j=0;j<M.MY+2;j++)
  {
   Da=sqrt(ρ[j+1][i]/ρ[j][i]);
   ua=(u[j][i]+Da*u[j+1][i])/(1+Da);
   va=(v[j][i]+Da*v[j+1][i])/(1+Da);
   ha=(h[j][i]+Da*h[j+1][i])/(1+Da);
   aa=sqrt((M.ϒ-1)*(ha-0.5*(ua*ua+va*va)));
   Calculateλy(λ[j],va,aa);
   CalculateLRy(L[j],R[j],ua,va,aa);  
  }
  for(j=0;j<M.MY+2;j++)
  {
   for(n=0;n<4;n++)
   {
    α[j][n]=0;
	g[j][n]=0;
    SumRΦ[j][n]=0;
   }
  }
  for(j=0;j<M.MY+2;j++)
  {
   for(n=0;n<4;n++)
   {
    for(m=0;m<4;m++){α[j][n]=α[j][n]+L[j][n][m]*(M.U[m][j+1][i]-M.U[m][j][i]);}
   }
  }
  for(j=1;j<M.MY+1;j++)
  {
   for(n=0;n<4;n++)
   {
   g[j][n]=minmod(α[j-1][n],α[j][n],α[j+1][n]);
   }
  }
  for(j=1;j<M.MY+1;j++)
  {
   for(n=0;n<4;n++)
   {
	for(m=0;m<4;m++){SumRΦ[j][n]=SumRΦ[j][n]+R[j][n][m]*(-1/r)*(pow(r*λ[j][m],2)*g[j][m]+Q(r*λ[j][m])*(α[j][m]-g[j][m]));}
   }
  }
  for(n=0;n<4;n++)
  {
   for(j=2;j<M.MY+1;j++)
   {
	if(InWedge(M,"c",j,i)){continue;}
    GL=0.5*(M.G[n][j-1][i]+M.G[n][j][i]+SumRΦ[j-1][n]);
	GR=0.5*(M.G[n][j][i]+M.G[n][j+1][i]+SumRΦ[j][n]);
	M.U[n][j][i]=M.U[n][j][i]-r*(GR-GL);
   }
  }
 }
}
/*************************************************************
 * output result data
 **************************************************************/
void Output(MODEL M,double t)
{
 double **ρ=M.Ru[0],**u=M.Ru[1],**v=M.Ru[2],**p=M.Ru[3];
 fprintf(M.File,"VARIABLES=\"x\" ,\"y\" , \"density\" , \"u\" , \"v\", \"vel\", \"p\"\n");
 fprintf(M.File,"ZONE T=\"%.6lf s\", I=%d, J=%d, DATAPACKING=POINT\n",t,M.MX+1,M.MY+1);
 double den=0,mu=0,mv=0,mvel=0,mp=0;
 int i=0,j=0;
 for(j=1;j<=M.MY+1;j++)
 {
  for(i=1;i<=M.MX+1;i++)
  {
   if(InWedge(M,"o",j,i)){den=1,mu=0;mv=0;mvel=0;mp=1;}
   else{den=ρ[j][i],mu=u[j][i];mv=v[j][i];mvel=sqrt(mu*mu+mv*mv);mp=p[j][i];}
   fprintf(M.File,"%-20.4lf%-20.4lf%-20.4lf%-20.4lf%-20.4lf%-20.4lf%-20.4lf\n",(i-1)*M.dx,(j-1)*M.dy,den,mu,mv,mvel,mp);
  }
 }
}
/*************************************************************
 * output geometry of wedges
 **************************************************************/
void GeoGraph(MODEL M)
{
 double *x,*y,*a,*ac,xx,zero=0;int i=0;
 fprintf(M.File,"TEXT X=%.2lf,Y=%.2lf,CS=GRID,HU=FRAME,H=2.5,C=BLUE,BX=FILLED,BXF=WHITE,BXO=BLACK,BXM=80,LS=1.5\n",M.LX/6,1.2*M.LY);
 fprintf(M.File,"T=\"Lx=%.2lf, Ly=%.2lf, MeshX=%d, MeshY=%d\\\\nMa=%.2lf, CFL=%.4lf, Wedge Num=%d\"\n",M.LX,M.LY,M.MX,M.MY,M.Ma,M.CFL,M.WN);
 for(i=0;i<M.WN;i++)
 {
  x=M.Wxy[0][i]+0;a=M.Wxy[0][i]+1;y=M.Wxy[0][i]+2;ac=M.Wxy[0][i]+3;
  if(a[0]<ac[0])
  {
   fprintf(M.File,"GEOMETRY X=%.2lf,Y=%.2lf,T=LINE,C=BLACK,FC=BLACK,CS=GRID\n",zero,zero);
   fprintf(M.File,"1\n3\n%-20.4lf%-20.4lf\n%-20.4lf%-20.4lf\n%-20.4lf%-20.4lf\n",x[0],zero,M.LX,y[0],M.LX,zero);
  }
  else
  {
   xx=M.LY/tan(a[0])+x[0];
   fprintf(M.File,"GEOMETRY X=%.2lf,Y=%.2lf,T=LINE,C=BLACK,FC=BLACK,CS=GRID\n",zero,zero);
   fprintf(M.File,"1\n4\n%-20.4lf%-20.4lf\n%-20.4lf%-20.4lf\n%-20.4lf%-20.4lf\n%-20.4lf%-20.4lf\n",x[0],zero,xx,M.LY,M.LX,M.LY,M.LX,zero);
  }
 }
}
/*************************************************************
 * dynamic storage allocate for variables
 **************************************************************/
double *** AssignStorage(MODEL M,int Kd,int Jd,int Id)
{
 double ***temp=0;
 if((temp=(double ***)malloc(Kd*sizeof(double**)))==0){Error(M,"Not enough storage");}
 int k=0,j=0;
 for(k=0;k<Kd;k++)
 {
  if((temp[k]=(double **)malloc(Jd*sizeof(double*)))==0){Error(M,"Not enough storage");}
 }
 for(k=0;k<Kd;k++)
 {
  for(j=0;j<Jd;j++)
  {
   if((temp[k][j]=(double *)malloc(Id*sizeof(double)))==0){Error(M,"Not enough storage");}
  }
 }
 return temp;
}
/*************************************************************
 * storage space release
 **************************************************************/
void RetrieveStorage(double ***dva,int Kd,int Jd)
{
 int k=0,j=0;
 if(dva==0){return;}
 for(k=0;k<Kd;k++)
 {
  for(j=0;j<Jd;j++)
  {
   free((void*)dva[k][j]);
  }
 }
 for(k=0;k<Kd;k++){free((void*)dva[k]);}
 free((void*)dva);
}
/*************************************************************
 * error control
 **************************************************************/
void Error(MODEL M,const char *state)
{
 printf("\n(>_<)! oops! %s\n\n",state);
 PostProcess(M);exit(1);
}
/*************************************************************
 * user input values read in and check
 **************************************************************/
void Input(MODEL M,double inf,double sup,const char *format, ...)
{
 va_list arg_p;
 double temp=0;int buff=0;
 va_start(arg_p,format);
 if(strcmp(format,"%lf")==0){if(!scanf("%lf",va_arg(arg_p,double *))){Error(M,"Illegal input");}}
 else {if(!scanf("%d",va_arg(arg_p,int *))){Error(M,"Illegal input");}}
 va_end(arg_p);
 va_start(arg_p,format);
 if(strcmp(format,"%lf")==0)
 {
  temp=*va_arg(arg_p,double *);
  if((temp>sup)||(temp<inf)){Error(M,"Wrong parametr");}
 }
 else
 {
  buff=*va_arg(arg_p,int *);
  if((buff>sup)||(buff<inf)){Error(M,"Wrong parametr");}
 }
 va_end(arg_p);
}
