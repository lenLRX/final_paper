
# include<math.h>
# include<stdio.h>
//#include "stdafx.h"
//#include<stdio.h>
#include<stdlib.h>

#include<math.h>


#define Q         9      //离散速度的总个数
#define Diameter  20.8    //圆柱直径
#define NX        520//25*Diameter    //x方向格子数
#define NY        208//10*Diameter    //y方向格子数
#define CenterX   208//10*Diameter    //圆柱中心X坐标
#define CenterY   (NY/2.0)//    //圆柱中心Y坐标


#define U         0.160256    //均匀来流速度
#define	Re        100    //流动雷诺数
#define TotNumbs  200000   //迭代次数

extern int e[Q][2];
extern double w[Q];
extern int opp[Q];
extern double rho[NX+1][NY+1],u[NX+1][NY+1],v[NX+1][NY+1],u0[NX+1][NY+1],v0[NX+1][NY+1],f[NX+1][NY+1][Q],ft[NX+1][NY+1][Q];

extern double c,cs2,dx,dy,Lx,Ly,dt,rho0,tau_f, niu,error;

extern int none,fluid,inlet,outlet,wall1,wall2;            //格点标示
extern int flag[NX+1][NY+1];                         //格点标示矩阵
extern double stream[NX+1][NY+1],vorticity[NX+1][NY+1];

void init();
double feq(int k,double rho,double u,double v);
  void collision();
  void streaming(); 
  void boundary();
  void caculateForce(long m);
  void macro();
  void ERROR();



  void output(long m);



