
# include<math.h>
# include<stdio.h>
//#include "stdafx.h"
//#include<stdio.h>
#include<stdlib.h>

#include<math.h>


#define Q         9      //��ɢ�ٶȵ��ܸ���
#define Diameter  20.8    //Բ��ֱ��
#define NX        520//25*Diameter    //x���������
#define NY        208//10*Diameter    //y���������
#define CenterX   208//10*Diameter    //Բ������X����
#define CenterY   (NY/2.0)//    //Բ������Y����


#define U         0.160256    //���������ٶ�
#define	Re        100    //������ŵ��
#define TotNumbs  200000   //��������

extern int e[Q][2];
extern double w[Q];
extern int opp[Q];
extern double rho[NX+1][NY+1],u[NX+1][NY+1],v[NX+1][NY+1],u0[NX+1][NY+1],v0[NX+1][NY+1],f[NX+1][NY+1][Q],ft[NX+1][NY+1][Q];

extern double c,cs2,dx,dy,Lx,Ly,dt,rho0,tau_f, niu,error;

extern int none,fluid,inlet,outlet,wall1,wall2;            //����ʾ
extern int flag[NX+1][NY+1];                         //����ʾ����
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



