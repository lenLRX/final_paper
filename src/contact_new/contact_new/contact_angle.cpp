//����ָ����ʽ��Ч�ܶȣ�EMD������Gong�����������Ӽ�������,�����У�w=2,h=6,�����=3��
//#include "StdAfx.h"
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
#include <time.h>
#include <omp.h>
using namespace std;

//#define H_r                                6 //΢���߶�
#define R                                 20.0//initial(��ʼҺ�ΰ뾶)

#define LX                                101
#define LY                                101
#define LZ                                101
#define Q                                 19

#define X0                                 int((LX-1)/2.0) //Һ�γ�ʼλ��
#define Y0                                 int((LX-1)/2.0)
#define Z0                                 int(R+1)//

#define X1                                 (X0-R-1) //��Һ��
#define Y1                                 Y0
#define Z1                                 Z0
#define X2                                 (X0+R+1) //��Һ��
#define Y2                                 Y0
#define Z1                                 Z0
#define Tcr                                0.0729 //PR���̶�Ӧ���ٽ��¶�
#define T                                  (0.85*Tcr) //��Ӧ�ܶȱȶ�Ӧ���¶�
#define Beta                               1.16//����SC����ʱ���������Ӽ�������������ϲ���
 
#define tau0                               0.60 //��ԥʱ��

#define G                                  -1.0 //���������������ҪС���ٽ�G=-0.222����Ӧ���¶�=-1/G

#define GW0                                -1.56//��������ǿ��

//#define GR                                 -5e-6 //����




#define fluid                            0  //initial_geo
#define solid                            1
#define PI                               3.1415926


#define n0_in                              6.6293//5.9079  //Һ������ٽ��ܶ�ln2=0.693
#define n0_out                             0.3413//0.5801 ////�����ٽ��ܶ�ln2=0.693
#define ERR                                1.0e-6//�����ж�
double R_const=1.0, a=2.0/49.0, b=2.0/21.0, w_pian=0.344, //ˮ��ƫ������		   
alfa_T=(1.0+(0.37464+1.54226*w_pian-0.26992*w_pian*w_pian)*(1.0-sqrt(T/Tcr)))*(1.0+(0.37464+1.54226*w_pian-0.26992*w_pian*w_pian)*(1.0-sqrt(T/Tcr)));
//#define max(x,y)	                      (((x) <= (y)) ? (y) : (x))

/*
int flag[LX][LY][LZ];
double s[LX][LY][LZ], e[Q][3],matrix[Q][Q],
        f0[LX][LY][LZ][Q],n0[LX][LY][LZ],
		veq0x[LX][LY][LZ],veq0y[LX][LY][LZ],veq0z[LX][LY][LZ],vx[LX][LY][LZ],vy[LX][LY][LZ],vz[LX][LY][LZ],
		u0x[LX][LY][LZ],u0y[LX][LY][LZ],u0z[LX][LY][LZ],vx0[LX][LY][LZ],vy0[LX][LY][LZ],vz0[LX][LY][LZ],
		vetical_vel[100000],vetical_positive[100000],vetical_negative[100000],horizontal_vel[100000],horizontal_positive[100000],horizontal_negative[100000],
		A1[100000],A2[100000],bridge_Radius[100000],contactangle[100000],
		f_u[100000],f_v[100000],f_vv[100000],
		g0[LX][LY][LZ][Q],F0x[LX][LY][LZ],F0y[LX][LY][LZ],F0z[LX][LY][LZ],error,mon[LX][LY][LZ][Q],meq[LX][LY][LZ][Q],
	    pressure[LX][LY][LZ];
*/
double  e[Q][3], matrix[Q][Q],contactangle[100000];
double   ****g0, ****mon, ****meq, ****f0, ***n0,***s,***flag,
***psi0,***veq0x, ***veq0y, ***veq0z, ***vx, ***vy, ***vz,
***u0x, ***u0y, ***u0z,
 ***F0x, ***F0y, ***F0z, error,
***pressure;

void init_space()
{
	g0 = new double***[LX];
    mon = new double***[LX];
	meq = new double***[LX];	
	f0 = new double***[LX];
	
	n0 = new double**[LX];
	flag = new double**[LX];
	s = new double**[LX];
	psi0 = new double**[LX];
	veq0x = new double**[LX];
	veq0y = new double**[LX];
	veq0z = new double**[LX];
	vx = new double**[LX];
	vy = new double**[LX];
	vz = new double**[LX];
	u0x = new double**[LX];
	u0y = new double**[LX];
	u0z = new double**[LX];
	F0x = new double**[LX];
	F0y = new double**[LX];
	F0z = new double**[LX];
	pressure = new double**[LX];	

	for (int x = 0; x < LX; x++)
	{
		g0[x] = new double**[LY];
		mon[x] = new double**[LY];
		meq[x] = new double**[LY];		
		f0[x] = new double**[LY];
		
		n0[x] = new double*[LY];
		flag[x] = new double*[LY];
	    s[x] = new double*[LY];
		psi0[x] = new double*[LY];
		veq0x[x] = new double*[LY];
		veq0y[x] = new double*[LY];
		veq0z[x] = new double*[LY];
		vx[x] = new double*[LY];
		vy[x] = new double*[LY];
		vz[x] = new double*[LY];
		u0x[x] = new double*[LY];
		u0y[x] = new double*[LY];
		u0z[x] = new double*[LY];
		F0x[x] = new double*[LY];
		F0y[x] = new double*[LY];
		F0z[x] = new double*[LY];
		pressure[x] = new double*[LY];		

		for (int y = 0; y < LY; y++)
		{
			g0[x][y] = new double*[LZ];
			mon[x][y] = new double*[LZ];
			meq[x][y] = new double*[LZ];			
			f0[x][y] = new double*[LZ];
			

			n0[x][y] = new double[LZ];
			flag[x][y] = new double[LZ];
	        s[x][y] = new double[LZ];
			psi0[x][y] = new double[LZ];
			veq0x[x][y] = new double[LZ];
			veq0y[x][y] = new double[LZ];
			veq0z[x][y] = new double[LZ];
			vx[x][y] = new double[LZ];
			vy[x][y] = new double[LZ];
			vz[x][y] = new double[LZ];
			u0x[x][y] = new double[LZ];
			u0y[x][y] = new double[LZ];
			u0z[x][y] = new double[LZ];
			F0x[x][y] = new double[LZ];
			F0y[x][y] = new double[LZ];
			F0z[x][y] = new double[LZ];
			pressure[x][y] = new double[LZ];			

			for (int z = 0; z < LZ; z++)
			{
				g0[x][y][z] = new double[Q];
				mon[x][y][z] = new double[Q];
				meq[x][y][z] = new double[Q];				
				f0[x][y][z] = new double[Q];
			}
		}
	}
}

double Sv=1.0/tau0,

S0=0.0,
Se=0.3,
SE=0.7,
Sq=1.1,
Spi=1.1,
St=1.8;

/*
S0=0.0,
Se=8.0*(2.0*tau0-1)/(8.0*tau0-1),
SE=8.0*(2.0*tau0-1)/(8.0*tau0-1),
Sq=8.0*(2.0*tau0-1)/(8.0*tau0-1),
Spi=8.0*(2.0*tau0-1)/(8.0*tau0-1),
St=8.0*(2.0*tau0-1)/(8.0*tau0-1);
*/
//s={0.0,Se,SE,0.0,Sq,0.0,Sq,0.0,Sq,Sv,Spi,Sv,Spi,Sv,Sv,Sv,St,St,St}; 
//�ɳھ���S
double S_tau[Q][Q]={
{S0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, Se, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, SE, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, S0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, Sq, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, S0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Sq, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, S0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Sq, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,Sv, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, Spi, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, Sv, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, Spi, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, Sv, 0.0, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, Sv, 0.0, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Sv, 0.0, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, St, 0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, St, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, St}
};


double gravity0,angle,time0,time1;
double w[Q]={1.0/3,1.0/18,1.0/18,1.0/18,1.0/18,1.0/18,1.0/18,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36}; 
int r[Q]={0,2,1,4,3,6,5,10,9,8,7,14,13,12,11,18,17,16,15}; 

double M[Q][Q]={
{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,1.0, 1.0},
{-30.0, -11.0, -11.0, -11.0, -11.0, -11.0, -11.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0,8.0, 8.0},
{12.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,1.0, 1.0},
{0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0,1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0,0.0, 0.0},
{0.0, -4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0,1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0,0.0, 0.0},
{0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, 1.0,-1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0,1.0, -1.0},
{0.0, 0.0, 0.0, -4.0, 4.0, 0.0, 0.0, 1.0, 1.0,-1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0,1.0, -1.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0,0.0, 0.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0,-1.0, -1.0},
{0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 4.0, 0.0, 0.0,0.0, 0.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0,-1.0, -1.0},
{0.0, 2.0, 2.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0,1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -2.0, -2.0,-2.0,-2.0},
{0.0, -4.0, -4.0, 2.0, 2.0, 2.0, 2.0, 1.0, 1.0,1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -2.0, -2.0,-2.0,-2.0},
{0.0, 0.0, 0.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0,1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 0.0, 0.0,0.0,0.0},
{0.0, 0.0, 0.0, -2.0, -2.0, 2.0, 2.0, 1.0, 1.0,1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 0.0, 0.0,0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0,-1.0, 1.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0,0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0,-1.0, 1.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0,0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0,1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 0.0, 0.0,0.0, 0.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0,1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0,1.0,-1.0},
{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0,1.0, 1.0}
};

double M_inv[Q][Q]={ 
{1.0/19.0, -5.0/399.0, 1.0/21.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0},
{1.0/19.0, -11.0/2394.0, -1.0/63.0, 1.0/10.0, -1.0/10.0, 0.0, 0.0, 0.0, 0.0,1.0/18.0, -1.0/18.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0},
{1.0/19.0, -11.0/2394.0, -1.0/63.0, -1.0/10.0, 1.0/10.0, 0.0, 0.0, 0.0, 0.0,1.0/18.0, -1.0/18.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0},
{1.0/19.0, -11.0/2394.0, -1.0/63.0, 0.0, 0.0, 1.0/10.0, -1.0/10.0, 0.0, 0.0,-1.0/36.0, 1.0/36.0, 1.0/12.0, -1.0/12.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0},
{1.0/19.0, -11.0/2394.0, -1.0/63.0, 0.0, 0.0, -1.0/10.0,1.0/10.0, 0.0, 0.0,-1.0/36.0, 1.0/36.0, 1.0/12.0, -1.0/12.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0},
{1.0/19.0, -11.0/2394.0, -1.0/63.0, 0.0, 0.0, 0.0, 0.0,1.0/10.0, -1.0/10.0,-1.0/36.0, 1.0/36.0, -1.0/12.0, 1.0/12.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0},
{1.0/19.0, -11.0/2394.0, -1.0/63.0, 0.0, 0.0, 0.0, 0.0,-1.0/10.0, 1.0/10.0,-1.0/36.0, 1.0/36.0, -1.0/12.0,  1.0/12.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0},
{1.0/19.0, 4.0/1197.0, 1.0/252.0, 1.0/10.0, 1.0/40.0, 1.0/10.0,1.0/40.0,0.0, 0.0,1.0/36.0, 1.0/72.0, 1.0/12.0, 1.0/24.0, 1.0/4.0, 0.0, 0.0, 1.0/8.0,-1.0/8.0, 0.0},
{1.0/19.0, 4.0/1197.0, 1.0/252.0, -1.0/10.0, -1.0/40.0, 1.0/10.0,1.0/40.0,0.0, 0.0,1.0/36.0, 1.0/72.0, 1.0/12.0, 1.0/24.0, -1.0/4.0, 0.0, 0.0, -1.0/8.0,-1.0/8.0, 0.0},
{1.0/19.0, 4.0/1197.0, 1.0/252.0, 1.0/10.0, 1.0/40.0, -1.0/10.0,-1.0/40.0,0.0, 0.0,1.0/36.0, 1.0/72.0, 1.0/12.0, 1.0/24.0, -1.0/4.0, 0.0, 0.0, 1.0/8.0,1.0/8.0, 0.0},
{1.0/19.0, 4.0/1197.0, 1.0/252.0, -1.0/10.0, -1.0/40.0, -1.0/10.0,-1.0/40.0,0.0, 0.0,1.0/36.0, 1.0/72.0, 1.0/12.0,1.0/24.0, 1.0/4.0, 0.0, 0.0, -1.0/8.0,1.0/8.0, 0.0},
{1.0/19.0, 4.0/1197.0, 1.0/252.0, 1.0/10.0, 1.0/40.0, 0.0,0.0,1.0/10.0, 1.0/40.0,1.0/36.0, 1.0/72.0, -1.0/12.0, -1.0/24.0, 0.0, 0.0,1.0/4.0, -1.0/8.0,0.0, 1.0/8.0},
{1.0/19.0, 4.0/1197.0, 1.0/252.0, -1.0/10.0, -1.0/40.0, 0.0,0.0,1.0/10.0, 1.0/40.0,1.0/36.0, 1.0/72.0, -1.0/12.0,-1.0/24.0, 0.0, 0.0,-1.0/4.0, 1.0/8.0,0.0, 1.0/8.0},
{1.0/19.0, 4.0/1197.0, 1.0/252.0, 1.0/10.0, 1.0/40.0, 0.0,0.0,-1.0/10.0, -1.0/40.0,1.0/36.0, 1.0/72.0, -1.0/12.0, -1.0/24.0, 0.0, 0.0,-1.0/4.0, -1.0/8.0,0.0, -1.0/8.0},
{1.0/19.0, 4.0/1197.0, 1.0/252.0, -1.0/10.0, -1.0/40.0, 0.0,0.0,-1.0/10.0, -1.0/40.0,1.0/36.0, 1.0/72.0, -1.0/12.0, -1.0/24.0, 0.0, 0.0,1.0/4.0, 1.0/8.0,0.0, -1.0/8.0},
{1.0/19.0, 4.0/1197.0, 1.0/252.0, 0.0, 0.0, 1.0/10.0,1.0/40.0,1.0/10.0, 1.0/40.0,-1.0/18.0, -1.0/36.0, 0.0, 0.0, 0.0, 1.0/4.0,0.0, 0.0,1.0/8.0, -1.0/8.0},
{1.0/19.0, 4.0/1197.0, 1.0/252.0, 0.0, 0.0, -1.0/10.0, -1.0/40.0,1.0/10.0, 1.0/40.0,-1.0/18.0, -1.0/36.0, 0.0, 0.0, 0.0, -1.0/4.0,0.0, 0.0,-1.0/8.0, -1.0/8.0},
{1.0/19.0, 4.0/1197.0, 1.0/252.0, 0.0, 0.0, 1.0/10.0,1.0/40.0,-1.0/10.0, -1.0/40.0,-1.0/18.0, -1.0/36.0, 0.0, 0.0, 0.0, -1.0/4.0,0.0, 0.0,1.0/8.0, 1.0/8.0},
{1.0/19.0, 4.0/1197.0, 1.0/252.0, 0.0, 0.0, -1.0/10.0,-1.0/40.0,-1.0/10.0, -1.0/40.0,-1.0/18.0, -1.0/36.0, 0.0, 0.0,0.0, 1.0/4.0,0.0, 0.0,-1.0/8.0, 1.0/8.0}
};



double eq(int k,double n,double ux,double uy,double uz);
void ini_geo(void);
void initial(void);
void collision(void);
void stream(void);
void para(void);
void contact_angle(int m);
void evolve(void);
double Error();
void output(int m);
void output_contant_angle(int m);

void vetical_velocity(int m);
void output_vetical_velocity(int m);
void efficient_kinetic(int m);
void bridge_radius(int m);



double eq(int k,double n,double ux,double uy,double uz)
{
	double eu,uv,feq;
	eu=(e[k][0]*ux+e[k][1]*uy+e[k][2]*uz);
	uv=(ux*ux+uy*uy+uz*uz);
	feq=w[k]*n*(1.0+3.0*eu+4.5*eu*eu-1.5*uv);
	return feq;
}
//ini_geo
void ini_geo()
{
	int i,j,z,m,n;

	for(j=0;j<LY;j++)
		for(i=0;i<LX;i++)
		{                     
			
				flag[i][j][0]=solid; //�±���
				s[i][j][0]=1.0;
				flag[i][j][LZ-1]=solid;//�ϱ���
				s[i][j][LZ-1]=1.0;
			
		}   


for(z=1;z<LZ-1;z++)  //�������
	for(j=0;j<LY;j++)
		for(i=0;i<LX;i++)
		{
			if(flag[i][j][z]!=solid)
			{
				flag[i][j][z]=fluid;
				s[i][j][z]=0.0;
			}

		}


}

//initial
void initial()
{
	int k1, k2, k3;
	for(k1=0; k1<Q; k1++)
	{
		for(k2=0; k2<Q; k2++)
		{
			matrix[k1][k2] = 0.0;
			for(k3=0; k3<Q; k3++)
			{
				matrix[k1][k2] += M_inv[k1][k3] * S_tau[k3][k2];
			}
		}
	}


	int i,j,k,z;
#pragma omp parallel for private(i,j,z,k) num_threads(2)
 for(z=0;z<LZ;z++)
	for(j=0;j<LY;j++)
		for(i=0;i<LX;i++)//ע�� i �ķ�Χ(ֻ�����岿�ֳ�ʼ��)			
		{
						
				if((double((i-X0)*(i-X0)+(j-Y0)*(j-Y0)+(z-Z0)*(z-Z0)))<=R*R)
				{
					n0[i][j][z]=n0_in;
					
				}
			/*	else if((double((i-X0-2*R-15)*(i-X0-2*R-15)+(j-Y0)*(j-Y0)))<=R*R)
				{
					n0[i][j]=n0_in;
					
				}
			*/
				else
				{
					n0[i][j][z]=n0_out;
					
				}	
				veq0x[i][j][z]=0.0;
				veq0y[i][j][z]=0.0;
				veq0z[i][j][z]=0.0;
			    vx[i][j][z]=0.0; //����error�ж�
			    vy[i][j][z]=0.0;
				vz[i][j][z]=0.0;

				for(k=0;k<Q;k++)
				{
					f0[i][j][z][k]=eq(k,n0[i][j][z],veq0x[i][j][z],veq0y[i][j][z],veq0z[i][j][z]);

				}

		}
		//para();

}
//collision
void collision()
{
	int i,j,k,z,k1,k2;
    double temp0,temp1,temp2,jx,jy,jz,temp3;

#pragma omp parallel for private(i,j,z,k,jx,jy,jz) num_threads(2)
for(z=0;z<LZ;z++)
	for(j=0;j<LY;j++)
		for(i=0;i<LX;i++)
		{
			if(flag[i][j][z]==fluid)
			{
			    jx=n0[i][j][z]*veq0x[i][j][z];
                jy=n0[i][j][z]*veq0y[i][j][z];
                jz=n0[i][j][z]*veq0z[i][j][z];
                
                meq[i][j][z][0]=  n0[i][j][z];
                meq[i][j][z][1]=  -11.0*n0[i][j][z]+19.0/n0[i][j][z]*(jx*jx + jy*jy+jz*jz);
                meq[i][j][z][2]=  3.0*n0[i][j][z]-11.0/n0[i][j][z]/2.0*(jx*jx + jy*jy+jz*jz);
                meq[i][j][z][3]=  jx;
                meq[i][j][z][4]=  -2.0/3.0*jx;
                meq[i][j][z][5]=  jy;
                meq[i][j][z][6]=  -2.0/3.0*jy;
                meq[i][j][z][7]=  jz;
                meq[i][j][z][8]=  -2.0/3.0*jz;                
                meq[i][j][z][9]=  (2.0*jx*jx-jy*jy-jz*jz)/n0[i][j][z];
                meq[i][j][z][10]=  -0.5*(2.0*jx*jx-jy*jy-jz*jz)/n0[i][j][z];
                meq[i][j][z][11]=  (jy*jy-jz*jz)/n0[i][j][z];
                meq[i][j][z][12]=  -0.5*(jy*jy-jz*jz)/n0[i][j][z];
                meq[i][j][z][13]=  jx*jy/n0[i][j][z];
                meq[i][j][z][14]=  jy*jz/n0[i][j][z];
                meq[i][j][z][15]=  jx*jz/n0[i][j][z];
                meq[i][j][z][16]=  0.0;
                meq[i][j][z][17]=  0.0;
                meq[i][j][z][18]=  0.0;		
			}
			else
			{
				for(k=0;k<Q;k++)
				{
					meq[i][j][z][k]=0.0;
				}
			}
		}

#pragma omp parallel for private(i,j,z,k1,k2,temp2) num_threads(2)
for(z=0;z<LZ;z++)
	for(j=0;j<LY;j++)
		for(i=0;i<LX;i++)
		{
			if(flag[i][j][z]==fluid)
			{
		   	for(k1=0;k1<Q;k1++)		
				{	
					temp2 = 0.0;
					for (k2=0; k2<Q; k2++)
					{						
						temp2 = temp2 + M[k1][k2]*f0[i][j][z][k2];
					}
					mon[i][j][z][k1]=temp2;
				}
			}
			else
			{
				for(k1=0;k1<Q;k1++)
				{
					mon[i][j][z][k1]=0.0;
				}
			}
		}
#pragma omp parallel for private(i,j,z,k1,k2,temp0,temp1,temp3) num_threads(2)
  for(z=0;z<LZ;z++)
	for(j=0;j<LY;j++)
		for(i=0;i<LX;i++)
		{
			if(flag[i][j][z]==fluid)
			{
			for(k1=0;k1<Q;k1++)		
				{	
					temp3 = 0.0;
					temp0=eq(k1,n0[i][j][z],veq0x[i][j][z],veq0y[i][j][z],veq0z[i][j][z]);
					temp1=eq(k1,n0[i][j][z],u0x[i][j][z],u0y[i][j][z],u0z[i][j][z]); //EDM������������������ɵ��ܶȷֲ���������
					for (k2=0; k2<Q; k2++)
					{						
						temp3 += matrix[k1][k2] * (mon[i][j][z][k2] - meq[i][j][z][k2]);
					}
					g0[i][j][z][k1]=f0[i][j][z][k1]-temp3+temp1-temp0;
					
				}
			}
			else
			{
				for(k1=0;k1<Q;k1++)
				{
					g0[i][j][z][k1]=0.0;
				}
			}
		}

/*
#pragma omp parallel for private(i,j,z,k,temp0,temp1) num_threads(2)
  for(z=0;z<LZ;z++)
	for(j=0;j<LY;j++)
		for(i=0;i<LX;i++)
		{
			if(flag[i][j][z]==fluid)
			{
				for(k=0;k<Q;k++)
				{
					temp0=eq(k,n0[i][j][z],veq0x[i][j][z],veq0y[i][j][z],veq0z[i][j][z]);
					temp1=eq(k,n0[i][j][z],u0x[i][j][z],u0y[i][j][z],u0z[i][j][z]); //EDM������������������ɵ��ܶȷֲ���������

					g0[i][j][z][k]=f0[i][j][z][k]-(f0[i][j][z][k]-temp0)/tau0+temp1-temp0;

				}
			}
			else
			{
				for(k=0;k<Q;k++)
				{
					g0[i][j][z][k]=0.0;
				}
			}
		}
*/		

}
//stream
void stream()
{
	int i,j,z,k,id,jd,zd;
//Ǩ�Ʋ����ò���
#pragma omp parallel for private(i,j,z,k,id,jd,zd) num_threads(2)
 for(z=0;z<LZ;z++)
	for(i=0;i<LX;i++)  
		for(j=0;j<LY;j++)
		{
			if(flag[i][j][z]==fluid)	
			{
				for(k=0;k<Q;k++)
				{
					id=i-int(e[k][0]);jd=j-int(e[k][1]);zd=z-int(e[k][2]);
					if(id>LX-1) id=0;	if(id<0) id=LX-1;
					if(jd>LY-1) jd=0;	if(jd<0) jd=LY-1;

					if(flag[id][jd][zd]==solid)	
						f0[i][j][z][k]=g0[i][j][z][r[k]];
					else if(flag[id][jd][zd]!=solid)	
						f0[i][j][z][k]=g0[id][jd][zd][k];					
				}
			}
		}
	
}
//para
void para()
{

		int i,j,k,z;
		double temp0;
		double psi0[LX][LY][LZ];

#pragma omp parallel for private(i,j,z,k,temp0) num_threads(2)
		for(z=0;z<LZ;z++)
			for(j=0;j<LY;j++)
				for(i=0;i<LX;i++)
				{	
					if(flag[i][j][z]==fluid)
					{
					temp0=0.0;
					for(k=0;k<Q;k++)
					{
						temp0+=f0[i][j][z][k];

					}
					n0[i][j][z]=temp0; //����ܶ�

				pressure[i][j][z]=n0[i][j][z]*R_const*T/(1.0-b*n0[i][j][z])-a*alfa_T*n0[i][j][z]*n0[i][j][z]/(1.0+2.0*b*n0[i][j][z]-b*b*n0[i][j][z]*n0[i][j][z]);  //gȡֵG������
				psi0[i][j][z]=sqrt(2.0*(pressure[i][j][z]-n0[i][j][z]/3.0)*3/G); 
					}
					else
					{
						n0[i][j][z]=0.0;
						psi0[i][j][z]=0.0;
						pressure[i][j][z]=0.0;
					}
				}		
 

	

	int id,jd,zd;
	double Fx_temp,Fy_temp,Fz_temp,Fxx_temp,Fyy_temp,Fzz_temp,sum_x,sum_y,sum_z;

#pragma omp parallel for private(i,j,k,z,Fx_temp,Fy_temp,Fz_temp,Fxx_temp,Fyy_temp,Fzz_temp,sum_x,sum_y,sum_z,id,jd,zd) num_threads(2)
	for(z=1;z<LZ-1;z++)
		for(j=0;j<LY;j++)
			for(i=0;i<LX;i++)
		{
		Fx_temp=0.0;
        Fy_temp=0.0;
        Fz_temp=0.0;
        Fxx_temp=0.0;
        Fyy_temp=0.0;
        Fzz_temp=0.0;
		sum_x=0.0;
        sum_y=0.0;
		sum_z=0.0;

			if(flag[i][j][z]==fluid)
			{
				for(k=1;k<Q;k++)
				{
					id=i+int(e[k][0]);jd=j+int(e[k][1]);zd=z+int(e[k][2]);
					if(id>LX-1) id=0; if(id<0) id=LX-1;
					if(jd>LY-1) jd=0; if(jd<0) jd=LY-1;					

					  Fx_temp= Fx_temp+w[k]*e[k][0]*psi0[id][jd][zd];
                      Fxx_temp= Fxx_temp+w[k]*e[k][0]*psi0[id][jd][zd]*psi0[id][jd][zd];

                      Fy_temp= Fy_temp+w[k]*e[k][1]*psi0[id][jd][zd];
                      Fyy_temp= Fyy_temp+w[k]*e[k][1]*psi0[id][jd][zd]*psi0[id][jd][zd];

					  Fz_temp= Fz_temp+w[k]*e[k][2]*psi0[id][jd][zd];
                      Fzz_temp= Fzz_temp+w[k]*e[k][2]*psi0[id][jd][zd]*psi0[id][jd][zd];

					  sum_x=sum_x+w[k]*e[k][0]*s[id][jd][zd];
					  sum_y=sum_y+w[k]*e[k][1]*s[id][jd][zd];
					  sum_z=sum_z+w[k]*e[k][2]*s[id][jd][zd];		
               
			   }
			F0x[i][j][z]=-Beta*G*psi0[i][j][z]*Fx_temp-(1-Beta)/2.0*Fxx_temp*G-GW0*n0[i][j][z]*sum_x;
			F0y[i][j][z]=-Beta*G*psi0[i][j][z]*Fy_temp-(1-Beta)/2.0*Fyy_temp*G-GW0*n0[i][j][z]*sum_y; 
			F0z[i][j][z]=-Beta*G*psi0[i][j][z]*Fz_temp-(1-Beta)/2.0*Fzz_temp*G-GW0*n0[i][j][z]*sum_z;  //+GR*(n0[i][j][z]-0.0596);  //F0y[i][j]=f0y+gravity0*n0[i][j]*m0; ������
			}
		}



double mome0x,mome0y,mome0z;

#pragma omp parallel for private(i,j,k,z,mome0x,mome0y,mome0z) num_threads(2)
	for(z=0;z<LZ;z++)
		for(j=0;j<LY;j++)
			for(i=0;i<LX;i++)
		{
			if(flag[i][j][z]==fluid)
			{
				//vx0[i][j][z]=vx[i][j][z]; //������һʱ����ʵ�ٶȣ����������ж�
				//vy0[i][j][z]=vy[i][j][z];
				//vz0[i][j][z]=vz[i][j][z];
				mome0x=mome0y=mome0z=0.0;
				for(k=0;k<Q;k++)
				{
					mome0x+=f0[i][j][z][k]*e[k][0]; //x������
					mome0y+=f0[i][j][z][k]*e[k][1]; //y������
					mome0z+=f0[i][j][z][k]*e[k][2]; //z������
					
				}
				//mome0x=mome0x;mome0y=mome0y;mome0z=mome0z;				

				//u0x[i][j]=(mome0x/den0[i][j]); // SC���������ƽ��̬�ٶ�
				//u0y[i][j]=(mome0y/den0[i][j]);
				veq0x[i][j][z]=(mome0x/n0[i][j][z]); // EDM���������ƽ��̬�ٶ�
				veq0y[i][j][z]=(mome0y/n0[i][j][z]);
				veq0z[i][j][z]=(mome0z/n0[i][j][z]);

				//veq0x[i][j]=(u0x[i][j]+tau0*F0x[i][j]/den0[i][j]); // SC����������
				//veq0y[i][j]=(u0y[i][j]+tau0*F0y[i][j]/den0[i][j]);
				u0x[i][j][z]=(veq0x[i][j][z]+F0x[i][j][z]/n0[i][j][z]); //EDM����������
				u0y[i][j][z]=(veq0y[i][j][z]+F0y[i][j][z]/n0[i][j][z]);
				u0z[i][j][z]=(veq0z[i][j][z]+F0z[i][j][z]/n0[i][j][z]);

				//vx[i][j]=(u0x[i][j]+0.5*F0x[i][j]/den0[i][j]); //SC�����������ʵ�ٶ�
				//vy[i][j]=(u0y[i][j]+0.5*F0y[i][j]/den0[i][j]);
				vx[i][j][z]=(veq0x[i][j][z]+0.5*F0x[i][j][z]/n0[i][j][z]); //EDM�����������ʵ�ٶ�
				vy[i][j][z]=(veq0y[i][j][z]+0.5*F0y[i][j][z]/n0[i][j][z]);
				vz[i][j][z]=(veq0z[i][j][z]+0.5*F0z[i][j][z]/n0[i][j][z]);
			}
			else 
			{
				veq0x[i][j][z]=0.0;
				veq0y[i][j][z]=0.0;
				veq0z[i][j][z]=0.0;
				vx[i][j][z]=0.0;
				vy[i][j][z]=0.0;
				vz[i][j][z]=0.0;
			}
		}

			

}

/*
void bridge_radius(int m)
{
	double temp1,temp2;
	int j;
    double Rho_interface=(n0_in+n0_out)/2;
	for(j=1;j<LY-2;j++)
	 if(flag[X0][j]==fluid)
	{
		temp1=n0[X0][j];
		temp2=n0[X0][j+1];
		
		if(temp1<Rho_interface&&temp2>=Rho_interface)
			A1[m]=double(j)+(Rho_interface-temp1)/(temp2-temp1);
		if(temp1>Rho_interface&&temp2<=Rho_interface)
			A2[m]=double(j)+(Rho_interface-temp1)/(temp2-temp1);
	}
	bridge_Radius[m]=A2[m]-A1[m];

}

void vetical_velocity(int m) //������ֱ�����ٶ�
{
	double temp1=0.0,temp2=0.0;//,temp3=0.0,temp4=0.0,temp5=0.0,temp6=0.0,temp7=0.0;
	int i,j,z;

for(z=0;z<LZ;z++)
	for(j=0;j<LY;j++) 
		for(i=0;i<LX;i++)
		{
			if(flag[i][j][z]==fluid)
			{
				if(n0[i][j][z]>=(n0_in+n0_out)/2)
				{
					temp1+=n0[i][j][z]*vz[i][j][z];
					//temp5+=n0[i][j]*vx[i][j];
					temp2+=n0[i][j][z];
					//{if(vy[i][j]>=0) temp3+=n0[i][j]*vy[i][j];	else if(vy[i][j]<0) temp4+=n0[i][j]*vy[i][j];}
					//{if(vx[i][j]>=0) temp6+=n0[i][j]*vx[i][j];	else if(vx[i][j]<0) temp7+=n0[i][j]*vx[i][j];}
					
				}
			}
		}
	vetical_vel[m]=temp1/temp2;
	//vetical_positive[m]=temp3/temp2;
	//vetical_negative[m]=temp4/temp2;

	//horizontal_vel[m]=temp5/temp2;
	//horizontal_positive[m]=temp6/temp2;
	//horizontal_negative[m]=temp7/temp2;

}


void efficient_kinetic(int m) //������Ч����
{
	double temp1=0.0,temp2=0.0,temp3=0.0;
	int i,j;
	for(j=0;j<LY;j++) 
		for(i=0;i<LX;i++)
		{
			if(flag[i][j]==fluid)
			{
				if(n0[i][j]>=(n0_in+n0_out)/2)
				{
					temp1+=n0[i][j]*vx[i][j]*vx[i][j];
					temp2+=n0[i][j]*vy[i][j]*vy[i][j];
					if(vy[i][j]>0) temp3+=n0[i][j]*vy[i][j]*vy[i][j];
				}
			}
		}
		f_u[m]=temp1/(temp1+temp2);
		f_v[m]=temp2/(temp1+temp2);
		f_vv[m]=temp3/(temp1+temp2);		

}
*/
void contact_angle(int m)  //����Ӵ���
{

	double temp1,temp2;
	int i,z;
	double b0,b1=0.0,b2=0.0,a0,a1=0.0,Rho_interface=(n0_in+n0_out)/2,  //���ܶȵ��ڣ�Rho_Һ+Rho_����/2��Ϊ���棬�˴���j=2ʱ���Ŵ���R_interface�ĵ�
		Radius;
	//double angle_former=angle,ratio; //��¼��һ�νӴ��ǣ������жϱ仯

	for(z=1;z<LZ-2;z++)
	{
		temp1=n0[X0][Y0][z];
		temp2=n0[X0][Y0][z+1];

		if(temp1>Rho_interface&&temp2<=Rho_interface)
			a1=double(z)+(Rho_interface-temp1)/(temp2-temp1);
	}
	a0=a1-2.0;  //��ڸ߶�

	for(i=0;i<LX-1;i++)
	{
		temp1=n0[i][Y0][2];
		temp2=n0[i+1][Y0][2];

		if(temp1<Rho_interface&&temp2>=Rho_interface)
			b1=double(i)+(Rho_interface-temp1)/(temp2-temp1);
		if(temp1>Rho_interface&&temp2<=Rho_interface)
			b2=double(i)+(Rho_interface-temp1)/(temp2-temp1);
	}
	b0=b2-b1; //��y=Y0�ص�����


	Radius=a0/2.0+b0*b0/(8.0*a0);

	angle=atan(b0/(2.0*Radius-2.0*a0));

	if (angle<0)
		angle=PI+angle;

	angle=angle/PI*180; //תΪ�Ƕ�
	contactangle[m]=angle;
	//ratio=fabs((angle-angle_former)/angle_former);  //�Ӵ��Ǳ仯��

	//printf("�뾶=%f\t\t�Ӵ���=%f\t\t�仯��=%f%%\t\n",Radius,angle,ratio*100);

}

//evolve
void evolve()
{
	collision();
	stream();
	para();
	
}

/*
double Error()
{
	int i,j,z;
	double temp1,temp2,err;
	temp1=0;
	temp2=0;

//#pragma omp parallel for reduction(+:temp1,temp2)
	for(z=1;z<LZ-1;z++)
		for(i=0;i<LX;i++)
			for(j=0;j<LY;j++)
			{
				temp1+=(vx[i][j][z]-vx0[i][j][z])*(vx[i][j][z]-vx0[i][j][z])+(vy[i][j][z]-vy0[i][j][z])*(vy[i][j][z]-vy0[i][j][z])+(vz[i][j][z]-vz0[i][j][z])*(vz[i][j][z]-vz0[i][j][z]);
				temp2+=(vx[i][j][z]*vx[i][j][z]+vy[i][j][z]*vy[i][j][z]+vz[i][j][z]*vz[i][j][z]);
			}
			temp1=sqrt(temp1);
			temp2=sqrt(temp2);
			err=temp1/(temp2);
			return err;
}
*/

void output(int m)  //���
{
	int i,j,z,Rho_wall_fake=-2;  //�������ܶ����¸�ֵ��Ϊ�����ܶ�ͼ����ʾ����

for(z=0;z<LZ;z++)
	for(j=0;j<LY;j++)
		for(i=0;i<LX;i++)//ע�� i �ķ�Χ(ֻ�����岿�ֳ�ʼ��)
		{
			if(flag[i][j][z]==solid&&z==0)
				n0[i][j][z]=Rho_wall_fake;
			if(flag[i][j][z]==solid&&z==LZ-1)
				n0[i][j][z]=n0_out ;
		}

	ostringstream name;	
	name<<"result_"<<m<<".dat";
	ofstream out(name.str().c_str());
	out<<"Title=\"Laplace Equation\"\n"<<"VARIABLES=\"X\",\"Y\",\"Z\",\"V1\",\"V2\",\"V3\",\"Rho\",\"Pressure\"\n"<<"ZONE T=\"BOX\",I="
		<<LX<<",J="<<LY<<",K="<<LZ<<",F=POINT"<<endl;

	for(z=0;z<LZ;z++)
		for(j=0;j<LY;j++)
			for(i=0;i<LX;i++)
			{
				out<<i<<"   "<<j<<"   "<<z<<"   "<<vx[i][j][z]<<"   "<<vy[i][j][z]<<"   "<<vz[i][j][z]<<"   "<<n0[i][j][z]<<"   "<<pressure[i][j][z]<<endl;
			}

}


void output_vetical_velocity(int m)
{
	/*ostringstream name;
    name<<"vetical_velocity_"<<m<<".txt";
    ofstream out(name.str().c_str());
    
        for(int i=0;i<=m;i++)
            out<<i<<"   "<<vetical_vel[i]<<"   "<<vetical_positive[i]<<"   "<<vetical_negative[i]<<"   "<<horizontal_vel[i]<<"   "<<horizontal_positive[i]<<"   "<<horizontal_negative[i]<<"   "<<f_u[i]<<"   "<<f_v[i]<<"   "<<f_vv[i]<<"   "<<A1[i]<<"   "<<A2[i]<<"   "<<bridge_Radius[i]<<endl;
    */

	ostringstream name;
	name<<"contact_angle"<<m<<".txt";
    ofstream out(name.str().c_str());
    
        for(int i=0;i<=m;i++)
            out<<i<<"   "<<contactangle[i]<<endl;

	/*ostringstream name;
	int z;
	name<<"pressure_"<<m<<".txt";
	ofstream out(name.str().c_str());	
		for(z=0;z<LZ;z++)
		{
			out<<z<<"   "<<n0[X0][Y0][z]<<"   "<<pressure[X0][Y0][z]<<endl;
		}
		*/


}

//main
void main()
{
	using namespace std;
	init_space();
	e[0][0]=0.0;  e[0][1]=0.0;  e[0][2]=0.0;
	e[1][0]=1.0;  e[1][1]=0.0;  e[1][2]=0.0;
	e[2][0]=-1.0;  e[2][1]=0.0;  e[2][2]=0.0;
	e[3][0]=0.0; e[3][1]=1.0;  e[3][2]=0.0;
	e[4][0]=0.0;  e[4][1]=-1.0;  e[4][2]=0.0;

	e[5][0]=0.0;  e[5][1]=0.0;  e[5][2]=1.0;
	e[6][0]=0.0; e[6][1]=0.0;  e[6][2]=-1.0;
	e[7][0]=1.0; e[7][1]=1.0;  e[7][2]=0.0;
	e[8][0]=-1.0;  e[8][1]=1.0;  e[8][2]=0.0;
	e[9][0]=1.0;  e[9][1]=-1.0;  e[9][2]=0.0;

	e[10][0]=-1.0;  e[10][1]=-1.0;  e[10][2]=0.0;
	e[11][0]=1.0;  e[11][1]=0.0;  e[11][2]=1.0;
	e[12][0]=-1.0;  e[12][1]=0.0;  e[12][2]=1.0;
	e[13][0]=1.0;  e[13][1]=0.0;  e[13][2]=-1.0;
	e[14][0]=-1.0;  e[14][1]=0.0;  e[14][2]=-1.0;

	e[15][0]=0.0;  e[15][1]=1.0;  e[15][2]=1.0;
	e[16][0]=0.0;  e[16][1]=-1.0;  e[16][2]=1.0;
	e[17][0]=0.0;  e[17][1]=1.0;  e[17][2]=-1.0;
	e[18][0]=0.0;  e[18][1]=-1.0;  e[18][2]=-1.0;

	time0=clock();
	ini_geo();
	initial();
	angle=180.0; //��ʼ�Ӵ��ǣ������жϼ�������нӴ��Ǳ仯��
	printf("tau=%f\n",tau0);
	printf("Start calculation, wait...\n");
	int m; //ѭ����ֵ
	//output(-1); //�����ʼʱ�̵��ܶȷֲ�
	//vetical_velocity(-1);
	//double temp1,temp2,temp3,temp4;

	for(m=0;m<=15000 ;m++)
	{
		//temp1=n0[X0][7]; //����뿪����ʱ��
		//temp3=vetical_vel[m-1];
		evolve();
		contact_angle(m);
		//temp2=n0[X0][7];
		//vetical_velocity(m); //������ֱ�����ٶ�	
		//temp4=vetical_vel[m];
		//efficient_kinetic(m); //���㵱������
		//bridge_radius(m);  //����Һ�Ű뾶�仯
		time1=clock();
		if(m%5==0)
		{
			//error=Error();
			cout<<"The"<<m<<"th computation result:"<<endl<<"The densities 4,3,2,1 are:"
				<<setprecision(6)<<n0[X0][Y0][4]<<","<<n0[X0][Y0][3]<<","<<n0[X0][Y0][2]<<","<<n0[X0][Y0][1]<<endl;
			//cout<<"The max relative error of uv is:"<<setiosflags(ios::scientific)<<error<<endl;
			printf("�Ӵ���=%f\t\tʱ��=%f\n",contactangle[m],(time1-time0)/1000.0);
			//printf("��ֱ�����ٶ�=%f\t\n",vetical_vel[m]);
			//printf("f_u=%f\t\tf_v=%f\t\tf_vv=%f%%\t\n",f_u[m],f_v[m],f_vv[m]*100);
			//printf("a1=%f\t\ta2=%f\t\tbridge radius=%f\t\n",A1[m],A2[m],bridge_Radius[m]);

			//contact_angle(m); //����Ӵ���
			//output_contant_angle(m);
			
		}

		if(m%1000==0)
		{
			
			output(m); //�������
			output_vetical_velocity(m);

		}
		/*	
		if(m<=2000&&m%10==0)
		{
			output_vetical_velocity(m);
			output(m); //�������

		}

		if(m>2000)
		{
			if(temp3>=0&&temp4<0) {output_vetical_velocity(m);output(m); break;}
			if(m%1000==0) {output_vetical_velocity(m);output(m); } //output_contant_angle(m);	
			if(error<ERR) {output(m); break;}
		 }

		if(temp1>=(n0_in+n0_out)/2&&temp2<(n0_in+n0_out)/2)
		{
			output_vetical_velocity(m);
			output(m);
		}
		*/
		
		
	}

	

    //printf("�������error=%e\t����ѭ������IT=%d\n",error,m);

  system("pause");

}
