//����ָ����ʽ��Ч�ܶȣ�EMD������Gong�����������Ӽ�������
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

using namespace std;


#define R                                  25.0//initial(��ʼҺ�ΰ뾶)
#define X0                                 50 //Һ�γ�ʼλ��
#define Y0                                 50
#define Z0                                 26//

#define Beta                               0.886//����SC����ʱ���������Ӽ�������������ϲ���
 
#define tau0                               1.0 //��ԥʱ��
#define G                                  -0.3175 //���������������ҪС���ٽ�G=-0.222����Ӧ���¶�=-1/G
#define GW0                                -0.13 //���������������ǿ��



#define LX                                101
#define LY                                101
#define LZ                                101
#define Q                                 19

#define fluid                            0  //initial_geo
#define solid                            1
#define PI                               3.1415926


#define n0_in                              2.509  //Һ������ٽ��ܶ�ln2=0.693
#define n0_out                             0.0824 ////�����ٽ��ܶ�ln2=0.693



int flag[LX][LY][LZ];
double s[LX][LY][LZ], e[Q][3],
        f0[LX][LY][LZ][Q],n0[LX][LY][LZ],
		veq0x[LX][LY][LZ],veq0y[LX][LY][LZ],veq0z[LX][LY][LZ],vx[LX][LY][LZ],vy[LX][LY][LZ],vz[LX][LY][LZ],
		u0x[LX][LY][LZ],u0y[LX][LY][LZ],u0z[LX][LY][LZ],vx0[LX][LY][LZ],vy0[LX][LY][LZ],vz0[LX][LY][LZ],
		vetical_vel[100000],vetical_positive[100000],vetical_negative[100000],horizontal_vel[100000],horizontal_positive[100000],horizontal_negative[100000],
		A1[100000],A2[100000],bridge_Radius[100000],
		f_u[100000],f_v[100000],f_vv[100000],
		g0[LX][LY][LZ][Q],F0x[LX][LY][LZ],F0y[LX][LY][LZ],F0z[LX][LY][LZ],error,
	    pressure[LX][LY][LZ],p[LX][LY][LZ],
		C[LX][LY][LZ];
double gravity0,angle;
double w[Q]={1.0/3,1.0/18,1.0/18,1.0/18,1.0/18,1.0/18,1.0/18,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36}; 
int r[Q]={0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17}; 


double eq(int k,double n,double ux,double uy,double uz);
void ini_geo(void);
void initial(void);
void collision(void);
void stream(void);
void para(void);
void contact_angle(void);
void evolve(void);
double Error();
void output(int m);



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
	int i,j,k,z;

 for(z=0;z<LZ;z++)
	for(j=0;j<LY;j++)
		for(i=0;i<LX;i++)//ע�� i �ķ�Χ(ֻ�����岿�ֳ�ʼ��)			
		{
						
				if((double((i-X0)*(i-X0)+(j-Y0)*(j-Y0)+(z-Z0)*(z-Z0)))<=R*R)
				{
					n0[i][j][z]=n0_in;
					
				}
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
		para();

}
//collision
void collision()
{
	int i,j,k,z;
    double temp0,temp1;

#pragma omp parallel for
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
		

}
//stream
void stream()
{
	int i,j,z,k,id,jd,zd;
#pragma omp parallel for
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
					n0[i][j][z]=temp0; 
					psi0[i][j][z]=1-exp(-n0[i][j][z]);	//��Ч�ܶȣ�SC-EOS
					pressure[i][j][z]=n0[i][j][z]/3+0.5*G*6*(1-exp(-n0[i][j][z]))*(1-exp(-n0[i][j][z])); //ѹ��
					}
					else
					{
						n0[i][j][z]=0.0;
						psi0[i][j][z]=0.0;
						pressure[i][j][z]=0.0;
					}
				}		
 

	

	int id,jd,zd;
	double f0x,f0y,f0z,g,gw0;

	for(z=1;z<LZ-1;z++)
		for(j=0;j<LY;j++)
			for(i=0;i<LX;i++)
		{
			f0x=f0y=f0z=0.0;
			if(flag[i][j][z]==fluid)
			{
				for(k=1;k<Q;k++)
				{
					id=i+int(e[k][0]);jd=j+int(e[k][1]);zd=z+int(e[k][2]);
					if(id>LX-1) id=0; if(id<0) id=LX-1;
					if(jd>LY-1) jd=0; if(jd<0) jd=LY-1;					

					if(k>=1&&k<=6) 	 
					{g=G;gw0=GW0; }
					else if(k>=7&&k<=18) 
					{g=G/2.0; gw0=GW0/2.0; }
							
		
					f0x+=-Beta*psi0[i][j][z]*psi0[id][jd][zd]*g*e[k][0]-(1-Beta)/2*psi0[id][jd][zd]*psi0[id][jd][zd]*g*e[k][0]-psi0[i][j][z]*s[id][jd][zd]*gw0*e[k][0];
					f0y+=-Beta*psi0[i][j][z]*psi0[id][jd][zd]*g*e[k][1]-(1-Beta)/2*psi0[id][jd][zd]*psi0[id][jd][zd]*g*e[k][1]-psi0[i][j][z]*s[id][jd][zd]*gw0*e[k][1];
					f0z+=-Beta*psi0[i][j][z]*psi0[id][jd][zd]*g*e[k][2]-(1-Beta)/2*psi0[id][jd][zd]*psi0[id][jd][zd]*g*e[k][2]-psi0[i][j][z]*s[id][jd][zd]*gw0*e[k][2];
			   }
			F0x[i][j][z]=f0x;F0y[i][j][z]=f0y;F0z[i][j][z]=f0z; 
			}
		}



double mome0x,mome0y,mome0z;


	for(z=0;z<LZ;z++)
		for(j=0;j<LY;j++)
			for(i=0;i<LX;i++)
		{
			if(flag[i][j][z]==fluid)
			{
				vx0[i][j][z]=vx[i][j][z]; //������һʱ����ʵ�ٶȣ����������ж�
				vy0[i][j][z]=vy[i][j][z];
				vz0[i][j][z]=vz[i][j][z];
				mome0x=mome0y=mome0z=0.0;
				for(k=0;k<Q;k++)
				{
					mome0x+=f0[i][j][z][k]*e[k][0]; //x������
					mome0y+=f0[i][j][z][k]*e[k][1]; //y������
					mome0z+=f0[i][j][z][k]*e[k][2]; //z������
					
				}
				mome0x=mome0x;mome0y=mome0y;mome0z=mome0z;				

			
				veq0x[i][j][z]=(mome0x/n0[i][j][z]); // EDM���������ƽ��̬�ٶ�
				veq0y[i][j][z]=(mome0y/n0[i][j][z]);
				veq0z[i][j][z]=(mome0z/n0[i][j][z]);

			
				u0x[i][j][z]=(veq0x[i][j][z]+F0x[i][j][z]/n0[i][j][z]); //EDM����������
				u0y[i][j][z]=(veq0y[i][j][z]+F0y[i][j][z]/n0[i][j][z]);
				u0z[i][j][z]=(veq0z[i][j][z]+F0z[i][j][z]/n0[i][j][z]);

				
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


void contact_angle()  //����Ӵ���
{

	double temp1,temp2;
	int i,z;
	double b0,b1,b2,a0,a1,Rho_interface=(n0_in+n0_out)/2,  //���ܶȵ��ڣ�Rho_Һ+Rho_����/2��Ϊ���棬�˴���j=2ʱ���Ŵ���R_interface�ĵ�
		Radius;
	double angle_former=angle,ratio; //��¼��һ�νӴ��ǣ������жϱ仯

	for(z=1;z<LZ-2;z++)
	{
		temp1=n0[X0][Y0][z];
		temp2=n0[X0][Y0][z+1];

		if(temp1>Rho_interface&&temp2<=Rho_interface)
			a1=double(z)+(Rho_interface-temp1)/(temp2-temp1);
	}
	a0=a1-3.0;  //��ڸ߶�

	for(i=0;i<LX-1;i++)
	{
		temp1=n0[i][Y0][3];
		temp2=n0[i+1][Y0][3];

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

	ratio=fabs((angle-angle_former)/angle_former);  //�Ӵ��Ǳ仯��

	printf("�뾶=%f\t\t�Ӵ���=%f\t\t�仯��=%f%%\t\n",Radius,angle,ratio*100);

}

//evolve
void evolve()
{
	collision();
	stream();
	para();
	
}


double Error()
{
	int i,j,z;
	double temp1,temp2,err;
	temp1=0;
	temp2=0;

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


void output(int m)  //���
{
	/*
	int i,j,z,Rho_wall_fake=-2;  //�������ܶ����¸�ֵ��Ϊ�����ܶ�ͼ����ʾ

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
	*/
	int i, j, z;
	ostringstream name;
	name << "cavity_" << m << ".data";
	ofstream out(name.str().c_str(), ofstream::binary);
	for (j = 0; j<LY; j++)
		for (i = 0; i<LX; i++){
			out.write((char*)&n0[i][j][26], sizeof(double));
		}

}



//main
void main()
{
	using namespace std;
	e[0][0]=0.0;  e[0][1]=0.0;  e[0][2]=0.0;
	e[1][0]=1.0;  e[1][1]=0.0;  e[1][2]=0.0;
	e[2][0]=-1.0;  e[2][1]=0.0;  e[2][2]=0.0;
	e[3][0]=0.0; e[3][1]=1.0;  e[3][2]=0.0;
	e[4][0]=0.0;  e[4][1]=-1.0;  e[4][2]=0.0;
	e[5][0]=0.0;  e[5][1]=0.0;  e[5][2]=1.0;
	e[6][0]=0.0; e[6][1]=0.0;  e[6][2]=-1.0;
	e[7][0]=1.0; e[7][1]=1.0;  e[7][2]=0.0;
	e[8][0]=-1.0;  e[8][1]=-1.0;  e[8][2]=0.0;
	e[9][0]=1.0;  e[9][1]=-1.0;  e[9][2]=0.0;
	e[10][0]=-1.0;  e[10][1]=1.0;  e[10][2]=0.0;
	e[11][0]=1.0;  e[11][1]=0.0;  e[11][2]=1.0;
	e[12][0]=-1.0;  e[12][1]=0.0;  e[12][2]=-1.0;
	e[13][0]=1.0;  e[13][1]=0.0;  e[13][2]=-1.0;
	e[14][0]=-1.0;  e[14][1]=0.0;  e[14][2]=1.0;
	e[15][0]=0.0;  e[15][1]=1.0;  e[15][2]=1.0;
	e[16][0]=0.0;  e[16][1]=-1.0;  e[16][2]=-1.0;
	e[17][0]=0.0;  e[17][1]=1.0;  e[17][2]=-1.0;
	e[18][0]=0.0;  e[18][1]=-1.0;  e[18][2]=1.0;

	ini_geo();
	initial();
	angle=180.0; //��ʼ�Ӵ��ǣ������жϼ�������нӴ��Ǳ仯��
	printf("Start calculation, wait...\n");
	int m; //ѭ����ֵ
	output(-1); //�����ʼʱ�̵��ܶȷֲ�

	for(m=0; ;m++)
	{

		evolve();

		if(m%100==0)
		{
			
			cout<<"The"<<m<<"th computation result:"<<endl;

			contact_angle(); //����Ӵ���			
			
		}

		//if(m%4000==0)
		{
			
			output(m); //�������
			//output_vetical_velocity(m);

		}

		
	}

	   

  system("pause");

}
