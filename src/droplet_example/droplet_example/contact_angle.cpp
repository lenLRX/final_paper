//采用指数形式有效密度，EMD添加力项，Gong方法计算粒子间作用力
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


#define R                                  25.0//initial(初始液滴半径)
#define X0                                 50 //液滴初始位置
#define Y0                                 50
#define Z0                                 26//

#define Beta                               0.886//采用SC方程时，计算粒子间作用力方程耦合参数
 
#define tau0                               1.0 //弛豫时间
#define G                                  -0.3175 //相间作用力参数，要小于临界G=-0.222，对应的温度=-1/G
#define GW0                                -0.13 //流体与壁面作用力强度



#define LX                                101
#define LY                                101
#define LZ                                101
#define Q                                 19

#define fluid                            0  //initial_geo
#define solid                            1
#define PI                               3.1415926


#define n0_in                              2.509  //液相高于临界密度ln2=0.693
#define n0_out                             0.0824 ////低于临界密度ln2=0.693



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
			
				flag[i][j][0]=solid; //下壁面
				s[i][j][0]=1.0;
				flag[i][j][LZ-1]=solid;//上壁面
				s[i][j][LZ-1]=1.0;
			
		}   


for(z=1;z<LZ-1;z++)  //标记流体
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
		for(i=0;i<LX;i++)//注意 i 的范围(只给流体部分初始化)			
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
			    vx[i][j][z]=0.0; //用于error判断
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
					temp1=eq(k,n0[i][j][z],u0x[i][j][z],u0y[i][j][z],u0z[i][j][z]); //EDM处理外力，由外力造成的密度分布函数差异

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
					psi0[i][j][z]=1-exp(-n0[i][j][z]);	//有效密度，SC-EOS
					pressure[i][j][z]=n0[i][j][z]/3+0.5*G*6*(1-exp(-n0[i][j][z]))*(1-exp(-n0[i][j][z])); //压力
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
				vx0[i][j][z]=vx[i][j][z]; //保存上一时刻真实速度，用于收敛判断
				vy0[i][j][z]=vy[i][j][z];
				vz0[i][j][z]=vz[i][j][z];
				mome0x=mome0y=mome0z=0.0;
				for(k=0;k<Q;k++)
				{
					mome0x+=f0[i][j][z][k]*e[k][0]; //x方向动量
					mome0y+=f0[i][j][z][k]*e[k][1]; //y方向动量
					mome0z+=f0[i][j][z][k]*e[k][2]; //z方向动量
					
				}
				mome0x=mome0x;mome0y=mome0y;mome0z=mome0z;				

			
				veq0x[i][j][z]=(mome0x/n0[i][j][z]); // EDM处理外力项，平衡态速度
				veq0y[i][j][z]=(mome0y/n0[i][j][z]);
				veq0z[i][j][z]=(mome0z/n0[i][j][z]);

			
				u0x[i][j][z]=(veq0x[i][j][z]+F0x[i][j][z]/n0[i][j][z]); //EDM处理外力项
				u0y[i][j][z]=(veq0y[i][j][z]+F0y[i][j][z]/n0[i][j][z]);
				u0z[i][j][z]=(veq0z[i][j][z]+F0z[i][j][z]/n0[i][j][z]);

				
				vx[i][j][z]=(veq0x[i][j][z]+0.5*F0x[i][j][z]/n0[i][j][z]); //EDM处理外力项，真实速度
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


void contact_angle()  //计算接触角
{

	double temp1,temp2;
	int i,z;
	double b0,b1,b2,a0,a1,Rho_interface=(n0_in+n0_out)/2,  //将密度等于（Rho_液+Rho_气）/2定为界面，此处当j=2时，才存在R_interface的点
		Radius;
	double angle_former=angle,ratio; //记录上一次接触角，用于判断变化

	for(z=1;z<LZ-2;z++)
	{
		temp1=n0[X0][Y0][z];
		temp2=n0[X0][Y0][z+1];

		if(temp1>Rho_interface&&temp2<=Rho_interface)
			a1=double(z)+(Rho_interface-temp1)/(temp2-temp1);
	}
	a0=a1-3.0;  //球冠高度

	for(i=0;i<LX-1;i++)
	{
		temp1=n0[i][Y0][3];
		temp2=n0[i+1][Y0][3];

		if(temp1<Rho_interface&&temp2>=Rho_interface)
			b1=double(i)+(Rho_interface-temp1)/(temp2-temp1);
		if(temp1>Rho_interface&&temp2<=Rho_interface)
			b2=double(i)+(Rho_interface-temp1)/(temp2-temp1);
	}
	b0=b2-b1; //被y=Y0截掉长度


	Radius=a0/2.0+b0*b0/(8.0*a0);

	angle=atan(b0/(2.0*Radius-2.0*a0));

	if (angle<0)
		angle=PI+angle;

	angle=angle/PI*180; //转为角度

	ratio=fabs((angle-angle_former)/angle_former);  //接触角变化率

	printf("半径=%f\t\t接触角=%f\t\t变化率=%f%%\t\n",Radius,angle,ratio*100);

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


void output(int m)  //输出
{
	/*
	int i,j,z,Rho_wall_fake=-2;  //将壁面密度重新赋值，为了在密度图中显示

	for(z=0;z<LZ;z++)
	for(j=0;j<LY;j++)
	for(i=0;i<LX;i++)//注意 i 的范围(只给流体部分初始化)
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
	angle=180.0; //初始接触角，用于判断计算过程中接触角变化率
	printf("Start calculation, wait...\n");
	int m; //循环赋值
	output(-1); //输出初始时刻的密度分布

	for(m=0; ;m++)
	{

		evolve();

		if(m%100==0)
		{
			
			cout<<"The"<<m<<"th computation result:"<<endl;

			contact_angle(); //计算接触角			
			
		}

		//if(m%4000==0)
		{
			
			output(m); //输出数据
			//output_vetical_velocity(m);

		}

		
	}

	   

  system("pause");

}
