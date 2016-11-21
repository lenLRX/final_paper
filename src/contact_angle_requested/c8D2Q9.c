#include "c8D2Q9.h"

int e[Q][2]={{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};
double w[Q]={4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};
int opp[Q]={0,3,4,1,2,7,8,5,6};
double rho[NX+1][NY+1],u[NX+1][NY+1],v[NX+1][NY+1],u0[NX+1][NY+1],v0[NX+1][NY+1],f[NX+1][NY+1][Q],ft[NX+1][NY+1][Q];
double c,cs2,dx=1.0,dy=1.0,Lx,Ly,dt=1.0,rho0=1.0,tau_f, niu,error,Cd,Cl;
int none=0,fluid=1,inlet=2,outlet=3,wall1=4,wall2=5;     //格点标示
int flag[NX+1][NY+1];
double fx[NX+1][NY+1][Q];
double stream[NX+1][NY+1],vorticity[NX+1][NY+1];   //流函数,涡函数
 
int i,j,k,ip,jp,iq,jq,id,jd;
main()
{
	long n;
	init();
	for(n=0;n<=TotNumbs;n++)
	{
		collision();
	    streaming();

	    boundary();
		caculateForce(n);
		macro();
		ERROR();
		if(n%100==0)
			{
				printf("The %ldth computation\'s max relative error of velocity is: %f\n",n,error);
				printf("拽力系数=%f\t\n",Cd);	
		}
		//if(n%10000==0) output(n);
	if(error<1.0e-6)
	{
		output(n);
        printf("step=%d \n",n);		
        //return 0;
        system("pause");
	}

	}
    return 0;

}

/**************************************************************/
void init()
{
	Lx=dx*NX;
    Ly=dy*NY;
	c=dx/dt;
	cs2=c*c/3.0;
	niu=U*Diameter/Re;

	tau_f=niu/(dt*cs2)+0.5;
	printf("tauf=%f\n",tau_f);

	for(i=0;i<=NX;i++)
		for(j=0;j<=NY;j++)
		{
			u[i][j]=0;
            v[i][j]=0;
			stream[i][j]=0;
			vorticity[i][j]=0;
			rho[i][j]=rho0;
			flag[i][j]=fluid;
		}

 /*   for(i=0;i<=NX;i++) 
	{
		flag[i][NY]=wall1;
	    flag[i][0]=wall2;
	}
*/
	for(j=0;j<=NY;j++)
	{
		flag[0][j]=inlet;
		u[0][j]=U;
		flag[NX][j]=outlet;
	}
	
	for(i=0;i<=NX;i++)
		for(j=0;j<=NY;j++)
			if((i-CenterX*1.0)*(i-CenterX*1.0)+(j-CenterY*1.0)*(j-CenterY*1.0)<=0.25*Diameter*Diameter)
			{
				flag[i][j]=none;
				rho[i][j]=0;
			}

	for(i=0;i<=NX;i++)
		for(j=0;j<=NY;j++)
			for(k=0;k<Q;k++)
				f[i][j][k]=feq(k,rho[i][j],u[i][j],v[i][j]);


}


/**************************************************************/ 
double feq(int k,double rho,double u,double v)
{
	double euv,uuvv,feq;
	euv=(e[k][0]*u+e[k][1]*v);
	uuvv=(u*u+v*v);
	feq=w[k]*rho*(1.0+3.0*euv+4.5*euv*euv-1.5*uuvv);
	return feq;
}


/**************************************************************/ 
void collision()
{
	for(i=0;i<=NX;i++)
		for(j=0;j<=NY;j++)
			if(flag[i][j]==fluid)
			{ 
			  for(k=0;k<Q;k++)
              f[i][j][k]=f[i][j][k]+(feq(k,rho[i][j],u[i][j],v[i][j])-f[i][j][k])/tau_f;
			}
}


/**************************************************************/ 
void streaming()
{	
	for(i=0;i<=NX;i++)            
		for(j=0;j<=NY;j++)
			if(flag[i][j]==fluid)
			{ 
				for(k=0;k<Q;k++)
				{
					ip=i-e[k][0];
			        jp=j-e[k][1];
		            if(ip>=0&&ip<=NX)
					{
						if(jp>NY) jp=0;
			            if(jp<0) jp=NY;
						if(flag[ip][jp]!=none) ft[i][j][k]=f[ip][jp][k];
					}
				}
			}


	for(i=0;i<=NX;i++)            
		for(j=0;j<=NY;j++)
			if(flag[i][j]!=none)
			{
			    for(k=0;k<Q;k++)
				{
				   ip=i+e[k][0];
				   jp=j+e[k][1];
				   if(ip>=0&&ip<=NX)
				   {
					  if(jp>NY) jp=0;
			          if(jp<0) jp=NY;
					  if(flag[ip][jp]==inlet||flag[ip][jp]==outlet)        //ip,jp为边界格点
					  ft[ip][jp][k]=f[i][j][k];
				   }
				}
			}


}




/**************************************************************/
void boundary()
{
	double a,b;//流体、固体格点连线与圆的交点
	double a1,a2;//连线与圆的交点两个可能得x值
	double b1,b2;//连线与圆的交点两个可能得y值
	double ijp;//m方向上两点距离
	double k1,k2,k3,k4,k5;
	double d1,d2;
	double delt;
	double ux,vx,xx;

/*
	for(i=0;i<=NX;i++)
		for(j=0;j<=NY;j++)
			if(flag[i][j]==fluid)
			{
				for(k=0;k<Q;k++)
				{
					ip=i+e[k][0];
				    jp=j+e[k][1];
					if(ip>=0&&ip<=NX&&jp>=0&&jp<=NY)
					{
						switch(flag[ip][jp])
						{
						case 4:ft[ip][jp][k]=f[i][j][k];break;
						case 5:ft[ip][jp][k]=f[i][j][k];break;
						case 2:ft[ip][jp][k]=f[i][j][k];break;
//						case 3:ft[ip][jp][k]=f[i][j][k];break;
//						case 6:ft[ip][jp][k]=f[i][j][k];ft[ip][jp][opp[k]]=ft[ip][jp][k];break;
						}
					}
				}
			}
*/
//进口，速度入口
    for(i=0;i<=NX;i++)            
		for(j=0;j<=NY;j++)
			if(flag[i][j]==inlet)
			{
				rho[i][j]=(ft[i][j][0]+ft[i][j][2]+ft[i][j][4]+2.0*(ft[i][j][3]+ft[i][j][6]+ft[i][j][7]))*c/(c-U);
                ft[i][j][1]=ft[i][j][3]+2.0*rho[i][j]*U/(3*c);
				ft[i][j][5]=ft[i][j][7]-0.5*(ft[i][j][2]-ft[i][j][4])+rho[i][j]*U/6.0;
				ft[i][j][8]=ft[i][j][6]+0.5*(ft[i][j][2]-ft[i][j][4])+rho[i][j]*U/6.0;
			}
//出口，自由出流
           
	for(j=0;j<=NY;j++)
	{			
		rho[NX][j]=rho[NX-1][j];
		u[NX][j]=u[NX-1][j];
		v[NX][j]=v[NX-1][j];
		for(k=0;k<Q;k++)
		{
			ft[NX][j][k]=ft[NX-1][j][k];			
		}
	}
	     
/*	for(i=0;i<=NX;i++) //出口，非平衡外推条件           
		for(j=0;j<=NY;j++)	
		{		
			if(flag[i][j]==outlet)
			{
				rho[i][j]=rho[i-1][j];
				u[i][j]=u[i-1][j];
				v[i][j]=v[i-1][j];

				for(k=0;k<Q;k++)
				ft[i][j][k]=feq(k,rho[i][j],u[i][j],v[i][j])+(1.0-1.0/tau_f)*(ft[i-1][j][k]-feq(k,rho[i-1][j],u[i-1][j],v[i-1][j]));
			}
		}

	for(j=0;j<=NY;j++)   //出口边界，速度边界
	{			
		//rho[NX][j]=rho[NX-1][j];
		u[NX][j]=u[NX-1][j];
		v[NX][j]=0;		
		 rho[NX][j]=(ft[NX][j][0]+ft[NX][j][2]+ft[NX][j][4]+2.0*(ft[NX][j][1]+ft[NX][j][5]+ft[NX][j][8]))/(1+u[NX][j]);
         ft[NX][j][3]=ft[NX][j][1]-2*rho[NX][j]/3*u[NX][j];
         ft[NX][j][6]=ft[NX][j][8]+1/2*(ft[NX][j][4]-ft[NX][j][2])-rho[NX][j]*u[NX][j]/6;
	     ft[NX][j][7]=ft[NX][j][5]+1/2*(ft[NX][j][2]-ft[NX][j][4])-rho[NX][j]*u[NX][j]/6;
	}
*/			
		
 

/*********FH格式***********************************/

	for(i=0;i<=NX;i++)
		for(j=0;j<=NY;j++)
			if(flag[i][j]==fluid)
			{
				for(k=0;k<Q;k++)
				{
				   ip=i+e[k][0];
				   jp=j+e[k][1];
				   iq=i-e[k][0];
				   jq=j-e[k][1];
				   id=i-e[k][0]-e[k][0];
				   jd=j-e[k][1]-e[k][1];
				   if(ip>=0&&ip<=NX&&iq>=0&&iq<=NX&&id>=0&&id<=NX)
				   {
					  if(jp>NY) jp=0;
			          if(jp<0) jp=NY;
					  if(jq>NY) jq=0;
					  if(jq<0) jq=NY;
					  if(flag[ip][jp]==none)        //ip,jp为固体格点
					  {
						  if(k==2)
						  {
							  a=i;
							  b=CenterY-sqrt(0.25*Diameter*Diameter-(i-CenterX)*(i-CenterX));
							  delt=sqrt((a-i)*(a-i)+(b-j)*(b-j));
						  }
						  else if(k==4)
						  {
							  a=i;
							  b=CenterY+sqrt(0.25*Diameter*Diameter-(i-CenterX)*(i-CenterX));
							  delt=sqrt((a-i)*(a-i)+(b-j)*(b-j));
						  }
						  else if(k==1||k==3||k==5||k==6||k==7||k==8)
						  {
						      k1=(jp-j)/(ip-i);
						      k2=-k1*ip+jp-CenterY;
						      k3=1+k1*k1;
						      k4=2.0*k1*k2-2.0*CenterX;
						      k5=CenterX*CenterX+k2*k2-0.25*Diameter*Diameter;
						      a1=(-k4+sqrt(k4*k4-4.0*k3*k5))/(2.0*k3);
						      a2=(-k4-sqrt(k4*k4-4.0*k3*k5))/(2.0*k3);
						      b1=k1*(a1-ip)+jp;
						      b2=k1*(a2-ip)+jp;
						      d1=sqrt((a1-i)*(a1-i)+(b1-j)*(b1-j))+sqrt((a1-ip)*(a1-ip)+(b1-jp)*(b1-jp));
						      d2=sqrt((a2-i)*(a2-i)+(b2-j)*(b2-j))+sqrt((a2-ip)*(a2-ip)+(b2-jp)*(b2-jp));
						      ijp=sqrt((double)((ip-i)*(ip-i))+(double)((jp-j)*(jp-j)));
						      if(d2>=d1)
							  {
							     a=a1;
							     b=b1;
							  }
						      if(d2<d1)
							  {
							     a=a2;
							     b=b2;
							  }

						      delt=sqrt((a-i)*(a-i)+(b-j)*(b-j))/ijp;
						  }
						/*  
						  if(delt>0&&delt<0.5)
						  {
							  ux=u[i][j];
							  vx=v[i][j];
							  xx=(2.0*delt-1.0)/(tau_f-1.0);
							  fx[ip][jp][k]=w[k]*rho[i][j]*(1.0+3.0*(e[k][0]*ux+e[k][1]*vx)+4.5*(e[k][0]*u[i][j]+e[k][1]*v[i][j])*(e[k][0]*u[i][j]+e[k][1]*v[i][j])-1.5*(u[i][j]*u[i][j]+v[i][j]*v[i][j]));
						      fx[ip][jp][opp[k]]=(1-xx)*f[i][j][k]+xx*fx[ip][jp][k];
						  }
						  else if(delt>=0.5&&delt<=1.0)
						  {
							  ux=u[i][j]*(delt-1.0)/delt;
							  vx=v[i][j]*(delt-1.0)/delt;
							  xx=(2.0*delt-1.0)/(tau_f);
							  fx[ip][jp][k]=w[k]*rho[i][j]*(1.0+3.0*(e[k][0]*ux+e[k][1]*vx)+4.5*(e[k][0]*u[i][j]+e[k][1]*v[i][j])*(e[k][0]*u[i][j]+e[k][1]*v[i][j])-1.5*(u[i][j]*u[i][j]+v[i][j]*v[i][j]));
						      fx[ip][jp][opp[k]]=(1-xx)*f[i][j][k]+xx*fx[ip][jp][k];
						  }

						  ft[i][j][opp[k]]=fx[ip][jp][opp[k]];
						  */
						  if(delt>0&&delt<0.5)
						  {
							   ft[i][j][opp[k]] =delt*(1+2*delt)*f[i][j][k]+(1-4*delt*delt)*f[iq][jq][k]-delt*(1-2*delt)*f[id][jd][k];
						  }
						  else if(delt>=0.5&&delt<=1.0)
						  {
							  ft[i][j][opp[k]] =f[i][j][k]/(2*delt+1)/delt+(2*delt-1)/delt*f[i][j][opp[k]]+(1-2*delt)/(1+2*delt)*f[iq][jq][opp[k]];
						  }
						  // ft[i][j][opp[k]] = ((1-delt)*f[iq][jq][k]+delt*f[i][j][k]+delt*f[i][j][opp[k]])/(1+delt); //Yu的一次线性插值

					  }
				   }
				}
			}
/**************************************************/




		
}




/**************************************************************/
void caculateForce(long m)
{
	double cl=0.0,cd=0.0;
	FILE *fp1,*fp2;

	for(i=0;i<=NX;i++)            
		for(j=0;j<=NY;j++)
			if(flag[i][j]==fluid)
			{
				for(k=0;k<Q;k++)
				{
				   ip=i+e[k][0];
				   jp=j+e[k][1];
				   if(ip>=0&&ip<=NX)
				   {
					   if(jp>NY) jp=0;
					   if(jp<0) jp=NY;
					  if(flag[ip][jp]==none)        //ip,jp为边界格点
					  {
						  cl+=(f[i][j][k]+ft[i][j][opp[k]])*e[k][1];
					      cd+=(f[i][j][k]+ft[i][j][opp[k]])*e[k][0];
					  }
				   }				
				}
			}
	
	Cl=cl/(0.5*rho0*U*U*Diameter);
    Cd=cd/(0.5*rho0*U*U*Diameter);

    if(m==0)  
	{
	    if((fp1=fopen("cl.plt","w"))==NULL)  printf("cannot open the cl.plt file\n");
	    fprintf(fp1,"variables=\"step\",\"cl\"\n");
		fprintf(fp1,"%ld%,%f\n",m,Cl);

        if((fp2=fopen("cd.plt","w"))==NULL)  printf("cannot open the cd.plt file\n");
	    fprintf(fp2,"variables=\"step\",\"cd\"\n");
		fprintf(fp2,"%ld%,%f\n",m,Cd);

	    fclose(fp1);
	    fclose(fp2);
	}
	
	
	
    else
	{
    	if((fp1=fopen("cl.plt","a"))==NULL)  printf("cannot open the cl.plt file\n");
    	fprintf(fp1,"%ld%,%f\n",m,Cl);

        if((fp2=fopen("cd.plt","a"))==NULL)  printf("cannot open the cd.plt file\n");
	    fprintf(fp2,"%ld%,%f\n",m,Cd);


	    fclose(fp1);
	    fclose(fp2);
	}


}
/**************************************************************/
void macro()
{
	
//ft值传给f
	for(i=0;i<=NX;i++)            
		for(j=0;j<=NY;j++)
			if(flag[i][j]!=none)
			{
				for(k=0;k<Q;k++)
				f[i][j][k]=ft[i][j][k];
			}
	for(i=0;i<=NX;i++)            //计算宏观量
		for(j=0;j<=NY;j++)
			if(flag[i][j]==fluid)
			{
			      u0[i][j]=u[i][j];
                  v0[i][j]=v[i][j];
			      rho[i][j]=0;
                  u[i][j]=0;
		 	      v[i][j]=0;
			      for(k=0;k<Q;k++)
				  {
				      rho[i][j]+=f[i][j][k];
                      u[i][j]+=e[k][0]*f[i][j][k];
                      v[i][j]+=e[k][1]*f[i][j][k];
				  }
			 u[i][j]/=rho[i][j];
             v[i][j]/=rho[i][j];
			}


/*
	for(i=0;i<=NX;i++)            
		for(j=0;j<=NY;j++)
			if(flag[i][j]==fluid)
			{
				for(k=0;k<=j;k++) stream[i][j]+=u[i][k];
				vorticity[i][j]=(v[i+1][j]-v[i-1][j])/2.0-(u[i][j+1]-u[i][j-1])/2.0;
			}
*/

/*
//上下边界，非平衡态外推
	for(i=0;i<=NX;i++)            
		for(j=0;j<=NY;j++)	
		{
			if(flag[i][j]==wall1)
			{
				rho[i][j]=rho[i][j-1];
				f[i][j][4]=feq(4,rho[i][j],u[i][j],v[i][j])+(1.0-1.0/tau_f)*(f[i][j-1][4]-feq(4,rho[i][j-1],u[i][j-1],v[i][j-1]));
				f[i][j][7]=feq(7,rho[i][j],u[i][j],v[i][j])+(1.0-1.0/tau_f)*(f[i][j-1][7]-feq(7,rho[i][j-1],u[i][j-1],v[i][j-1]));
				f[i][j][8]=feq(8,rho[i][j],u[i][j],v[i][j])+(1.0-1.0/tau_f)*(f[i][j-1][8]-feq(8,rho[i][j-1],u[i][j-1],v[i][j-1]));
			}
			else if(flag[i][j]==wall2)
			{
                rho[i][j]=rho[i][j+1];
				f[i][j][2]=feq(2,rho[i][j],u[i][j],v[i][j])+(1.0-1.0/tau_f)*(f[i][j+1][2]-feq(2,rho[i][j+1],u[i][j+1],v[i][j+1]));
				f[i][j][5]=feq(5,rho[i][j],u[i][j],v[i][j])+(1.0-1.0/tau_f)*(f[i][j+1][5]-feq(5,rho[i][j+1],u[i][j+1],v[i][j+1]));
				f[i][j][6]=feq(6,rho[i][j],u[i][j],v[i][j])+(1.0-1.0/tau_f)*(f[i][j+1][6]-feq(6,rho[i][j+1],u[i][j+1],v[i][j+1]));
			}
		}
	f[0][0][1]=f[0][0][3];
	f[0][0][8]=f[0][0][6];

	f[NX][0][3]=f[NX][0][1];
	f[NX][0][7]=f[NX][0][5];

    f[0][NY][1]=f[0][NY][3];
	f[0][NY][5]=f[0][NY][7];

	f[NX][NY][3]=f[NX][NY][1];
	f[NX][NY][6]=f[NX][NY][8];
*/

}
/**************************************************************/

void ERROR()
{
	double temp1,temp2;
	temp1=0.0;
	temp2=0.0;
	for(i=0;i<=NX;i++)            
		for(j=0;j<=NY;j++)
			if(flag[i][j]==fluid)
			{
			   temp1+=((u[i][j]-u0[i][j])*(u[i][j]-u0[i][j])+(v[i][j]-v0[i][j])*(v[i][j]-v0[i][j]));
			   temp2+=(u[i][j]*u[i][j]+v[i][j]*v[i][j]);
			}
		temp1=sqrt(temp1);
		temp2=sqrt(temp2);
		error=temp1/(temp2);

}

/**************************************************************/
void output(long m)   //输出
{
    char filename[100];
	FILE *fp;

	sprintf(filename,"%ld_%f.plt",m,error);
    if((fp=fopen(filename,"w"))==NULL)  printf("cannot open the %s file\n",filename);
	fputs("Title=\"LBM Unsteady Flow Past a Cylinder\"\n",fp);
    //fputs("VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"STR\",\"VOR\"\n",fp);
	fputs("VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"Rho\"\n",fp);
	fprintf(fp,"ZONE I=%d,J=%d,F=POINT\n",NX+1,NY+1);
		
	for(j=0;j<=NY;j++)		
	    for(i=0;i<=NX;i++)
		{
			//fprintf(fp,"%d,%d,%f,%f,%f,%f\n",i,j,u[i][j],v[i][j],stream[i][j],vorticity[i][j]);
			fprintf(fp,"%d,%d,%f,%f,%f\n",i,j,u[i][j],v[i][j],rho[i][j]);
		}
     fclose(fp);
}

/**************************************************************/

