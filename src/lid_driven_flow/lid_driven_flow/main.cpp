#include<iostream>
#include<cmath>
#include<cstdlib>
#include<fstream>
#include<iostream>
#include<sstream>
#include<iomanip>
#include<string>

using namespace std;
const int Q = 9;
const int NX = 256;
const int NY = 256;
const double U = 0.1;

int e[Q][2] = { { 0, 0 }, { 1, 0 }, { 0, 1 }, { -1, 0 }, { 0, -1 }, { 1, 1 }, { -1, 1 }, { -1, -1 }, { 1, -1 } };
double w[Q] = { 4.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 36, 1.0 / 36, 1.0 / 36, 1.0 / 36 };

double rho[NX + 1][NY + 1], u[NX + 1][NY + 1][2], u0[NX + 1][NY + 1][2], f[NX + 1][NY + 1][Q], F[NX + 1][NY + 1][Q];

int i, j, k, ip, jp, n;
double c, Re, dx, dy, Lx, Ly, dt, rho0, P0, tau_f, niu, error;

void init();
double feq(int k, double rho, double u[2]);
void evolution();
void output(int m);
void Error();

int main(){
	init();
	for (n = 0;; n++){
		evolution();
		if (n % 100 == 0){
			Error();
			cout << sizeof(double) << endl;
			cout << "The " << n << "th computation result:" << endl
				<< "The u,v of point (NX/2,NY/2) is: " << setprecision(6)
				<< u[NX / 2][NY / 2][0] << " , " << u[NX / 2][NY / 2][1] << endl;
			cout << "The relative error of uv is: "
				<< setiosflags(ios::scientific) << error << endl;
		}

		if (n % 1000 == 0)
			output(n);
	}

	return 0;
}

void init(){
	dx = 1.0;
	dy = 1.0;
	Lx = dx * double(NX);
	Ly = dy * double(NY);
	dt = dx;
	c = dx / dt;
	rho0 = 1.0;
	Re = 1000;
	niu = U * Lx / Re;
	tau_f = 3.0 * niu + 0.5;
	cout << "tau_f: " << tau_f << endl;

	for (i = 0; i <= NX; i++)
		for (j = 0; j <= NY; j++){
			u[i][j][0] = 0;
			u[i][j][1] = 0;
			rho[i][j] = rho0;
			u[i][NY][0] = U;
			for (k = 0; k < Q; k++)
				f[i][j][k] = feq(k, rho[i][j], u[i][j]);
		}
}

double feq(int k, double rho, double u[2]){
	double eu, uv, feq;
	eu = e[k][0] * u[0] + e[k][1] * u[1];
	uv = u[0] * u[0] + u[1] * u[1];
	feq = w[k] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * uv);
	return feq;
}

void evolution(){
#pragma omp parallel for
	for (i = 1; i < NX; i++)
		for (j = 1; j < NY; j++)
			for (k = 0; k < Q; k++){
				ip = i - e[k][0];
				jp = j - e[k][1];
				F[i][j][k] = f[ip][jp][k] + (feq(k, rho[ip][jp], u[ip][jp]) - f[ip][jp][k]) / tau_f;
			}
#pragma omp parallel for
	for (i = 1; i < NX; i++)
		for (j = 1; j < NY; j++){
			u0[i][j][0] = u[i][j][0];
			u0[i][j][1] = u[i][j][1];

			rho[i][j] = 0;
			u[i][j][0] = 0;
			u[i][j][1] = 0;

			for (k = 0; k < Q; k++){
				f[i][j][k] = F[i][j][k];
				rho[i][j] += f[i][j][k];

				u[i][j][0] += e[k][0] * f[i][j][k];
				u[i][j][1] += e[k][1] * f[i][j][k];
			}

			u[i][j][0] /= rho[i][j];
			u[i][j][1] /= rho[i][j];
		}

#pragma omp parallel for
	for (j = 1; j < NY; j++)
		for (k = 0; k < Q; k++){
			rho[NX][j] = rho[NX - 1][j];
			f[NX][j][k] = feq(k, rho[NX][j], u[NX][j]) + 
				f[NX - 1][j][k] - feq(k, rho[NX - 1][j], u[NX - 1][j]);
			rho[0][j] = rho[1][j];
			f[0][j][k] = feq(k, rho[0][j], u[0][j]) + f[1][j][k] 
				- feq(k, rho[1][j], u[1][j]);
		}
#pragma omp parallel for
	for (i = 1; i < NX; i++)
		for (k = 0; k < Q; k++){
			rho[i][0] = rho[i][1];
			f[i][0][k] = feq(k,rho[i][0],u[i][0]) 
				+ f[i][1][k] - feq(k,rho[i][1],u[i][1]);

			rho[i][NY] = rho[i][NY - 1];
			u[i][NY][0] = U;
			f[i][NY][k] = feq(k, rho[i][NY], u[i][NY]) 
				+ f[i][NY - 1][k] - feq(k, rho[i][NY - 1], u[i][NY - 1]);

		}
}

void output(int m)
{
	ostringstream name;
	name << "cavity_" << m << ".data";
	ofstream out(name.str().c_str(),ofstream::binary);

	//out << "Title = \"LBM Lid Driven Flow\"\n" << endl;

	for (j = 0; j <= NY; j++)
	for (i = 0; i <= NX; i++)
		{
			/*
			double ipos = double(i) / Lx;
			double jpos = double(i) / Lx;
			out.write((char*)&ipos, sizeof(double));
			out.write((char*)&jpos, sizeof(double));
			*/
			out.write((char*)&u[i][j][0],sizeof(double));
			out.write((char*)&u[i][j][1], sizeof(double));
			/*
			out << u[i][j][0]
				<< u[i][j][1] << endl;
				*/
		}
}

void Error(){
	double temp1, temp2;
	temp1 = 0;
	temp2 = 0;

	for (i = 1; i < NX; i++)
		for (j = 1; j < NY; j++){
			temp1 += ((u[i][j][0] - u0[i][j][0]) * (u[i][j][0] - u0[i][j][0]) + (u[i][j][1] - u0[i][j][1]) * (u[i][j][1] - u0[i][j][1]));
			temp2 += (u[i][j][0] * u[i][j][0] + u[i][j][1] * u[i][j][1]);
		}


	temp1 = sqrt(temp1);
	temp2 = sqrt(temp2);

	cout << "temp1: " << temp1 << endl;

	cout << "temp2: " << temp2 << endl;

	error = temp1 / (temp2 + 1e-30);
}