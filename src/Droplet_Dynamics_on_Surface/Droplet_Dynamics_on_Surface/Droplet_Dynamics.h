#ifndef __DROPLET_DYNAMICS_H__
#define __DROPLET_DYNAMICS_H__
#include "Lattice4D.h"
#include "Lattice3D.h"
#include "Lattice2D.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

int e[19][3];
const int r[19] = { 0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17 };

//rho0 = 1
//Tc = 4.5rho0 = 4.5
//T = Tc x 0.65
//T = 2.925
//g = - 1 / T

//TOD0: add gravity

class Droplet_Dynamics
{
public:
	Droplet_Dynamics(int xdim = 80, int ydim = 40, int zdim = 30)
		:xdim(xdim), ydim(ydim), zdim(zdim),
		U(xdim, ydim, zdim), U0(xdim, ydim, zdim), UF(xdim, ydim, zdim),
		V(xdim, ydim, zdim), V0(xdim, ydim, zdim), VF(xdim, ydim, zdim),
		W(xdim, ydim, zdim), W0(xdim, ydim, zdim), WF(xdim, ydim, zdim),
		rho(xdim, ydim, zdim), Force(xdim,ydim,zdim,3),F(xdim, ydim, zdim, 19),
		F0(xdim, ydim, zdim, 19), FC(xdim, ydim, zdim, 19),gs(xdim, ydim){
		w[0] = 1.0 / 3.0;
		fill(w + 1, w + 7, 1.0 / 18.0);
		fill(w + 7, w + 19, 1.0 / 36.0);

		e[0][0] = 0.0;  e[0][1] = 0.0;  e[0][2] = 0.0;//
		e[1][0] = 1.0;  e[1][1] = 0.0;  e[1][2] = 0.0;//
		e[2][0] = -1.0;  e[2][1] = 0.0;  e[2][2] = 0.0;//
		e[3][0] = 0.0; e[3][1] = 1.0;  e[3][2] = 0.0;//
		e[4][0] = 0.0;  e[4][1] = -1.0;  e[4][2] = 0.0;//
		e[5][0] = 0.0;  e[5][1] = 0.0;  e[5][2] = 1.0;//
		e[6][0] = 0.0; e[6][1] = 0.0;  e[6][2] = -1.0;//
		e[7][0] = 1.0; e[7][1] = 1.0;  e[7][2] = 0.0;//
		e[8][0] = -1.0;  e[8][1] = -1.0;  e[8][2] = 0.0;//
		e[9][0] = 1.0;  e[9][1] = -1.0;  e[9][2] = 0.0;//
		e[10][0] = -1.0;  e[10][1] = 1.0;  e[10][2] = 0.0;//
		e[11][0] = 1.0;  e[11][1] = 0.0;  e[11][2] = 1.0;//
		e[12][0] = -1.0;  e[12][1] = 0.0;  e[12][2] = -1.0;//
		e[13][0] = 1.0;  e[13][1] = 0.0;  e[13][2] = -1.0;//
		e[14][0] = -1.0;  e[14][1] = 0.0;  e[14][2] = 1.0;//
		e[15][0] = 0.0;  e[15][1] = 1.0;  e[15][2] = 1.0;//
		e[16][0] = 0.0;  e[16][1] = -1.0;  e[16][2] = -1.0;//
		e[17][0] = 0.0;  e[17][1] = 1.0;  e[17][2] = -1.0;//
		e[18][0] = 0.0;  e[18][1] = -1.0;  e[18][2] = 1.0;//

		double T = 2.925;
		//g = -1 / T;
		g = -0.33898;
		g1 = g;
		g2 = g / 2;
		Fint_const = -(1 - beta) / 2;
		Fint_const_g1 = Fint_const * g1;
		Fint_const_g2 = Fint_const * g2;

		xdim_1 = xdim - 1;
		ydim_1 = ydim - 1;
		zdim_1 = zdim - 1;

		init();
	}

	void init(){
		//init solids;
#pragma omp parallel for
		for (int y = 0; y < ydim; y++){
			for (int x = 0; x < xdim;x++){
				if (x < (xdim - 12) / 2 || x >(xdim + 12))
					gs.at(x, y) = -0.25;
				else
					gs.at(x, y) = -0.1;
			}
		}
		//init droplet
#pragma omp parallel for
		for (int x = 0; x < xdim; x++){
			for (int y = 0; y < ydim; y++){
				for (int z = 0; z < zdim; z++){
					if (sqrt((x - xdim / 2) * (x - xdim / 2)
						+ (y - ydim / 2) * (y - ydim / 2)
						+ (z - zdim / 2) * (z - zdim / 2)) <= 12){
						rho.at(x, y, z) = 2.745;
					}
					else
						rho.at(x, y, z) = 0.0617;
				}
			}
		}

		//init F
#pragma omp parallel for
		for (int x = 0; x < xdim; x++)
			for (int y = 0; y < ydim; y++)
				for (int z = 0; z < zdim; z++)
					for (int k = 0; k < Q; k++)
						F.at(x, y, z, k) = feq(k, x, y, z);

	}

	inline double effective_mass(double _rho){
		return sqrt(2 * (calc_pressure(_rho) - _rho / 3) / (c0 * g));
		//return 1 - exp(-_rho);
	}

	inline double calc_pressure(double _rho){
		double _1_m_exp = 1 - exp(-_rho);
		return 1.0 / 3.0 *_rho + 0.5 * g * c0 * _1_m_exp * _1_m_exp;
	}

	inline double Fint(int k,int x, int y, int z){
		int xp = x + e[k][0];
		int yp = y + e[k][1];
		int zp = z + e[k][2];

		if (zp < 0)
			return 0;
		if (zp > zdim_1)
			return 0;

		if (xp < 0)
			xp = xdim_1;
		if (xp > xdim_1)
			xp = 0;
		if (yp < 0)
			yp = ydim_1;
		if (yp > ydim_1)
			yp = 0;

		double temp1;
		double temp2;
		double m_x = effective_mass(rho.at(xp, yp, zp));
		if(k < 7){
			temp1 = g1 * m_x;
			temp2 = Fint_const_g1 * m_x * m_x;
		}else{
			temp1 = g2 * m_x;
			temp2 = Fint_const_g2 * m_x * m_x;
		}

		temp1 = temp1 * neg_beta * effective_mass(rho.at(x,y,z));

		return temp1 + temp2;
	}

	inline double Fg(double _rho){
		return 9.8 * _rho;
	}

	inline double Fs(int x,int y,int z){
		static double w_const = w[6] + w[12] + w[13] + w[16] + w[17];
		return -effective_mass(rho.at(x, y, z))
			* w_const * gs.at(x, y);
	}

	void update(){
		
#pragma omp parallel for
		for (int x = 0; x < xdim; x++)
			for (int y = 0; y < ydim; y++)
				for (int z = 0; z < zdim; z++){
					for (int k = 0; k < Q; k++){
						double temp1 = feq(k, x, y, z);
						/*
						FC.at(x, y, z, k) = F.at(x, y, z, k)
							+ (temp1 - F.at(x, y, z, k)) / tau
						
							+ feq(k, x, y, z,
							UF.at(x, y, z),
							VF.at(x, y, z),
							WF.at(x, y, z))
							- temp1;
							*/
						FC.at(x, y, z, k) = feq(k, x, y, z,
							UF.at(x, y, z),
							VF.at(x, y, z),
							WF.at(x, y, z));
							
					}
				}

#pragma omp parallel for
		for (int x = 0; x < xdim; x++)
			for (int y = 0; y < ydim; y++)
				for (int z = 0; z < zdim; z++){

					for (int k = 0; k < Q; k++){
						int xp = x - e[k][0];
						int yp = y - e[k][1];
						int zp = z - e[k][2];

						if (zp < 0 || zp > zdim_1){
							F.at(x, y, z, k) = FC.at(x, y, z, r[k]);
						}
						else{
							if (xp < 0)
								xp = xdim_1;
							if (xp > xdim_1)
								xp = 0;

							if (yp < 0)
								yp = ydim_1;
							if (yp > ydim_1)
								yp = 0;

							F.at(x, y, z, k) = FC.at(xp, yp, zp, k);
						}

						
					}
				}
		rho.clean();
#pragma omp parallel for
		for (int x = 0; x < xdim; x++)
			for (int y = 0; y < ydim; y++)
				for (int z = 0; z < zdim; z++){
					for (int k = 0; k < Q; k++){
						rho.at(x, y, z) += F.at(x, y, z, k);
					}

					if (rho.at(x, y, z) < 0){
						cout << "x: " << x << " y: " << y << " z:" << z << endl;
						cout << "rho: " << rho.at(x, y, z) << endl;
						cout << "U: " << U.at(x, y, z)
							<< " V: " << V.at(x, y, z)
							<< " W: " << W.at(x, y, z) << endl;
						cout << "UF: " << UF.at(x, y, z)
							<< " VF: " << VF.at(x, y, z)
							<< " WF: " << WF.at(x, y, z) << endl;

						cout << "Force: " << Force.at(x, y, z, 0)
							<< " Force: " << Force.at(x, y, z, 1)
							<< " Force: " << Force.at(x, y, z, 2) << endl;
						for (int k = 0; k < Q; k++)
							cout << "k: " << k << " F: " << F.at(x, y, z, k) << endl;
						cout << "error" << endl;
					}
				}
					

		
		double max = -1;
		int max_x, max_y, max_z;
		double max_rho;
		Force.clean();
#pragma omp parallel for
		for (int x = 0; x < xdim; x++)
			for (int y = 0; y < ydim; y++)
				for (int z = 0; z < zdim; z++){
					for (int k = 1; k < Q; k++){

						int xp = x + e[k][0];
						int yp = y + e[k][1];
						int zp = z + e[k][2];

						if (zp < 0)
							continue;
						if (zp > zdim_1)
							continue;

						if (xp < 0)
							xp = xdim_1;
						if (xp > xdim_1)
							xp = 0;
						if (yp < 0)
							yp = ydim_1;
						if (yp > ydim_1)
							yp = 0;

						double temp1;
						double temp2;
						double m_x = effective_mass(rho.at(xp, yp, zp));
						if (k < 7){
							temp1 = g1 * m_x;
							temp2 = Fint_const_g1 * m_x * m_x;
						}
						else{
							temp1 = g2 * m_x;
							temp2 = Fint_const_g2 * m_x * m_x;
						}

						temp1 = temp1 * neg_beta * effective_mass(rho.at(x, y, z));

						Force.at(x, y, z, 0) += e[k][0] * (temp1 + temp2);
						Force.at(x, y, z, 1) += e[k][1] * (temp1 + temp2);
						Force.at(x, y, z, 2) += e[k][2] * (temp1 + temp2);

					}
					if (Force.at(x, y, z, 0) > max){
						max = Force.at(x, y, z, 0);
						max_x = x;
						max_y = y;
						max_z = z;
						max_rho = rho.at(x,y,z);
					}
						
				}
		/*
		cout << "max force: " << max << endl;
		cout << " rho: "<< max_rho <<" x: " << max_x << " y: " << max_y << " z: " << max_z << endl;
		cout << "rho: " << rho.at(max_x, max_y, max_z) << endl;
		cout << "U: " << U.at(max_x, max_y, max_z)
			<< " V: " << V.at(max_x, max_y, max_z)
			<< " W: " << W.at(max_x, max_y, max_z) << endl;
		cout << "UF: " << UF.at(max_x, max_y, max_z)
			<< " VF: " << VF.at(max_x, max_y, max_z)
			<< " WF: " << WF.at(max_x, max_y, max_z) << endl;

		cout << "Force: " << Force.at(max_x, max_y, max_z, 0)
			<< " Force: " << Force.at(max_x, max_y, max_z, 1)
			<< " Force: " << Force.at(max_x, max_y, max_z, 2) << endl;
		for (int k = 0; k < Q; k++)
			cout << "k: " << k << " F: " << F.at(max_x, max_y, max_z, k) << endl;
			*/
		/*
		//bottom boundary
#pragma omp parallel for
		for (int x = 0; x < xdim; x++)
			for (int y = 0; y < ydim; y++){
				int z = 0;
				Force.at(x, y, z, 2) += Fs(x, y, z); 
			}
		//bottom boundary ends
		*/


		
		U.swap(U0);
		V.swap(V0);
		W.swap(W0);
		U.clean();
		V.clean();
		W.clean();
		UF.clean();
		VF.clean();
		WF.clean();
		
#pragma omp parallel for
		for (int x = 0; x < xdim; x++)
			for (int y = 0; y < ydim; y++)
				for (int z = 0; z < zdim; z++){
					double temp_u = 0;
					double temp_v = 0;
					double temp_w = 0;
					if (x == 52 && y == 20 && z == 15)
						cout << "pause" << endl;
					for (int k = 0; k < Q; k++){
						temp_u += e[k][0] * F.at(x, y, z, k);
						temp_v += e[k][1] * F.at(x, y, z, k);
						temp_w += e[k][2] * F.at(x, y, z, k);
					}

					

					double _rho = rho.at(x, y, z);
					
					
						
					U.at(x, y, z) = temp_u / _rho;
					V.at(x, y, z) = temp_v / _rho;
					W.at(x, y, z) = temp_w / _rho;

					
					if (false && abs(U.at(x, y, z)) > 10){
						cout << "x: " << x << " y: " << y << " z:" << z << endl;
						cout << "rho: " << rho.at(x, y, z) << endl;
						cout << "U: " << U.at(x, y, z)
							<< " V: " << V.at(x, y, z)
							<< " W: " << W.at(x, y, z) << endl;
						cout << "UF: " << UF.at(x, y, z)
							<< " VF: " << VF.at(x, y, z)
							<< " WF: " << WF.at(x, y, z) << endl;

						cout << "Force: " << Force.at(x, y, z, 0)
							<< " Force: " << Force.at(x, y, z, 1)
							<< " Force: " << Force.at(x, y, z, 2) << endl;
						for (int k = 0; k < Q; k++)
							cout << "k: " << k << " F: " <<F.at(x, y, z, k) << endl;
						cout << "error" << endl;
					}
					
					
					

					UF.at(x, y, z) = U.at(x, y, z) + Force.at(x, y, z, 0) / _rho;
					VF.at(x, y, z) = V.at(x, y, z) + Force.at(x, y, z, 1) / _rho;
					WF.at(x, y, z) = W.at(x, y, z) + Force.at(x, y, z, 2) / _rho;

					/*
					if (x == 52 && y == 20 && z == 15){
						cout << "error" << endl;
						cout << "x:" << x << " y:" << y << " z:" << z << endl;
						cout << "U: (" << U.at(x, y, z) << " , "
							<< V.at(x, y, z) << " , "
							<< W.at(x, y, z) << " ) " << endl;
						cout << "UF: (" << UF.at(x, y, z) << " , "
							<< VF.at(x, y, z) << " , "
							<< WF.at(x, y, z) << " ) " << endl;
						for (int k = 0; k < Q; k++)
							cout << "k: " << k << " F: " << F.at(x, y, z, k) << endl;
						cout << "pause" << endl;
					}
					*/
				}
				
	}

	inline double feq(int k, int x, int y, int z){
		double eu, uv, feq;
		eu = e[k][0] * U.at(x, y, z) 
			+ e[k][1] * V.at(x, y, z) + e[k][2] * W.at(x, y, z);
		uv = U.at(x, y, z) * U.at(x, y, z)
			+ V.at(x, y, z) * V.at(x, y, z)
			+ W.at(x, y, z) * W.at(x, y, z);
		feq = w[k] * rho.at(x, y, z) * 
			(1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * uv);

		/*
		if (feq > 1){
			cout << "uv: " << uv << endl;
			cout << "eu: " << eu << endl;
			cout << "rho: " << rho.at(x, y, z) << endl;
			cout << "w: " << w[k] << endl;
			cout << "U: " << U.at(x, y, z)
				<< " V: " << V.at(x, y, z) 
				<< " W: " << W.at(x, y, z) << endl;
			cout << "UF: " << UF.at(x, y, z)
				<< " VF: " << VF.at(x, y, z)
				<< " WF: " << WF.at(x, y, z) << endl;

			cout << "Force: " << Force.at(x, y, z, 0)
				<< " Force: " << Force.at(x, y, z, 1)
				<< " Force: " << Force.at(x, y, z, 2) << endl;

			cout << "error" << endl;
		}
		*/
			
		return feq;
	}

	inline double feq(int k, int x, int y, int z,double u,double v,double w){
		double eu, uv, feq;
		eu = e[k][0] * u
			+ e[k][1] * v + e[k][2] * w;
		uv = u * u
			+ v * v
			+ w * w;
		feq = this->w[k] * rho.at(x, y, z) *
			(1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * uv);
		return feq;
	}

	double error(){
		double temp1 = 0;
		double temp2 = 0;
//#pragma omp parallel for reduction (+:temp1,temp2)
		for (int x = 0; x < xdim; x++)
			for (int y = 0; y < ydim; y++)
				for (int z = 0; z < zdim; z++){
					temp1
						+= (U0.at(x, y, z) - U.at(x, y, z)) * (U0.at(x, y, z) - U.at(x, y, z))
						+ (V0.at(x, y, z) - V.at(x, y, z)) * (V0.at(x, y, z) - V.at(x, y, z))
						+ (W0.at(x, y, z) - W.at(x, y, z)) * (W0.at(x, y, z) - W.at(x, y, z));
					temp2 
						+= (U.at(x, y, z) + U.at(x, y, z)) * (U.at(x, y, z) + U.at(x, y, z))
						+ (V.at(x, y, z) + V.at(x, y, z)) * (V.at(x, y, z) + V.at(x, y, z))
						+ (W.at(x, y, z) + W.at(x, y, z)) * (W.at(x, y, z) + W.at(x, y, z));
				}
		cout << "temp1: " << temp1 << endl;
		cout << "temp2: " << temp2 << endl;
		return temp1 / (temp2 + 1e-30);
	}

	void output(int m){
		ostringstream name;
		name << "cavity_" << m << ".data";
		ofstream out(name.str().c_str(), ofstream::binary);

		//out << "Title = \"LBM Lid Driven Flow\"\n" << endl;

		/*
		for (int j = 0; j < ydim; j++)
			for (int i = 0; i < xdim; i++)
			{

				out.write((char*)&U.at(i,j,1), sizeof(double));
				out.write((char*)&V.at(i,j,1), sizeof(double));
			}
			*/
		for (int j = 0; j < ydim; j++)
			for (int i = 0; i < xdim; i++)
				out.write((char*)&rho.at(i, j, 15),sizeof(double));
	}

//private:
	

	double w[19];

	double beta = 0.886;
	double neg_beta = -0.886;
	double Fint_const;
	double Fint_const_g1;
	double Fint_const_g2;
	double g;
	double g1;
	double g2;
	double tau = 1.0;
	double c0 = 6.0;

	int Q = 19;

	int xdim, ydim, zdim;
	int xdim_1, ydim_1, zdim_1;


	Lattice3D U;//xdirection
	Lattice3D U0;//xdirection
	Lattice3D UF;//xdirection
	Lattice3D V;//ydirection
	Lattice3D V0;//ydirection
	Lattice3D VF;//ydirection
	Lattice3D W;//zdirection
	Lattice3D W0;//zdirection
	Lattice3D WF;//zdirection
	Lattice3D rho;
	Lattice4D Force;
	Lattice4D F;
	Lattice4D FC;
	Lattice4D F0;
	Lattice2D gs;
};

#endif//__DROPLET_DYNAMICS_H__