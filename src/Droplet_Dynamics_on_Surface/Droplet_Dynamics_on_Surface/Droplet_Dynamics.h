#ifndef __DROPLET_DYNAMICS_H__
#define __DROPLET_DYNAMICS_H__
#include "Lattice4D.h"
#include "Lattice3D.h"
#include "Lattice2D.h"
#include <cmath>
#include <iostream>

extern int e[19][3];

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
		U(xdim, ydim, zdim), U0(xdim, ydim, zdim),
		V(xdim, ydim, zdim), V0(xdim, ydim, zdim),
		W(xdim, ydim, zdim), W0(xdim, ydim, zdim),
		rho(xdim, ydim, zdim), Force(xdim,ydim,zdim,3),F(xdim, ydim, zdim, 19),
		F0(xdim, ydim, zdim, 19), gs(xdim, ydim){
		w[0] = 1.0 / 3.0;
		fill(w + 1, w + 7, 1.0 / 18.0);
		fill(w + 7, w + 19, 1.0 / 36.0);

		double T = 2.925;
		g = -1 / T;
		g1 = g;
		g2 = g / 2;
		Fint_const = (1 - beta) / 2;
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
						+ (z - zdim / 2) * (z - zdim / 2)) <= 12)
						rho.at(x, y, z) = 2.745;
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
	}

	inline double calc_pressure(double _rho){
		double _1_m_exp = 1 - exp(-_rho);
		return 1.0 / 3.0 *_rho + 0.5 * g * _1_m_exp * _1_m_exp;
	}

	inline double Fint(int k,int x, int y, int z){
		int xp = x + e[k][0];
		int yp = y + e[k][1];
		int zp = z + e[k][2];

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
		
		F.swap(F0);
		F.clean();
		Force.clean();

#pragma omp parallel for
		for (int x = 0; x < xdim; x++)
			for (int y = 0; y < ydim; y++)
				for (int z = 1; z < zdim_1; z++)
					for (int k = 0; k < Q; k++){
						for (int dir = 0; dir < 3; dir++){
							if (e[k][dir] != 0)
								Force.at(x, y, z, dir) += e[k][dir] * Fint(k, x, y, z);
						}
					}

		//top boundary
#pragma omp parallel for
		for (int x = 0; x < xdim; x++)
			for (int y = 0; y < ydim; y++){
				int z = zdim_1;
				for (int k = 0; k < Q; k++){
					int xp = x + e[k][0];
					int yp = y + e[k][1];
					int zp = z + e[k][2];

					if (zp > z)
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
					double temp = temp1 + temp2;

					for (int dir = 0; dir < 3; dir++){
						if (e[k][dir] != 0)
							Force.at(x, y, z, dir) += e[k][dir] * temp;
					}
				}
					
			}
		//top

		//bottom boundary
#pragma omp parallel for
		for (int x = 0; x < xdim; x++)
			for (int y = 0; y < ydim; y++){
				int z = 0;
				for (int k = 0; k < Q; k++){
					int xp = x + e[k][0];
					int yp = y + e[k][1];
					int zp = z + e[k][2];

					if (zp < 0)
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
					double temp = temp1 + temp2;

					for (int dir = 0; dir < 3; dir++){
						if (e[k][dir] != 0)
							Force.at(x, y, z, dir) += e[k][dir] * temp;
					}
				}
				Force.at(x, y, z, 2) += Fs(x, y, z);
			}
		//bottom boundary ends
		
#pragma omp parallel for
		for (int x = 0; x < xdim; x++)
			for (int y = 0; y < ydim; y++)
				for (int z = 1; z < zdim_1; z++){

					for (int k = 0; k < Q; k++){
						int xp = x - e[k][0];
						int yp = y - e[k][1];
						int zp = z - e[k][2];

						if (xp < 0)
							xp = xdim_1;
						if (xp > xdim_1)
							xp = 0;

						if (yp < 0)
							yp = ydim_1;
						if (yp > ydim_1)
							yp = 0;

						F.at(x, y, z, k) = F0.at(xp, yp, zp, k)
							+ (feq(k, xp, yp, zp) - F0.at(xp, yp, zp, k)) / tau + feq(k, xp, yp, zp,
							U.at(xp, yp, zp) + Force.at(xp, yp, zp, 0) / rho.at(xp, yp, zp),
							V.at(xp, yp, zp) + Force.at(xp, yp, zp, 1) / rho.at(xp, yp, zp),
							W.at(xp, yp, zp) + Force.at(xp, yp, zp, 2) / rho.at(xp, yp, zp))
							- feq(k, xp, yp, zp);
					}
				}
		U.swap(U0);
		V.swap(V0);
		W.swap(W0);
		rho.clean();
		U.clean();
		V.clean();
		W.clean();

#pragma omp parallel for
		for (int x = 0; x < xdim; x++)
			for (int y = 0; y < ydim; y++)
				for (int z = 1; z < zdim_1; z++){
					for (int k = 0; k < Q; k++){
						rho.at(x, y, z) += F.at(x, y, z, k);
						U.at(x, y, z) += e[k][0] * F.at(x, y, z, k);
						V.at(x, y, z) += e[k][1] * F.at(x, y, z, k);
						W.at(x, y, z) += e[k][2] * F.at(x, y, z, k);
					}
					double _rho = rho.at(x, y, z);
					U.at(x, y, z) /= _rho;
					V.at(x, y, z) /= _rho;
					W.at(x, y, z) /= _rho;
				}

		//top bottom boundary
#pragma omp parallel for
		for (int x = 0; x < xdim; x++)
			for (int y = 0; y < ydim; y++){
				for (int z = 0; z < 1; z++){
					rho.at(x, y, z) = rho.at(x, y, z + 1);
					for (int k = 0; k < Q; k++){
						int xp = x - e[k][0];
						int yp = y - e[k][1];
						int zp = z - e[k][2];

						if (zp < 0)
							continue;

						if (xp < 0)
							xp = xdim_1;
						if (xp > xdim_1)
							xp = 0;

						if (yp < 0)
							yp = ydim_1;
						if (yp > ydim_1)
							yp = 0;


						F.at(x, y, z, k) = F0.at(xp, yp, zp, k)
							+ (feq(k, xp, yp, zp) - F0.at(xp, yp, zp, k)) / tau
							+ feq(k, xp, yp, zp,
							U.at(xp, yp, zp) + Force.at(xp, yp, zp, 0) / rho.at(xp, yp, zp),
							V.at(xp, yp, zp) + Force.at(xp, yp, zp, 1) / rho.at(xp, yp, zp),
							W.at(xp, yp, zp) + Force.at(xp, yp, zp, 2) / rho.at(xp, yp, zp))
							- feq(k, xp, yp, zp);
					}
				}

				for (int z = zdim_1; z < zdim; z++){
					rho.at(x, y, z) = rho.at(x, y, z - 1);
					for (int k = 0; k < Q; k++){
						int xp = x - e[k][0];
						int yp = y - e[k][1];
						int zp = z - e[k][2];

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

						F.at(x, y, z, k) = F0.at(xp, yp, zp, k)
							+ (feq(k, xp, yp, zp) - F0.at(xp, yp, zp, k)) / tau
							+ feq(k, xp, yp, zp,
							U.at(xp, yp, zp) + Force.at(xp, yp, zp, 0) / rho.at(xp, yp, zp),
							V.at(xp, yp, zp) + Force.at(xp, yp, zp, 1) / rho.at(xp, yp, zp),
							W.at(xp, yp, zp) + Force.at(xp, yp, zp, 2) / rho.at(xp, yp, zp))
							- feq(k, xp, yp, zp);
					}
				}
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
#pragma omp parallel for reduction (+:temp1,temp2)
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
		
		return temp1 / (temp2 + 1e-30);
	}

private:
	

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
	Lattice3D V;//ydirection
	Lattice3D V0;//ydirection
	Lattice3D W;//zdirection
	Lattice3D W0;//zdirection
	Lattice3D rho;
	Lattice4D Force;
	Lattice4D F;
	Lattice4D F0;
	Lattice2D gs;
};

#endif//__DROPLET_DYNAMICS_H__