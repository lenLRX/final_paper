#ifndef __LATTICE3D_H__
#define __LATTICE3D_H__
#include <cassert>
#include <algorithm>
#include "LatticeMaterial.h"
/*
z  y
| /
|/
O-------x
*/

using namespace std;
template<typename T>
class Lattice3D
{
public:
	Lattice3D(int xdim, int ydim, int zdim) :xdim(xdim), ydim(ydim), zdim(zdim){
		total_size = xdim * ydim * zdim;
		xy_block = xdim * ydim;
		buffer = new T[total_size];
		clean();
	}
	inline T& at(int x, int y, int z){
		assert(x < xdim);
		assert(y < ydim);
		if (z >= zdim)
			throw 0;
		return *(buffer + z * xy_block + y * xdim + x);
	}

	inline void swap(Lattice3D& other){
		std::swap(buffer, other.buffer);
	}

	inline void clean(){
		fill(buffer, buffer + total_size, 0);
	}

	~Lattice3D(){
		delete[] buffer;
	}
private:
	int xdim;
	int ydim;
	int zdim;
	int xy_block;
	int total_size;
	T* buffer;
};

template<>
inline void Lattice3D<Material>::clean(){
	fill(buffer, buffer + total_size, Material::solid);
}
#endif//__LATTICE3D_H__