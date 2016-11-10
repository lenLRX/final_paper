#ifndef __LATTICE2D_H__
#define __LATTICE2D_H__
#include <cassert>
#include <algorithm>
/*
z  y
| /
|/
O-------x
*/

using namespace std;

template<typename T>
class Lattice2D
{
public:
	Lattice2D(int xdim, int ydim) :xdim(xdim), ydim(ydim){
		total_size = xdim * ydim;
		buffer = new T[total_size];
		clean();
	}
	inline T& at(int x, int y){
		assert(x < xdim);
		assert(y < ydim);
		return *(buffer + y * xdim + x);
	}

	inline void clean(){
		fill(buffer, buffer + total_size, 0);
	}

	~Lattice2D(){
		delete[] buffer;
	}
private:
	int xdim;
	int ydim;
	int total_size;
	T* buffer;
};
#endif//__LATTICE2D_H__