#ifndef __LATTICE4D_H__
#define __LATTICE4D_H__
#include <cassert>
#include <algorithm>
/*
z  y
| /
|/
O-------x
*/

using namespace std;

class Lattice4D
{
public:
	Lattice4D(int dim1, int dim2, int dim3, int dim4) 
		:dim1(dim1), dim2(dim2), dim3(dim3), dim4(dim4){
		block1 = dim4;
		block2 = block1 * dim3;
		block3 = block2 * dim2;
		total_size = block3 * dim1;
		buffer = new double[total_size];
		if (buffer == nullptr)
			throw 0;
		clean();
	}
	inline void swap(Lattice4D& other){
		std::swap(buffer, other.buffer);
	}

	inline double& at(int x1, int x2, int x3, int x4){
		assert(x1 < dim1);
		assert(x2 < dim2);
		assert(x3 < dim3);
		assert(x4 < dim4);
		return *(buffer + x1 * block3 
			+ x2 * block2 + x3 * block1 + x4);
	}

	inline void clean(){
		fill(buffer, buffer + total_size, 0);
	}

	~Lattice4D(){
		delete[] buffer;
	}
private:
	int dim1;
	int dim2;
	int dim3;
	int dim4;
	int block1;
	int block2;
	int block3;
	int total_size;
	double* buffer;
};
#endif//__LATTICE4D_H__