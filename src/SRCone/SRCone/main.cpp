#include "Droplet_Dynamics.h"
#include <iostream>

/*
const int e[19][3] = {
	{ 0, 0, 0 }, //
	{ 1, 0, 0 },//
	{ -1, 0, 0 }, //
	{ 0, 1, 0 },//
	{ 0, -1, 0 }, //
	{ 0, 0, 1 },//
	{ 0, 0, -1 }, //
	{ 1, 1, 0 },//
	{ -1, -1, 0 }, //
	{ 1, -1, 0 },//
	{ -1, 1, 0 }, //
	{ 1, 0, 1 },//
	{ -1, 0, -1 }, //
	{ 1, 0 - 1 },//
	{ -1, 0, 1 }, //
	{ 0, 1, 1 },//
	{ 0, -1, -1 }, //
	{ 0, 1, -1 },//
	{ 0, -1, 1 }//
};
*/

int main(){
	Droplet_Dynamics model;
	for (int i = 0; i < 100000; i++){
		model.update();
		cout << "round:" << i << endl;
		/*
		cout << "round:" << i << endl;
		cout << "speed ( " << model.U.at(52, 20, 15) << " , " << model.V.at(52, 20, 15) << " , " << model.W.at(52, 20, 15) << " )"
			<< "rho: " << model.rho.at(52, 20, 15) << endl;
		cout << "Force ( " << model.Force.at(52, 20, 15,0) << " , " << model.Force.at(52, 20, 15,1) 
			<< " , " << model.Force.at(52, 20, 15,2) << " )"
			<< endl;
		for (int k = 0; k < 19; k++)
			cout << "k: " << k << " F: " << model.F.at(52, 20, 15, k) << endl;
			*/
		if (i % 5 == 0)
		model.output(i);
		//if (i % 100 == 0){
			cout << model.error() << endl;
		//}
	}
	return 0;
}