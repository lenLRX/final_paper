#include "Droplet_Dynamics.h"
#include <iostream>

int e[19][3] = {
	{ 0, 0, 0 }, { 1, 0, 0 },
	{ -1, 0, 0 }, { 0, 1, 0 },
	{ 0, -1, 0 }, { 0, 0, 1 },
	{ 0, 0, -1 }, { 1, 1, 0 },
	{ -1, -1, 0 }, { 1, -1, 0 },
	{ -1, 1, 0 }, { 1, 0, 1 },
	{ -1, 0, -1 }, { 1, 0 - 1 },
	{ -1, 0, 1 }, { 0, 1, 1 },
	{ 0, -1, -1 }, { 0, 1, -1 },
	{ 0, -1, 1 }
};

int main(){
	Droplet_Dynamics model;
	for (int i = 0; i < 100000; i++){
		model.update();
		cout << "round:" << i << endl;
		//if (i % 100 == 0){
			cout << model.error() << endl;
		//}
	}
	return 0;
}