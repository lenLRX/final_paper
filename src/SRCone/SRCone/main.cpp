#include "Droplet_Dynamics.h"
#include <iostream>
#include <iomanip>

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
	ofstream result_file("contact_angle.record");
	int total = 20;
	for (int r = 0; r < total; r++){
		double gs =  -2.5 - r * 0.1;
		Droplet_Dynamics model(50, 50, 40,gs);
		double last_contact_angle = 0;
		double last_model_error = 0;
		for (int i = 0;; i++){
			model.update();
			if (i % 200 == 0)
			{
				cout << "gs: " << gs << " round:" << i << endl;
				double model_error = model.error();
				double model_relative_error = abs(last_model_error - model_error) / model_error;
				cout << setprecision(12) << "model_error: " << model_error
				<< " relative error: " << model_relative_error << endl;
				last_model_error = model_error;

				double the_contact_angle = model.calc_contact_angle();
				double error = abs(last_contact_angle - the_contact_angle) / the_contact_angle;
				cout << setprecision(12) << "contact angle: " << the_contact_angle 
					<< " last contact angle: " << last_contact_angle << " relative error: " << error << endl;
				last_contact_angle = the_contact_angle;
				if (model_error<0.01 && error < 1e-5&& !isnan(model_relative_error) && model_error != 0)
				{
					result_file << gs << " , " << the_contact_angle << endl;
					break;
				}
					
			}
			if (i % 1000 == 0)
            {
				model.output(i);
			}
		}
		
	}
	return 0;
}