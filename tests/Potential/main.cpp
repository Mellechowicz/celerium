#include <iostream>
#include <string>
#include <functional>
#include <utility>
#include <vector>
#include <fstream>
#include "../../headers/potential.h"

float string_to_float(const std::string&& s){
	return std::stof(s);
}


int main(){
	// Formal definition of potential generator
	celerium::Potential<float,float(*)(const std::string&&)> potential_formal("Cr.UPF", string_to_float);

	// Alternatively - using default std::function<T(const std::string&&)>:
	std::function F= [](const std::string& s)->double{return std::stod(s);};
	  celerium::Potential<double> potential_double("Cr.UPF", F);

	std::function G= [](const std::string& s)->float{return std::stof(s);};
	  celerium::Potential<float>  potential_float ("Cr.UPF", G);

	std::function H= [](const std::string& s)->int{return std::stoi(s);};
	  celerium::Potential<int>    potential_int   ("Cr.UPF", H);

	// Special possibility for double-based generators:
	celerium::Potential<double> potential_double_short("Cr.UPF");	// constructor without converter works only for double
	celerium::Potential<double> potential_double_shorter;		// can not be used until file is provided 
	celerium::LocalPotential    potential_double_shortest;		// can not be used until file is provided 

	potential_double_shortest.input("Cr.UPF");				// providing input file (works also for std::strings)

	// Getting data from potential	
	std::vector<double> mesh,local,rho,volume;
	potential_double_shortest.get_mesh(mesh);
	potential_double_shortest.get_local(local);
	potential_double_shortest.get_charge_density(rho);
	potential_double_shortest.get_volume(volume);

	auto q1 = potential_double_shortest.get_charge(1.0);
	auto q0 = potential_double_shortest.get_charge();
	std::cout<<"Charge inside ball of r=1.0 A: "<<q1<<"\te"<<std::endl;
	std::cout<<"Total charge:                  "<<q0<<"\te"<<std::endl;

	std::ofstream out("data.dat");
	for(unsigned i = 0U; i< mesh.size(); ++i)
		out<<mesh[i]<<" "<<volume[i]<<" "<<local[i]<<" "<<rho[i]<<std::endl;
	out.close();

	// Error: broken file
/*
 *  celerium::Potential<double> broken_file("broken.upf",F);
 *  broken_file.get_mesh(mesh);
 *
 */

	// Error: no file
/*
 *  celerium::Potential<double> no_file("Coir.UPF");
 *  no_file.get_mesh(mesh);
 *
 */

	return 0;
}
