#ifndef _POTENTIAL_H
#define _POTENTIAL_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <regex>
#include <functional>
#include <stdexcept>
#include <algorithm>

#ifdef _VERBOSE
  #include <chrono>
#endif

namespace celerium{

template<typename T, typename stoT=std::function<double(const std::string&&)>,size_t RESERVE=2000>
class Potential{
protected:
	std::vector<T>  r;	// spherical coordinate: r
	std::vector<T> dr;	// spherical volume element 
	std::vector<T>  V;	// potential value along r
	std::vector<T>  q;	// charge density along r
	std::string file_name;	// name of data file
	std::ifstream file;	// input file
	bool was_file_read;
	bool was_file_set;
	stoT convert;

	// parser
	std::size_t read_file(){
	  if(was_file_set){
#ifdef _VERBOSE
	  std::cerr<<"celerium::Potential:\tReading..."<<std::flush;
	  auto start = std::chrono::system_clock::now();
#endif
	     r.clear();
	    dr.clear();
	     V.clear();
	     q.clear();
	     r.reserve(RESERVE);
	    dr.reserve(RESERVE);
	     V.reserve(RESERVE);
	     q.reserve(RESERVE);

	    std::regex   open_r( "<PP_R[ >]",       std::regex_constants::ECMAScript | std::regex_constants::icase);
	    std::regex  close_r("</PP_R[ >]",       std::regex_constants::ECMAScript | std::regex_constants::icase);
	    std::regex  open_dr( "<PP_RAB[ >]",     std::regex_constants::ECMAScript | std::regex_constants::icase);
	    std::regex close_dr("</PP_RAB[ >]",     std::regex_constants::ECMAScript | std::regex_constants::icase);
	    std::regex   open_V( "<PP_LOCAL[ >]",   std::regex_constants::ECMAScript | std::regex_constants::icase);
	    std::regex  close_V("</PP_LOCAL[ >]",   std::regex_constants::ECMAScript | std::regex_constants::icase);
	    std::regex   open_q( "<PP_RHOATOM[ >]", std::regex_constants::ECMAScript | std::regex_constants::icase);
	    std::regex  close_q("</PP_RHOATOM[ >]", std::regex_constants::ECMAScript | std::regex_constants::icase);

	    bool  in_r = false;
	    bool in_dr = false;
	    bool  in_V = false;
	    bool  in_q = false;

	    size_t  r_length = 0;
	    size_t dr_length = 0;
	    size_t  V_length = 0;
	    size_t  q_length = 0;

	    file = std::ifstream(file_name.c_str());

	    if (file.is_open()){
	      std::string line;

	      while (std::getline(file,line)){
		if     (std::regex_search( line,  close_r))  in_r = false;
		else if(std::regex_search( line,   open_r))  in_r =  true;
		else if(std::regex_search( line, close_dr)) in_dr = false;
		else if(std::regex_search( line,  open_dr)) in_dr =  true;
		else if(std::regex_search( line,  close_V))  in_V = false;
		else if(std::regex_search( line,   open_V))  in_V =  true;
		else if(std::regex_search( line,  close_q))  in_q = false;
		else if(std::regex_search( line,   open_q))  in_q =  true;
		else{
		  if      (in_r)  r_length += input_line( r,line);
		  else if(in_dr) dr_length += input_line(dr,line);
		  else if (in_V)  V_length += input_line( V,line);
		  else if (in_q)  q_length += input_line( q,line);
		}
	      }

    	    file.close();

	    } else throw std::runtime_error("celerium::Potential: Unable to open file."); 

	    if ( r_length != dr_length )
	      throw std::logic_error(std::string("celerium::Potential: Broken data:  r_length= ")
		+std::to_string( r_length)+std::string(" while dr_length= ")+std::to_string(dr_length));
	    if ( r_length !=  V_length )
	      throw std::logic_error(std::string("celerium::Potential: Broken data:  r_length= ")
		+std::to_string( r_length)+std::string(" while  V_length= ")+std::to_string( V_length));
	    if ( r_length !=  q_length )
	      throw std::logic_error(std::string("celerium::Potential: Broken data:  r_length= ")
		+std::to_string( r_length)+std::string(" while  q_length= ")+std::to_string( q_length));
	    if (dr_length !=  V_length )
	      throw std::logic_error(std::string("celerium::Potential: Broken data: dr_length= ")
		+std::to_string(dr_length)+std::string(" while  V_length= ")+std::to_string( V_length));
	    if (dr_length !=  q_length )
	      throw std::logic_error(std::string("celerium::Potential: Broken data: dr_length= ")
		+std::to_string(dr_length)+std::string(" while  q_length= ")+std::to_string( q_length));
	    if ( V_length !=  q_length )
	      throw std::logic_error(std::string("celerium::Potential: Broken data:  V_length= ")
		+std::to_string( V_length)+std::string(" while  q_length= ")+std::to_string( q_length));

	    was_file_read = true;
	  
#ifdef _VERBOSE
	  auto stop = std::chrono::system_clock::now();
	  std::cerr<<" success\n\t\t\tRead "<<r_length<<" records from file "<<file_name<<","<<std::endl;
	  std::cerr<<"\t\t\tin "<<std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()<<" ms."<<std::endl;
#endif
	    return r_length;
	  }

#ifdef _VERBOSE
	  std::cerr<<" failed\n\t\t\tUnable to open: file_name not given!"<<std::endl;
#endif
	  return -1;
	}

	// converting line of numbers to std::vector elements
	size_t input_line(std::vector<T>& data, std::string line){
	  std::regex floating_point("[\\d\\.eE+-]+");

	  auto begin = std::sregex_iterator(line.begin(), line.end(), floating_point);
	  auto end   = std::sregex_iterator();

	  for (std::sregex_iterator it = begin; it != end; ++it) 
	    data.push_back(convert(it->str()));

	  return std::distance(begin,end);
	}

public:
	Potential() = delete;
	Potential(const std::string&& input_file) = delete;
	Potential(const char* input_file) = delete;

	Potential(stoT&& converter):was_file_read(false),was_file_set(false){
		convert = converter;
	}
	Potential(const std::string&& input_file, stoT&& converter):was_file_read(false){
		convert = converter;
		file_name = input_file;
		was_file_set = true;
	}
	Potential(const char* input_file, stoT&& converter){
		convert = converter;
		file_name = std::string(input_file);
		was_file_set = true;
	}
	virtual ~Potential(){}

	/*
	 *
	 * GETTERS
	 *
	 */
	// mesh getter
	std::size_t get_mesh(std::vector<T>&  vec){
	  if(!was_file_read) read_file();

#ifdef _VERBOSE
	  if(vec.size()) std::cerr<<"celerium::Potential::get_mesh: you passed a non-zero vector!"<<std::endl;
	  std::cerr<<"celerium::Potential::get_mesh: Passing..."<<std::flush;
	  auto start = std::chrono::system_clock::now();
#endif
	  vec.clear();
	  vec = r;

#ifdef _VERBOSE
	  auto stop = std::chrono::system_clock::now();
	  std::cerr<<" success\n\t\t\tPassed "<<vec.size()<<" records,"<<std::endl;
	  std::cerr<<"\t\t\tin "<<std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()<<" µs."<<std::endl;
#endif
	  return vec.size();
	}

	// spherical volume getter
	std::size_t get_volume(std::vector<T>&  vec){
	  if(!was_file_read) read_file();

#ifdef _VERBOSE
	  if(vec.size()) std::cerr<<"celerium::Potential::get_volume: you passed a non-zero vector!"<<std::endl;
	  std::cerr<<"celerium::Potential::get_volume: Passing..."<<std::flush;
	  auto start = std::chrono::system_clock::now();
#endif
	  vec.clear();
	  vec = dr;

#ifdef _VERBOSE
	  auto stop = std::chrono::system_clock::now();
	  std::cerr<<" success\n\t\t\tPassed "<<vec.size()<<" records,"<<std::endl;
	  std::cerr<<"\t\t\tin "<<std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()<<" µs."<<std::endl;
#endif
	  return vec.size();
	}

	// local potential getter
	std::size_t get_local(std::vector<T>&  vec){
	  if(!was_file_read) read_file();

#ifdef _VERBOSE
	  if(vec.size()) std::cerr<<"celerium::Potential::get_local: you passed a non-zero vector!"<<std::endl;
	  std::cerr<<"celerium::Potential::get_local: Passing..."<<std::flush;
	  auto start = std::chrono::system_clock::now();
#endif
	  vec.clear();
	  vec = V;

#ifdef _VERBOSE
	  auto stop = std::chrono::system_clock::now();
	  std::cerr<<" success\n\t\t\tPassed "<<vec.size()<<" records,"<<std::endl;
	  std::cerr<<"\t\t\tin "<<std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()<<" µs."<<std::endl;
#endif
	  return vec.size();
	}

	// charge density getter
	std::size_t get_charge_density(std::vector<T>&  vec){
	  if(!was_file_read) read_file();

#ifdef _VERBOSE
	  if(vec.size()) std::cerr<<"celerium::Potential::get_charge_density: you passed a non-zero vector!"<<std::endl;
	  std::cerr<<"celerium::Potential::get_charge_density: Passing..."<<std::flush;
	  auto start = std::chrono::system_clock::now();
#endif
	  vec.clear();
	  vec = q;

#ifdef _VERBOSE
	  auto stop = std::chrono::system_clock::now();
	  std::cerr<<" success\n\t\t\tPassed "<<vec.size()<<" records,"<<std::endl;
	  std::cerr<<"\t\t\tin "<<std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()<<" µs."<<std::endl;
#endif
	  return vec.size();
	}

	/*
	 *
	 * SETTERS
	 *
	 */

	// changing input
	void input(const std::string&& input_file){
#ifdef _VERBOSE
	  std::cerr<<"celerium::Potential::input: Input file changed."<<std::endl;
#endif
		file_name = input_file;
		was_file_set  = true;
		was_file_read = false;
	}

	void input(const char* input_file){
#ifdef _VERBOSE
	  std::cerr<<"celerium::Potential::input: Input file changed."<<std::endl;
#endif
		file_name = std::string(input_file);
		was_file_set  = true;
		was_file_read = false;
	}

}; // end of class Potential

template<>
Potential<double>::Potential(){
		convert = [](const std::string&& s)->double{return std::stod(s);};
		was_file_set  = false;
		was_file_read = false;
}

template<>
Potential<double>::Potential(const char* input_file){
		convert = [](const std::string&& s)->double{return std::stod(s);};
		file_name = std::string(input_file);
		was_file_set  = true;
		was_file_read = false;
	}

template<>
Potential<double>::Potential(const std::string&& input_file){
		convert = [](const std::string&& s)->double{return std::stod(s);};
		file_name = input_file;
		was_file_set  = true;
		was_file_read = false;
	}

using LocalPotential = Potential<double,std::function<double(const std::string&&)>>;

} //end of namespace celerium

#endif
