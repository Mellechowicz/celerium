#ifndef ARITHMETICVECTOR_H
#define ARITHMETICVECTOR_H

#include <iostream>
#include <array>
#include <initializer_list>
#include <exception>
#include <cmath>
#include <type_traits>

namespace celerium{

//ArithmeticVectorN = N-dimensional vector of type T 
  template<size_t N,typename T = double>
  class ArithmeticVectorN{
  static_assert(std::is_same<T, double      >::value ||
		std::is_same<T, float       >::value ||
		std::is_same<T, int         >::value ||
		std::is_same<T, unsigned int>::value ||
		std::is_same<T, long        >::value ||
		std::is_same<T, long long   >::value,	"T must be arithmetic");
  protected:
	std::array<T,N> body;
  public:

    ArithmeticVectorN<N,T>(){
    }
    ArithmeticVectorN<N,T>(T value){
		body=std::array<T,N>();
		for(auto& x : body)
			x = value;
    }
    ArithmeticVectorN<N,T>(std::initializer_list<T>& list){
	for(size_t i = 0; i<N; ++i)
	  body[i] = static_cast<T>(list[i]);
    }
    ArithmeticVectorN<N,T>(std::array<T,N> input){
		body=input;
    }

    // FUNCTIONS
    
    size_t size() const{return N;}

    T* data() {return body.data();}
    
    ArithmeticVectorN versor(){
	auto mynorm = norm(*this);

	if(mynorm>1e-14)return (*this)/norm(*this);
	else		return ArithmeticVectorN<N,T>();
    }
    
    void normalize(){
     T mynorm=norm(*this);
     for (auto& it : body)
       it /= mynorm;
    }
    
    
    // OPERATORS
    
    T operator[](size_t n) const{
		if(n>=N || n<0){
		   throw std::domain_error("celerium::ArithmeticVectorN<N,T>::operator[]: there is no field to access!\n");
		}
     return body[n]; 
    }
    
    T& operator[](size_t n){
		if(n>=N || n<0){
		   throw std::domain_error("celerium::ArithmeticVectorN<N,T>::operator[]: there is no field to access!\n");
		}
     return body[n]; 
    }
    
    std::array<T,N> operator()(){
     return body; 
    }
    
     ArithmeticVectorN operator+(const ArithmeticVectorN & right) const {
        std::array<T,N>sum = body;
	for(unsigned int i=0; i<N; i++)
	    sum[i]+=right[i];
        return ArithmeticVectorN(sum);
    }
     ArithmeticVectorN operator-(const ArithmeticVectorN & right) const {
        std::array<T,N>sum = body;
	for(unsigned int i=0; i<N; i++)
	    sum[i]-=right[i];
        return ArithmeticVectorN(sum);
    }
    ArithmeticVectorN operator-() const{
      std::array<T,N>sum = body;
	for(unsigned int i=0; i<N; i++)
	    sum[i] *= -1;
        return ArithmeticVectorN(sum);
    }
    void operator+=(const ArithmeticVectorN & right) {
	for(unsigned int i=0; i<N; i++)
	    body[i]+=right[i];
    }
    void operator-=(const ArithmeticVectorN & right) {
	for(unsigned int i=0; i<N; i++)
	    body[i]-=right[i];
    }
     double operator*(const ArithmeticVectorN & right) const {
        double sum=0;
	for(unsigned int i=0; i<N; i++)
	    sum+=(body[i]*right[i]);
        return sum;
    }
     ArithmeticVectorN operator^(const ArithmeticVectorN & right) const {
	if(N!=3){
	    std::cerr<<("celerium::ArithmeticVectorN<N,T>::operator^: vector product exist only in 3 dimensions!\n");
	    return ArithmeticVectorN<N,T>();
	}
		
        ArithmeticVectorN product(3);
	
	product[0]=body[1]*right[2]-body[2]*right[1];
	product[1]=body[2]*right[0]-body[0]*right[2];
	product[2]=body[0]*right[1]-body[1]*right[0];
	
        return product;
     }
    template<typename T2>
    ArithmeticVectorN operator/(const T2 & scalar) const {
        std::array<T,N>sum = body;
	for(auto& x : sum)
	    x /= scalar;
        return ArithmeticVectorN(sum);
    }
    template<typename T2>
    ArithmeticVectorN operator*(const T2 & scalar) const {
        std::array<T,N>sum = body;
	for(auto& x : sum)
	    x *= scalar;
        return ArithmeticVectorN(sum);
    }
    template<typename T2>
    friend ArithmeticVectorN operator*(const T2 & scalar, const ArithmeticVectorN& right) {
        return right*scalar;
    }
    template<typename T2>
   void operator*=(const T & scalar) {
        for(auto& item : body)
	    item*=scalar;
    }
    template<typename T2>
    void operator/=(const T & scalar) {
        for(auto& item : body)
	    item/=scalar;
    }

    friend std::ostream& operator<<(std::ostream& stream, const ArithmeticVectorN& vector){
     stream<<"[";
     for (const auto& item : vector.body)
       stream<<item<<", ";
     stream<<"\b\b]";
     
     return stream;
    }
    
    T length_squared(){
      T localnorm=T(0);
      for(const auto& item : body){
        localnorm += item*item;
      }
      return localnorm;
    }
    
    T length(){
      return sqrt(length_squared());
    }
    
    // STATIC FUNCTIONS
    static T square_norm(const ArithmeticVectorN& vector){
      T localnorm=T(0);
      for(const auto& item : vector.body){
        localnorm += item*item;
      }
      return localnorm;
    }
    
    static T norm(const ArithmeticVectorN& vector){
      return sqrt(square_norm(vector));
    }
    
};	//end class ArithmeticVectorN
 

//********
// Alias      
//********
 
  using ArithmeticVector  = ArithmeticVectorN<3U,double>;
  using ArithmeticVectorF = ArithmeticVectorN<3U,float>;
  using ArithmeticVectorI = ArithmeticVectorN<3U,int>;
  using ArithmeticVectorU = ArithmeticVectorN<3U,unsigned int>;
  using ArithmeticVectorL = ArithmeticVectorN<3U,long>;
  using ArithmeticVectorLL= ArithmeticVectorN<3U,long long>;

} // end namespace celerium
#endif
