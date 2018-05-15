#include <iostream>
#include "../../headers/arithmeticvector.h"


int main(){
	celerium::ArithmeticVectorN<3 ,double> tst1;
	celerium::ArithmeticVectorN<4 ,double> tst2 = 1.0;
	celerium::ArithmeticVectorN<5 ,double> tst3({1.0,2.0,3.0,4.0,5.0});
	celerium::ArithmeticVectorN<10,   int> tst4({1,2,3,4,5});
	celerium::ArithmeticVector		tst5({1.0,2.0,3.0});
	celerium::ArithmeticVector		tst6({2.0,3.0,4.0});

	std::cout<<tst1<<std::endl;
	std::cout<<tst2<<std::endl;
	std::cout<<tst3<<std::endl;
	std::cout<<tst4<<std::endl;
	std::cout<<"Testing addition: "<<std::endl;
	std::cout<<tst5<<" + "<<tst6<<" = "<<tst5 + tst6<<std::endl;
	std::cout<<"Testing substraction: "<<std::endl;
	std::cout<<tst5<<" - "<<tst6<<" = "<<tst5 - tst6<<std::endl;
	std::cout<<"Testing scalar product: "<<std::endl;
	std::cout<<tst5<<" * "<<tst6<<" = "<<tst5 * tst6<<std::endl;
	std::cout<<"Testing cross product (operators precedence!): "<<std::endl;
	std::cout<<tst5<<" ^ "<<tst6<<" = "<<(tst5 ^ tst6)<<std::endl;
	std::cout<<"Testing division by a scalar: "<<std::endl;
	std::cout<<tst5<<" / 2.0 = "<<tst5/2.0<<std::endl;
	std::cout<<"Testing division by different scalar type: "<<std::endl;
	std::cout<<tst5<<" / 2 = "<<tst5/2<<std::endl;
	std::cout<<"Testing right multiplication: "<<std::endl;
	std::cout<<tst5<<" * 2.0 = "<<tst5*2.0<<std::endl;
	std::cout<<"Testing left multiplication: "<<std::endl;
	std::cout<<"2.0 * "<<tst5<<" = "<<2.0*tst5<<std::endl;
	std::cout<<"Testing multiplication by different scalar type: "<<std::endl;
	std::cout<<tst5<<" * 2 = "<<tst5*2<<std::endl;
	std::cout<<"Testing normalization: "<<std::endl;
	std::cout<<tst5<<".n = "<<tst5.versor()<<std::endl;
	std::cout<<"Testing length: "<<std::endl;
	std::cout<<"||"<<tst5<<"||**2 = "<<tst5.length_squared()<<std::endl;
	std::cout<<"Testing normalization length: "<<std::endl;
	std::cout<<"||"<<tst5<<".n|| = "<<tst6.versor().length()<<std::endl;

	
	celerium::ArithmeticVectorN<3,       int> tst7({1,2,3});
	celerium::ArithmeticVectorN<3,     float> tst8({1.f,2.f,3.f});
	celerium::ArithmeticVectorN<3,      long> tst9({1L,2L,3L});
	celerium::ArithmeticVectorN<3, long long> tst10({1LL,2LL,3LL});

	return 0;
}
