#include<tuple> 
#include<iostream> 
#include<cmath> 

template <typename F>
double trapezoidal_next(const F& f, double a, double b, double dx, double Sn){
    //given Sn = Sn_1, it's one lower than we want 
    //calculate the next value, we have the given dx value 
    double length = b-a; 
    double summy = 0; 
    double n = length/dx; 
    for(int i = 1; i < n; i = i+2){
        summy = summy + f(a + (i*dx));  
    }
     
    double Sn_next = (1.0/2.0)*Sn + (1.0/n)*length*summy; 
    return Sn_next; 
}


template <typename F>
std::tuple<bool, double> adaptive_simpson(const F& f, double a, double b, double tolerance){
    double Sn_1 = 0.5*(f(b)+f(a))*(b-a); 
    double Sn_1_prime = (4.0/3.0)*Sn_1;
    int maxIteration = 15; 
    double length = (b-a); 
    for(int i = 1; i < maxIteration; i++){
        double dx = (length)/(std::pow(2.0,i)); 
        double Sn = trapezoidal_next(f, a,b,dx, Sn_1);
        double Sn_prime = (4.0/3.0)*Sn - (1.0/3.0)* (Sn_1); 
        double error = (1.0/15.0)*(Sn_prime - Sn_1_prime); 
        if(std::abs(error/Sn_prime) < tolerance){
            return std::make_tuple(true, Sn_prime);;
        }
        Sn_1 = Sn; 
        Sn_1_prime = Sn_prime; 

    }
    return std::make_tuple(false, Sn_1_prime);
}




