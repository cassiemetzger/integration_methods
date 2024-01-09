#include<iostream> 
#include<tuple> 

template <typename F>
double gauss_quadrature_10(const F& f, double a, double b){
    //this method seems more accurate that simpson/adaptive
//x coords
    double x0 = -0.9739065285171717; 
    double x9 = -x0; 
    double x1 = -0.8650633666889845; 
    double x8 = -x1; 
    double x2 = -0.6794095682990244; 
    double x7 = -x2; 
    double x3 = -0.4333953941292472; 
    double x6 = -x3; 
    double x4 = -0.1488743389816312; 
    double x5 = -x4;
    //creating a list/array of values
    double xs[] = {x0, x1, x2, x3, x4, x5, x6, x7, x8, x9}; 
//weights
    double w0=  0.0666713443086881; 
    double w9 = w0; 
    double w1 = 0.1494513491505806; 
    double w8 = w1; 
    double w2 = 0.2190863625159821; 
    double w7 = w2; 
    double w3 = 0.2692667193099963; 
    double w6 = w3; 
    double w4 = 0.2955242247147529; 
    double w5 = w4; 
    //creating list of weight values so that they can be iterated through
    double weights[] = {w0, w1, w2, w3, w4, w5, w6, w7, w8, w9};
    
    double sum = 0.0; 
    for (int i = 0; i < 10; i++){
        double x = 0.5*(b-a)*xs[i] + 0.5*(b+a); 
        sum = sum + weights[i]*f(x); 

    }
    double result = 0.5*(b-a)*sum; 
    return result; 

}
