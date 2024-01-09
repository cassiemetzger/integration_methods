#pragma once
#include <iostream>
#include <vector>

class runge_kutta_4 {
public:
  runge_kutta_4(int num_equations) {//we tell it the number of equations we want to solve and resize our k as necessary 
    // constructor, initialize the data members
    n_eq = num_equations;
    y_tmp.resize(n_eq);
    k1.resize(n_eq);
    k2.resize(n_eq);
    k3.resize(n_eq);
    k4.resize(n_eq);
    //these are intermediate values, we need to define them in our constructor so that we can use them 
  }

  template <typename F, typename StopCondition>
  std::vector<double> integrate(const F &f, const StopCondition &stop, double h,
                                double x0, const std::vector<double> &y0) {
    // clear the arrays so that we start fresh every time we call this function
    xs.clear();
    result.clear();
//set initial conditions 
    double x = x0;
    std::vector<double> y = y0;
    while (stop(x, y) == false) {
      // append x and y to the result array
      xs.push_back(x);
      result.push_back(y);

      // compute the next step in y
      y = step(f, h, x, y);
      x += h;
    }
    return y;
  }

  template <typename F>
  std::vector<double> step(const F& f, double h, double x, const std::vector<double> &y) {
    //calculating k and y values for each function 
    k1 = f(x, y);
    for (int i = 0; i < n_eq; i++) {
      y_tmp[i] = y[i] + 0.5 * h * k1[i];
    }
    k2 = f(x + 0.5 * h, y_tmp);


    std::vector<double>y3(n_eq);
    for(int i = 0; i < n_eq; i++){
      y_tmp[i] = y[i] + 0.5 * h* k2[i];
    }
    k3 = f(x+0.5*h, y_tmp);


    std::vector<double>y4(n_eq);
    for(int i = 0; i < n_eq; i++){
      y_tmp[i] = y[i] + h*k3[i];
    }
    k4 = f(x + h, y_tmp);


  //computing averages of k1,k2,k3,k4
    std::vector<double> y_next(n_eq);
    for(int i = 0; i < n_eq; i++){
      y_next[i] = y[i] + h*(k1[i]/6.0 + k2[i]/3.0 + k3[i]/3.0 + k4[i]/6.0);
    }

    return y_next; 
  }

  int n_eq;
  std::vector<double> k1, k2, k3, k4, y_tmp;
  std::vector<double> xs;
  std::vector<std::vector<double>> result;
};
