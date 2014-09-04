#include <iostream>
#include <armadillo>
#include <math.h>
#include <fstream>
#include <string>
using namespace std;
using namespace arma;
#include "functions.h"

int main(){
  vec a,b,c,d,v,u_mid; //initializing arrays.
  double h,max_ard,eps;
  ofstream myfile;
  vec abs_rel_diff;
  myfile.open("h_mabrd.out");
  for (int n=10;n<=10000000;n*=10){ //Increasing number of grid points.
    cout << n<< endl;

    a = zeros(n-1) - 1;  //Assigning matrix elements
    b = zeros(n) + 2;
    c = zeros(n-1) - 1;
    d.set_size(n);   //Sizing vectors
    v.set_size(n+2);
    u_mid.set_size(n);

    v(0) = 0; //boundary conditions on v
    v(n+1) = 0; 
    
    h = 1./(n+1); //initializing and assigning step value h
    double x;
    for (double i=0;i<n;i++){//Assigning source  and exact values.
      x = i*h;
      d(i) = pow(h,2)*source_func(x);
      u_mid(i) = exact_func(x+h);} 
    
    vec v_mid; //initializing the solved vector values in the middle
    v_mid = alg_gen_diag(a,b,c,d,n); //Solving the linear equation.
    for (int i=1;i<=n;i++){v(i)=v_mid(i-1);}

    abs_rel_diff = abs((v_mid-u_mid)/u_mid);
    max_ard = max(abs_rel_diff);
    eps = log10(max_ard);
    myfile << h <<',' << eps << '\n';
    
  }
  myfile.close();
  return 0;
}


