#include <iostream>
#include <armadillo>
#include <math.h>
#include <fstream>
#include <string>
using namespace std;
using namespace arma;
#include "functions.h"

void write_to_file(int n, vec v,double h);

int main(){
  vec a,b,c,d,v; //initializing arrays.
  double h;
  for (int n=10;n<=1000;n*=10){ //Increasing number of grid points.

    a = zeros(n-1) - 1;  //Assigning matrix elements
    b = zeros(n) + 2;
    c = zeros(n-1) - 1;
    d.set_size(n);   //Sizing vectors
    v.set_size(n+2);

    v(0) = 0; //boundary conditions on v
    v(n+1) = 0; 

    h = 1./(n+1); //initializing and assigning step value h

    for (double i=0;i<n;i++){d(i) = pow(h,2)*source_func(i*h);} //Assigning source function values
    
    vec v_mid; //initializing the solved vector values in the middle

    v_mid = alg_gen_diag(a,b,c,d,n); //Solving the linear equation.
    for (int i=1;i<=n;i++){v(i)=v_mid(i-1);}
    
    write_to_file(n,v,h); //write result to file
    
  }
  return 0;
}


void write_to_file(int n, vec v,double h){
  ostringstream ostr;
  ostr << n;
  string numberstring = ostr.str();
  string filename = string("n") +numberstring + "diag_method_x_v.out";
  ofstream myfile;
  myfile.open(filename.c_str());
  //myfile << "x,v" << endl;
  for (int i=0;i<=n+1;i++){
    myfile << i*h << ',' << v(i) << endl;
  }
  myfile.close();
  return;
}
