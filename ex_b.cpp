#include <iostream>
#include <armadillo>
#include <math.h>
#include <fstream>
#include <string>
using namespace std;
using namespace arma;


double source_func(double x);
vec alg_gen_diag(vec a, vec b, vec c, vec d,int n);
void write_to_file(int n, vec v,double h);

int main(){
  for (int n=10;n<=1000;n*=10){ //Increasing number of grid points.
    vec a (n-1), b(n), c(n-1), d(n), v(n+2); //initializing arrays
    
    a = a*0 - 1;  //Assigning matrix elements.
    b = b*0 + 2;
    c = c*0 - 1;

    v(0) = 0; //boundary conditions on v
    v(n+1) = 0; 

    double h = 1./(n+1); //initializing and assigning step value h
    
    for (double i=0;i<n;i++){d(i) = source_func(i*h);} //Assigning source function values
    
    vec v_mid; //initializing the solved vector values in the middle

    v_mid = alg_gen_diag(a,b,c,d,n); //Solving the linear equation.
    for (int i=1;i<=n;i++){v(i)=v_mid(i-1);}
    
    write_to_file(n,v,h);
    
  }
  return 0;
}

double source_func(double x){
  return 100 * exp(-10*x);}

vec alg_gen_diag(vec a, vec b, vec c, vec d,int n){
  /*
    Solves the linear algebra problem Av = d where the nxn matrix A consists of
    three diagonal arrays a,b and c of length n-1, n and n-1 respectively. 
    
    Function returns the solution vector v. 
  */
  double k;
  vec v(n);
  for (int i = 1;i<=n-1;i++){
    k = a(i-1)/b(i-1);
    a(i-1) = 0;
    b(i) = b(i) - k*c(i-1);
    d(i) = d(i) - k*d(i-1);
  }
  v(n-1) = d(n-1)/b(n-1);
  for (int i=n-1;i>=1;i--){
    v(i-1) = (d(i-1)-c(i-1)*v(i))/b(i-1);
  }
  return v;
}

void write_to_file(int n, vec v,double h){
  ostringstream ostr;
  ostr << n;
  string numberstring = ostr.str();
  string filename = string("n") +numberstring + "diag_method.out";
  ofstream myfile;
  myfile.open(filename.c_str());
  myfile << "x,v \n";
  for (int i=0;i<=n+1;i++){
    myfile << i*h << ',' << v(i) << '\n';
  }
  myfile.close();
  return;
}