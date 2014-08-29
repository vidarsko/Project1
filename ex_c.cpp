#include <iostream>
#include <armadillo>
#include <math.h>
#include <fstream>
#include <string>
using namespace std;
using namespace arma;

double source_func(double x);
double exact_func(double x);
vec alg_gen_diag(vec a, vec b, vec c, vec d,int n);

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

double source_func(double x){
  return 100 * exp(-10*x);}

double exact_func(double x){
  return 1 - (1-exp(-10))*x -exp(-10*x);}

vec alg_gen_diag(vec a, vec b, vec c, vec d,int n){
  /*
    Solves the linear algebra problem Av = d where the nxn matrix A consists of
    three diagonal arrays a,b and c of length n-1, n and n-1 respectively. 
    
    Function returns the solution vector v. 
  */
  double k;
  vec v;
  v.set_size(n);
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
