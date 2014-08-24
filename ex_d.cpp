#include <iostream>
#include <math.h>
#include <fstream>
#include <time.h>
#include <armadillo>
using namespace std;
using namespace arma;
#include "functions.h"

int main(){
  mat A; //Initialize arrays and matrices.
  mat P,L,U;
  vec a,b,c,d,u_mid,vmid_LU,vmid,err_myalg,err_LU,err_myalg2,err_LU2;
  double x,h;
  long double time_LU,time_myalg;
  clock_t start_LU, end_LU, start_myalg, end_myalg;
  ofstream myfile;
  double max_ard,max_ard_LU,eps,eps_LU;
  vec abs_rel_diff, abs_rel_diff_LU;
  long double k_iterations = 5000; //for precision time.
  myfile.open("ex_d.out");
  
  for (int N=10;N<=1000;N*=10){ //Different matrix dimensions.
    cout << N << endl;
    
    //Construct diagonal vectors
    a = zeros(N-1) - 1;  
    b = zeros(N) + 2;
    c = zeros(N-1) - 1;
    
    //Allocate memory
    d.set_size(N);
    u_mid.set_size(N);
    
    //initializing and assigning step value h
    h = 1./(N+1); 

    //Assigning function values to d and exact function
    for (double i=0;i<N;i++){
      x = i*h;
      d(i) = pow(h,2)*source_func(x);
      u_mid(i) = exact_func(x+h);} 

    //Solving the problem with LU-decomposition
    A = zeros(N,N);
    
    for (int i=0;i<N;i++){ //Assigning matrix elements.
      if (i == 0){
	A(i,i) = b(i);
	A(i,i+1) = c(i);}
      else if (i==N-1){
	A(i,i-1) = a(i-1);
	A(i,i) = b(i);}
      else {
	A(i,i-1) = a(i-1);
	A(i,i) = b(i);
	A(i,i+1) = c(i);
      }
    }
    
    //LU-decomposition method
    
    start_LU = clock(); //Clock stamp

    for (int k=0;k<=k_iterations;k++){
    lu(L,U,P,A); //Decomp of A matrix. P = I because easy matrix.
    vmid_LU = solve(trimatu(U), solve(trimatl(L),d)); //corresponds to forwards and backwards substitution
    }

    end_LU = clock(); //Clock end stamp
    time_LU = (end_LU-start_LU)/(k_iterations*(long double) CLOCKS_PER_SEC);

    //My algorithm method
    
    start_myalg = clock(); //Clock stamp
    
    for (int k=0;k<=k_iterations;k++){
    vmid = alg_gen_diag(a,b,c,d,N); //Solving the equation using my algorithm
    }
    
    end_myalg = clock(); //Clock end stamp
    time_myalg = (end_myalg-start_myalg)/(k_iterations*(double) CLOCKS_PER_SEC);
    
    //Calculate errors (MARD)
    abs_rel_diff = abs((vmid-u_mid)/u_mid);
    max_ard = max(abs_rel_diff);
    eps = log10(max_ard);
    abs_rel_diff_LU = abs((vmid_LU-u_mid)/u_mid);
    max_ard_LU = max(abs_rel_diff_LU);
    eps_LU = log10(max_ard_LU);

    //Write to file
    myfile << "Matrix dimension N: " << N << ", Time LU: " << time_LU << "s, MARD LU: " << eps_LU << ", Time MyAlg: " << time_myalg << "s, MARD myalg: " << eps << endl;
    
  }
  myfile.close();
}
