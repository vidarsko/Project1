#include <iostream>
#include <armadillo>
#include <math.h>
#include <fstream>
#include <time.h>
using namespace std;
using namespace arma;
#include "functions.h"

int main(){
  mat A; //Initialize arrays and matrices.
  mat P,L,U;
  vec a,b,c,d,u,u_mid,v_LU,vmid_LU,v,vmid,err_myalg,err_LU,err_myalg2,err_LU2;
  double x,h,time_LU,time_myalg,mes_myalg,mes_LU;
  clock_t start_LU, end_LU, start_myalg, end_myalg;
  ofstream myfile;
  double abs_rel_diff, abs_rel_diff_LU,max_ard,max_ard_LU,eps,eps_LU;
  myfile.open("ex_d.out");
  
  for (int n=1;n<=3;n++){ //Different matrix dimensions.
    
    //Matrix dimension
    int N = pow(10,n);
    
    //Construct diagonal vectors
    a = zeros(N-1) - 1;  
    b = zeros(N) + 2;
    c = zeros(N-1) - 1;
    
    //Allocate memory
    d.set_size(N);
    v.set_size(N+2);
    v_LU.set_size(N+2);
    u.set_size(N+2);
    u_mid.set_size(N);

    //boundary conditions on v
    v(0) =  0;
    v_LU(0) = 0;
    v(N+1) = 0;
    v_LU(N+1) = 0;
    u(0) = 0;
    u(N+1) = 0;
    
    //initializing and assigning step value h
    h = 1./(n+1); 

    //Assigning function values to d and exact function
    for (double i=0;i<n;i++){
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

    lu(L,U,P,A); //Decomp of A matrix. P = I because easy matrix.
    vmid_LU = solve(trimatu(U), solve(trimatu(L),P*b)); //corresponds to forwards and backwards substitution.
    for (int i=1;i<=N;i++){v_LU(i)=vmid_LU(i-1);}
    
    end_LU = clock(); //Clock end stamp
    time_LU = (end_LU-start_LU)/CLOCKS_PER_SEC;

    //My algorithm method
    
    start_myalg = clock(); //Clock stamp
    
    vmid = alg_gen_diag(a,b,c,d,N); //Solving the equation using my algorithm
    for (int i=1;i<=N;i++){v(i)=vmid(i-1);}

    end_myalg = clock(); //Clock end stamp
    time_myalg = (end_myalg-start_myalg)/CLOCKS_PER_SEC;
    
    //Calculate errors (MARD)
    //abs_rel_diff = abs((v-u)/u);
    //max_ard = max(abs_rel_diff);
    //eps = log10(max_ard);
    abs_rel_diff_LU = abs((v_LU-u)/u);
    //max_ard_LU = max(abs_rel_diff_LU);
    //eps_LU = log10(max_ard_LU);

    //Write to file
    myfile << "Matrix dimension N: " << N << ", Time LU: " << time_LU << "s, MARD LU: " << eps_LU << ", Time MyAlg: " << time_myalg << "s, MARD myalg: " << eps << endl;
    
  }
  myfile.close();
}
