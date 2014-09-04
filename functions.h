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

double source_func(double x){
  return 100 * exp(-10*x);}

double exact_func(double x){
  return 1 - (1-exp(-10))*x -exp(-10*x);}


