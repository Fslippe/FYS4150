// Creates a tridiagonal matrix
arma::mat create_tridiagonal(int N, double a, double d, double e)
{
  arma::mat A = arma::mat(N, N, arma::fill::eye) * d;
  A(0,1) = e; // first superdiagonal element

  for (int i = 1; i < N-1; i++)
  {
    A(i,i+1) = e; // filling superdiagonal
    A(i,i-1) = a; // filling subdiagonal
  }

  A(N-1,N-2) = a; //last subdiagonal element
  return A;
}

// Creates a symmetric tridiagonal matrix
arma::mat create_symmetric_tridiagonal(int N, double a, double d)
{
  return create_tridiagonal(N, a, d, a);
}

// Calculates analytic eigenvalues
arma::vec analytic_eigenval(int N, double a, double d)
{
  const double pi = 2*acos(0.0);
  arma::vec eig_val = arma::vec(N);
  for(int i = 1; i <= N; ++i)
  {
    eig_val(i-1) = d + 2*a*cos(i*pi/(N+1));
  }
  return eig_val;
}

//Calculates analytic eigenvectors
arma::mat analytic_eigenvector(int N, double a, double d)
{
  const double pi = 2*acos(0.0);
  arma::mat eig_mat = arma::mat(N,N);
  for(int i = 1; i <= N; ++i)
  {
    arma::vec col_v = arma::vec(N);
    for(int j = 1; j <= N; ++j)
    {
      col_v(j-1) = sin((i*j*pi) / (N+1));
    }
    eig_mat.col(i-1) = col_v;
  }
  return eig_mat;
}

// Finding the max off-diagonal element of a matrix
double max_offdiag_symmetric(const arma::mat& A, int& k, int &l)
{
  int N = A.n_rows;
  assert(N > 1);
  assert(A.is_square());

  double maxval = std::abs(A(k,l)); //(i,j) (row, col)
  std::cout << maxval << "\n";

  int col_n = 1;

  for (int i = 1; i <= N-1; ++i) //Loops over rows
  {
    for (int j = 0; j < i; ++j) //loops over columns until the subdiagonal
    {
      if(std::abs(A(i,j)) > maxval)
      {
        k = i; l = j;
        maxval = std::abs(A(i,j));
      }
   }
  }
  return maxval;
 }

 void jacobi_rotate(arma::mat& A, arma::mat& R, int& k, int& l)
 {
    double tau = (A(l,l) - A(k,k)) / (2 * A(k,l));
    double A_ik;
    int N = A.n_rows;
    double t;
    if (tau > 0)
    {
       t = - tau + sqrt(1 + tau*tau);
    }
    else
    {
        t = - tau - sqrt(1 + tau*tau);
    }
    double c = 1 / sqrt(1+t*t);
    double s = c * t;

    A(k,k) = A(k,k)*c*c - 2*A(k,l)*c*s + A(l,l)*s*s;
    A(l,l) = A(l,l)*c*c - 2*A(k,l)*c*s + A(k,k)*s*s;
    A(k,l) = 0;
    A(l,k) = 0;

    for (int i = 0; i <= N-1; i++)
    {
        if (i != k && i != l)
        {
            A_ik = A(i,k);
            A(i,k) = A(i,k)*c -A(i,l)*s;
            A(k,i) = A(i,k);
            A(i,l) = A(i,l)*c + A_ik*s;
            A(l,i) = A(i,l);
        }
        R(i,k) = R(i,k)*c - R(i,l)*s;
        R(i,l) = R(i,l)*c + R(i,k)*s;
    }
        double maxval = max_offdiag_symmetric(A, k, l);
 }



 void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, const int maxiter, int& iterations, bool& converged)
 {
    int k=0; int l=1;
    double maxval = max_offdiag_symmetric(A, k, l);
    while (abs(A(k,l))> eps)
    {
        void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l);  
    }
    
 }
