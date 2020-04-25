#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int main() {

  mat A(2, 4, fill::zeros);
  cout << A.n_rows << "\n";
  cout << A.n_cols << "\n";
  cout << "I got here!\n";

}
