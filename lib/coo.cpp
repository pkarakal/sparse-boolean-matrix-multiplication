#include "coo.hpp"
#include <algorithm>

namespace bmm_lib {
  template<class I, class T>
  void coo_to_csr(const I n_row, const I n_col, const I nnz, const I Ai[],
				  const I Aj[], const T Ax[], I Bp[], I Bj[], T Bx[]) {
	// compute number of non-zero entries per row of A
	std::fill(Bp, Bp + n_row, 0);

	for (I n = 0; n < nnz; n++) {
	  Bp[Ai[n]]++;
	}

	// cumsum the nnz per row to get Bp[]
	for (I i = 0, cumsum = 0; i < n_row; i++) {
	  I temp = Bp[i];
	  Bp[i] = cumsum;
	  cumsum += temp;
	}
	Bp[n_row] = nnz;

	// write Aj,Ax into Bj,Bx
	for (I n = 0; n < nnz; n++) {
	  I row = Ai[n];
	  I dest = Bp[row];

	  Bj[dest] = Aj[n];
	  Bx[dest] = Ax[n];

	  Bp[row]++;
	}

	for (I i = 0, last = 0; i <= n_row; i++) {
	  I temp = Bp[i];
	  Bp[i] = last;
	  last = temp;
	}

	// now Bp,Bj,Bx form a CSR representation (with possible duplicates)
  }

  template<class I, class T>
  void coo_to_csc(const I n_row, const I n_col, const I nnz, const I Ai[],
				  const I Aj[], const T Ax[], I Bp[], I Bi[], T Bx[]) {
	coo_to_csr<I, T>(n_col, n_row, nnz, Aj, Ai, Ax, Bp, Bi, Bx);
  }
} // namespace bmm_lib
