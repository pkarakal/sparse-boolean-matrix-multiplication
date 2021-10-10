#include "bmm_lib.h"

namespace po = boost::program_options;
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
  copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
  return os;
}

void bmm_lib::parse_cli(int nargs, char** args,
						std::vector<std::string>& paths,
						std::vector<std::string>& dimensions) {
  po::options_description desc("Allowed options");
  desc.add_options()
	  ("help", "produce help message")
#ifdef USE_MMIO_MATRICES
	  ("path,p", po::value<std::vector<std::string>>(), "input path to file")
#endif
#ifdef USE_RANDOM_MATRICES
	  ("dimensions,d", po::value<std::vector<std::string>>(), "dimensions of matrices to generate")
#endif
	  ;
  po::variables_map vm;
  po::store(po::parse_command_line(nargs, args, desc), vm);
  po::notify(vm);
  if (vm.count("help")) {
	std::cout << desc << std::endl;
	return;
  }
#ifdef USE_MMIO_MATRICES
  paths = vm.count("path") ? vm["path"].as<std::vector<std::string>>()
						   : std::vector<std::string>();
#endif
#ifdef USE_RANDOM_MATRICES
  dimensions = vm.count("dimensions") ? vm["dimensions"].as<std::vector<std::string>>(): std::vector<std::string>();
#endif
}
#ifdef USE_MMIO_MATRICES
void bmm_lib::read_matrix(FILE* f, std::vector<uint32_t>& I,
  MM_typecode matcode;
  int M, N, nnz;
  uint32_t i;

  if (mm_read_banner(f, &matcode) != 0) {
	printf("Could not process Matrix Market banner.\n");
	exit(1);
  }

  /*  This is how one can screen matrix types if their application */
  /*  only supports a subset of the Matrix Market data types.      */

  if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
	  mm_is_sparse(matcode)) {
	printf("Sorry, this application does not support ");
	printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
	exit(1);
  }

  /* find out size of sparse matrix .... */

  if (mm_read_mtx_crd_size(f, &M, &N, &nnz) != 0)
	exit(1);

  I = std::vector<uint32_t>(nnz);
  J = std::vector<uint32_t>(nnz);
  val = std::vector<uint32_t>(nnz);

  /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
  /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
  /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

  /* Replace missing val column with 1s and change the fscanf to match pattern
   * matrices*/

  if (!mm_is_pattern(matcode)) {
	for (i = 0; i < nnz; i++) {
	  fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
	  I[i]--; /* adjust from 1-based to 0-based */
	  J[i]--;
	}
  } else {
	for (i = 0; i < nnz; i++) {
	  fscanf(f, "%d %d\n", &I[i], &J[i]);
	  val[i] = 1;
	  I[i]--; /* adjust from 1-based to 0-based */
	  J[i]--;
	}
  }

  if (f != stdin)
	fclose(f);

  if (M != N) {
	printf("COO matrix' columns and rows are not the same");
  }
}
#endif