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
						  std::vector<uint32_t>& J, std::vector<uint32_t>& val,
						  int& row_count, int& col_count) {
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

  I.resize(nnz);
  J.resize(nnz);
  val.resize(nnz);

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
  row_count = M;
  col_count = N;
}
#endif


int bmm_lib::fRand(int min, int max){
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(0, 1);
  int f = (int) dist(mt);
  return min + f * (max - min);
}


void bmm_lib::create_csc_matrices(std::vector<std::vector<uint32_t>>& matrix, bmm_lib::CSCMatrix& mat) {
  std::random_device rnd_device;
  // Specify the engine and distribution.
  std::mt19937 mersenne_engine {rnd_device()};  // Generates random integers
  std::uniform_int_distribution<int> dist {0, 1};

  auto gen = [&dist, &mersenne_engine](){
	return dist.min() + dist(mersenne_engine)* (dist.max()- dist.min());};

  long non_zero_count = 0;
  for(auto& row: matrix){
	std::generate(row.begin(), row.end(), gen);
	std::transform(row.begin(), row.end(), row.begin(), [&](uint32_t i){
	  if (i>0){
		++non_zero_count;
	  }
	  return i;
	});
  }
  mat.c_values.resize(non_zero_count);


  convert_dense_to_csc(matrix, mat);
}

/**
 * This takes two csc matrices and slices them up. As the operation is
 * matrix * ref_matrix matrix_col has to be equal to ref_matrix_row.
 * To achieve this, we compare them. If they don't match we have to remove
 * some elements from the vectors
 *
 * If matrix_col > ref_matrix_row, it searches for column entries greater than
 * ref_matrix_row and removes them completely. As csc formatted matrices are
 * column major, it is as easy as just looping until we find such a column. In
 * each loop the the row_position is incremented by the difference between current
 * index value and previous index value. Then every value after iter in cscColumn
 * and every value after row_pos in cscRow and c_values can be safely deleted from
 * the corresponding vectors
 *
 * If matrix_col < ref_matrix_row, rows from ref_matrix have to be deleted. To
 * achieve it, there is a loop over cscColumn, and in each loop we search between
 * the two columns there is value in a row greater than matrix_col. If there is
 * we loop over each such value and erase it from the corresponding vector. As a
 * result we have to decrease all the following values from the cscColumn items.
 * Because of some unwanted memory bugs, I ended up just swapping the matrices in
 * memory and calling the function again.
 *
 * Time complexity O(n^2)
 *
 * @param {CSCMatrix&} matrix
 * @param {CSCMatrix&} ref_matrix
 */
void bmm_lib::slice_matrix(bmm_lib::CSCMatrix& matrix,
						   bmm_lib::CSCMatrix& ref_matrix) {
  if (matrix.col_count == ref_matrix.row_count) {
	return;
  }
  if (matrix.col_count > ref_matrix.row_count) {
	long row_pos = 0;
	auto previous = matrix.cscColumn.begin();
	for(auto iter = ++(matrix.cscColumn.begin()); iter != matrix.cscColumn.end(); previous=iter, ++iter){
	  row_pos += *iter - *previous;
	  if (*iter > ref_matrix.row_count){
		matrix.cscColumn.erase(iter, matrix.cscColumn.end());
		matrix.cscRow.erase(matrix.cscRow.begin()+row_pos+1, matrix.cscRow.end());
		matrix.c_values.erase(matrix.c_values.begin()+row_pos+1, matrix.c_values.end());
		break;
	  }
	}
  } else {
	std::swap(matrix, ref_matrix);
	bmm_lib::slice_matrix(matrix, ref_matrix);
//	long row_pos_start;
//	long row_pos_end;
//	for(auto iter = ++(ref_matrix.cscColumn.begin()); iter != ref_matrix.cscColumn.end(); ++iter){
//	  row_pos_start = *std::prev(iter);
//	  row_pos_end = *std::next(iter);
//	  if(row_pos_start >= ref_matrix.cscRow.size() || row_pos_end >= ref_matrix.cscRow.size()){
//		break;
//	  }
//	  for(auto it= ref_matrix.cscRow.begin()+row_pos_start; it != ref_matrix.cscRow.begin()+row_pos_end; ++it ){
//		if(*it > matrix.col_count){
//		  long pos = std::distance(matrix.cscRow.begin(), matrix.cscRow.begin() + std::distance(ref_matrix.cscRow.begin(), it));
//		  if (pos < 0 || pos >= ref_matrix.cscRow.size()) {
//			break;
//		  }
//		  ref_matrix.cscRow.erase(ref_matrix.cscRow.begin() + pos);
//		  ref_matrix.c_values.erase(ref_matrix.c_values.begin() + pos);
//		  std::transform(iter, ref_matrix.cscColumn.end(), iter, [](uint32_t i){ return i > 0 ? --i: i;});
//		  --row_pos_end;
//		  if(row_pos_end == row_pos_start){
//			break;
//		  }
//		}
//	  }
////	  row_pos_start = *iter - *previous + 1;
//
////	  auto it = std::find_if(ref_matrix.cscRow.begin()+row_pos_start, ref_matrix.cscRow.begin()+row_pos_end, [&](uint32_t i) {return i > matrix.col_count;});
////	  while (it != matrix.cscRow.begin()+row_pos_end || it!= matrix.cscRow.begin() || it!=matrix.cscRow.end()){
////		std::transform(iter, ref_matrix.cscColumn.end(), iter, [](uint32_t i){ return --i;});
////		long pos = std::distance(matrix.cscRow.begin(), matrix.cscRow.begin() + std::distance(ref_matrix.cscRow.begin(), it));
////		if(pos < 0){
////		  break;
////		}
////		ref_matrix.cscRow.erase(ref_matrix.cscRow.begin()+ pos);
////		ref_matrix.c_values.erase(ref_matrix.c_values.begin()+ pos);
////		--row_pos_end;
////		it = std::find_if(ref_matrix.cscRow.begin()+row_pos_start, ref_matrix.cscRow.begin()+row_pos_end, [&](uint32_t i) {return i > matrix.col_count;});
////	  }
//	}
  }
}

void bmm_lib::convert_dense_to_csc(std::vector<std::vector<uint32_t>>& dense, bmm_lib::CSCMatrix& sparse){
  for(int i=0; i < dense.size(); ++i){
	for(int j=0; j < dense.at(i).size(); ++j){
	  if(dense.at(i).at(j) > 0){
		sparse.cscRow.push_back(i);
		sparse.cscColumn.push_back(j);
		sparse.c_values.push_back(dense.at(i).at(j));
	  }
	}
  }
}
