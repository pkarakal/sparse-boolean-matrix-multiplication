#include "../lib/bmm_lib.h"
#include "../lib/coo.hpp"
#include "../lib/definitions.h"
#include <algorithm>
#include <chrono>
#include <iostream>

//void multiply_matrices(bmm_lib::CSCMatrix& a, bmm_lib::CSCMatrix& b, bmm_lib::CSCMatrix& res){
//  std::vector<uint32_t> tmp{};
//  std::vector<uint32_t> resCols{};
//
//  for (uint32_t i = 0; i < a.row_count; i++) {
//	for (uint32_t j = a.cscRow[i]; j < a.cscRow[i + 1]; j++) {
//	  for (uint32_t k = b.cscRow[a.cscColumn[j]]; k < b.cscRow[a.cscColumn[j] + 1];
//		   k++) {
//		if (!std::binary_search(tmp.begin(), tmp.end(), b.cscColumn[k])) {
//		  tmp.push_back(b.cscColumn[k]);
//		  sort(tmp.begin(), tmp.end());
//		}
//	  }
//	}
//	resCols.insert(resCols.end(), tmp.begin(), tmp.end());
//	res.cscColumn[i + 1] = resCols.size();
//	tmp.clear();
//  }
//}

int main(int argc, char** argv) {
  try {
	std::vector<std::string> matrices{};
	std::vector<std::string> dimensions{};
	std::vector<bmm_lib::CSCMatrix> csc_matrices =
		std::vector<bmm_lib::CSCMatrix>(2);
	bmm_lib::parse_cli(argc, argv, matrices, dimensions);
#ifdef USE_MMIO_MATRICES
	std::vector<FILE*> input_matrices = std::vector<FILE*>(2);
	// check for capacity conflicts.
	if (matrices.size() != input_matrices.size()) {
	  std::cerr << "The number of matrices has to be equal to two."
				<< std::endl;
	  exit(-2);
	}
	// open input matrices in read mode only
	for (auto iter = matrices.begin(); iter != matrices.end(); ++iter) {
	  FILE* file = fopen(iter->c_str(), "r");
	  if (file == nullptr) {
		std::cerr << "File "
				  << matrices.at(std::distance(matrices.begin(), iter))
				  << " couldn't be opened or doesn't exist. Make sure the path "
					 "provided is correct"
				  << std::endl;
		exit(-1);
	  }
	  input_matrices.at(std::distance(matrices.begin(), iter)) = file;
	}
	matrices.clear();
	matrices.shrink_to_fit();

	for (auto iter = input_matrices.begin(); iter != input_matrices.end();
		 ++iter) {
	  std::vector<uint32_t> I{}, J{}, val{};
	  int rows{}, cols{};
	  bmm_lib::read_matrix(*iter, I, J, val, rows, cols);
	  // create symmetric values
	  std::vector<uint32_t> temp1 = std::vector<uint32_t>(I.begin(), I.end());
	  I.insert(std::end(I), std::begin(J), std::end(J));
	  J.insert(std::end(J), std::begin(temp1), std::end(temp1));
	  temp1.clear();
	  bmm_lib::CSCMatrix& cscMatrix =
		  csc_matrices.at(std::distance(input_matrices.begin(), iter));
	  cscMatrix.cscRow.resize(2 * I.size());
	  cscMatrix.cscColumn.resize(J.size() + 1);
	  cscMatrix.c_values.resize(2 * val.size());
	  cscMatrix.row_count = rows;
	  cscMatrix.col_count = cols;
	  if (I.at(0) < J.at(0)) {
		bmm_lib::coo_to_csc(uint32_t(I.size()), uint32_t(J.size()),
							uint32_t(2 * val.size()), I.data(), J.data(),
							val.data(), cscMatrix.cscColumn.data(),
							cscMatrix.cscRow.data(), cscMatrix.c_values.data());
	  } else {
		bmm_lib::coo_to_csc(uint32_t(J.size()), uint32_t(I.size()),
							uint32_t(2 * val.size()), J.data(), I.data(),
							val.data(), cscMatrix.cscColumn.data(),
							cscMatrix.cscRow.data(), cscMatrix.c_values.data());
	  }
	  cscMatrix.cscRow.shrink_to_fit();
	  cscMatrix.cscColumn.shrink_to_fit();
	  cscMatrix.c_values.shrink_to_fit();

	  I.clear();
	  J.clear();
	  val.clear();
	  I.shrink_to_fit();
	  J.shrink_to_fit();
	  val.shrink_to_fit();
	}
	bmm_lib::slice_matrix(csc_matrices.at(0), csc_matrices.at(1));
	std::cout << "After slice_matrix" << std::endl;
	for (const auto& matrix : csc_matrices) {
	  std::cout << matrix.cscColumn.size() << "\t" << matrix.cscRow.size()
				<< "\t" << matrix.c_values.size() << std::endl;
	}
#endif
#ifdef USE_RANDOM_MATRICES
	csc_matrices.at(0).row_count = std::atoi(dimensions.at(0).c_str());
	csc_matrices.at(0).col_count = std::atoi(dimensions.at(1).c_str());
	csc_matrices.at(1).row_count = std::atoi(dimensions.at(1).c_str());
	csc_matrices.at(1).col_count = std::atoi(dimensions.at(0).c_str());
	for (auto& matrix : csc_matrices) {
	  matrix.cscColumn = std::vector<uint32_t>(matrix.col_count);
	  matrix.cscRow = std::vector<uint32_t>(matrix.row_count);
	  matrix.c_values =
		  std::vector<uint32_t>(matrix.col_count * matrix.row_count);
	  std::vector<std::vector<uint32_t>> array =
		  std::vector<std::vector<uint32_t>>(matrix.row_count);
	  for (auto& col : array) {
		col = std::vector<uint32_t>(matrix.col_count);
	  }
	  bmm_lib::create_csc_matrices(array, matrix);
	  matrix.cscRow.shrink_to_fit();
	  matrix.cscColumn.shrink_to_fit();
	  matrix.c_values.shrink_to_fit();
	}
#endif
	bmm_lib::CSCMatrix result{};
	result.row_count = csc_matrices.at(0).row_count;
	result.col_count = csc_matrices.at(1).col_count;
	result.cscRow.resize(csc_matrices.at(0).row_count);
	result.cscColumn.resize(csc_matrices.at(1).col_count);
	result.c_values.resize(std::max(2 * csc_matrices.at(0).row_count +1,
									2 * (csc_matrices.at(1).col_count) + 1));
	auto start = std::chrono::high_resolution_clock::now();
	bmm_lib::csc_matmat(
		uint32_t(csc_matrices.at(0).row_count),
		uint32_t(csc_matrices.at(1).col_count),
		csc_matrices.at(0).cscRow.data(), csc_matrices.at(0).cscColumn.data(),
		csc_matrices.at(0).c_values.data(), csc_matrices.at(1).cscRow.data(),
		csc_matrices.at(1).cscColumn.data(), csc_matrices.at(1).c_values.data(),
		result.cscRow.data(), result.cscColumn.data(), result.c_values.data());
//	multiply_matrices(csc_matrices.at(0), csc_matrices.at(1), result);
//	for (auto& i : result.c_values) {
//	  i = 1;
//	}
	std::transform(result.c_values.begin(), result.c_values.end(), result.c_values.begin(), [](uint32_t i){ return 1;});
	auto stop = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> elapsed = stop - start;
	std::cout << "Took " << elapsed.count() * 10e-9 << "s" << std::endl;
	result.cscColumn.clear();
	result.cscRow.clear();
	result.c_values.clear();
	for (auto& m : csc_matrices) {
	  m.c_values.clear();
	  m.c_values.shrink_to_fit();
	  m.cscRow.clear();
	  m.cscRow.shrink_to_fit();
	  m.cscColumn.clear();
	  m.cscColumn.shrink_to_fit();
	}
	csc_matrices.clear();
	csc_matrices.shrink_to_fit();
  } catch (std::exception& ex) {
	std::cerr << ex.what() << std::endl;
  }
  return 0;
}
