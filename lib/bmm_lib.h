#ifndef SPARSE_BOOLEAN_MATRIX_MULTIPLICATION_BMM_LIB_H
#define SPARSE_BOOLEAN_MATRIX_MULTIPLICATION_BMM_LIB_H
#include "definitions.h"
#include <boost/program_options.hpp>
#include <iostream>
#include <ostream>
#include <string>
#include <vector>
#ifdef USE_MMIO_MATRICES
#include "mmio.h"
#endif

namespace bmm_lib {
  void parse_cli(int nargs, char** args, std::vector<std::string>& paths, std::vector<std::string>& dimensions);
#ifdef USE_MMIO_MATRICES
  void read_matrix(FILE* f, std::vector<uint32_t>& I, std::vector<uint32_t>& J,
#endif
} // namespace bmm_lib

#endif // SPARSE_BOOLEAN_MATRIX_MULTIPLICATION_BMM_LIB_H
