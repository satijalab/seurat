#include <Rinternals.h>

// helper to determine if external c++ pointer is valid
SEXP isnull(SEXP pointer) {
  return Rf_ScalarLogical(!R_ExternalPtrAddr(pointer));
}