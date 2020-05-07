/* plasmaState.i */ 
%module pyplasmastate
%include "stl.i"
%include "std_string.i"
%include "std_vector.i"
%include "cpointer.i"
%apply int *OUTPUT { int* icomp};
// Instantiate templates used by example
namespace std {
   %template(IntVector) vector<int>;
   %template(DoubleVector) vector<double>;
   %template(CharVector) vector<char>;
   %template(StringVector) vector<string>;
   %template(SizeTVector) vector<size_t>;
}
%{
#include "cxxps.h"
#include "multiarray.h"
#include <iostream>
%}
%include "cxxps.h"
%include "multiarray.h"
%apply int *OUTPUT { int* icomp};
%extend MultiArray<double> {
   const double& __getitem__(std::vector <size_t> iv){
      return ((*($self))[iv]);
   }
}
%extend MultiArray<int> {
   const int& __getitem__(std::vector <size_t> iv){
      return ((*($self))[iv]);
   }
}
%extend MultiArray<char> {
   const char& __getitem__(std::vector <size_t> iv){
      return ((*($self))[iv]);
   }
}
%extend MultiArray<string> {
   const char& __getitem__(std::vector <size_t> iv){
      return ((*($self))[iv]);
   }
}
