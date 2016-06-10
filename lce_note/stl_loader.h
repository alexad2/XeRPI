#include<vector>

#ifdef __CINT__
#pragma link C++ class vector<vector<float> >;
#else
template class std::vector<std::vector<float> >;
#endif
