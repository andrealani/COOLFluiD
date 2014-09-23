#ifndef convert_to_string_hh
#define convert_to_string_hh

#include <sstream>

template <class T>
std::string convert_to_string(T v)
{
  std::ostringstream oss;
  oss << v ;
  std::string sstr = oss.str();
  return sstr;
}

#endif
