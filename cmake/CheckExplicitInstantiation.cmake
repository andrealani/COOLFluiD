# Test whether the compiler supports explicit template instantiation.
# This actually creates a class template instantiation in one source
# file and tries to use it from another.  This approach checks that
# both the instantiation syntax and symbol linkage is handled
# properly.

 MAKE_DIRECTORY(${COOLFluiD_BINARY_DIR}/CMakeTmp/CheckExplicitInstantiation)
 STRING(ASCII 35 POUND)
 WRITE_FILE(
   ${COOLFluiD_BINARY_DIR}/CMakeTmp/CheckExplicitInstantiation/CMakeLists.txt
   "PROJECT(EXPLICIT)\n"
   "ADD_LIBRARY(A A.cxx)\n"
   "ADD_EXECUTABLE(B B.cxx)\n"
   "TARGET_LINK_LIBRARIES(B A)\n"
   )
 WRITE_FILE(
   ${COOLFluiD_BINARY_DIR}/CMakeTmp/CheckExplicitInstantiation/A.h
   "${POUND}ifndef A_h\n"
   "${POUND}define A_h\n"
   "template <class T> class A { public: static T Method(); };\n"
   "${POUND}endif\n"
   )
 WRITE_FILE(
   ${COOLFluiD_BINARY_DIR}/CMakeTmp/CheckExplicitInstantiation/A.cxx
   "${POUND}include \"A.h\"\n"
   "template <class T> T A<T>::Method() { return 0; }\n"
   "template class A<int>;"
   )
 WRITE_FILE(
   ${COOLFluiD_BINARY_DIR}/CMakeTmp/CheckExplicitInstantiation/B.cxx
   "${POUND}include \"A.h\"\n"
   "int main() { return A<int>::Method(); }\n"
   )
 TRY_COMPILE( CF_CXX_SUPPORTS_EXPLICIT_TEMPLATES
   ${COOLFluiD_BINARY_DIR}/CMakeTmp/CheckExplicitInstantiation/Build
   ${COOLFluiD_BINARY_DIR}/CMakeTmp/CheckExplicitInstantiation
   EXPLICIT OUTPUT_VARIABLE OUTPUT
   )
 IF(CF_CXX_SUPPORTS_EXPLICIT_TEMPLATES)
   LOG( "+++++  Checking support for C++ explicit template instantiation -- yes")
   SET(CF_CXX_SUPPORTS_EXPLICIT_TEMPLATES ON CACHE INTERNAL "Support for C++ explict templates")
   WRITE_FILE(${CMAKE_BINARY_DIR}/CMakeFiles/CMakeOutput.log
     "Determining if the C++ compiler supports explict template instantiation passed with the following output:\n"
     "${OUTPUT}\n" APPEND)
 ELSE()
   LOG( "+++++  Checking support for C++ explicit template instantiation -- no")
   SET(CF_CXX_SUPPORTS_EXPLICIT_TEMPLATES OFF CACHE INTERNAL "Support for C++ explict templates")
   WRITE_FILE(${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log
     "Determining if the C++ compiler supports explict template instantiation failed with the following output:\n"
     "${OUTPUT}\n" APPEND)
 ENDIF()
