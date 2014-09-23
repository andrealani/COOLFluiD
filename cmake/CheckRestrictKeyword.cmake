# This module checks if the C++ compiler supports the restrict keyword or
# some variant of it. The following variants are checked for in that order:
# 1. restrict            (The standard C99 keyword, not yet in C++ standard, Windows VS has it)
# 2. __restrict          (G++ has it)
# 3. __restrict__        (G++ has it too)
# 4. _Restrict           (seems to be used by Sun's compiler)
# These four cases seem to cover all existing variants; however some C++
# compilers don't support any variant, in which case the CF_RESTRICT_KEYWORD variable is set to nothing

SET(_CHECK_restrict_KEYWORD_SRC "
char f( const char * restrict x ){  return *x;}
int main(int argc, char *argv[]) { return 0; }
")

SET(_CHECK___restrict_KEYWORD_SRC "
char f( const char * __restrict x ){  return *x;}
int main(int argc, char *argv[]) { return 0; }
")

SET(_CHECK___restrict___KEYWORD_SRC "
char f( const char * __restrict__ x ){  return *x;}
int main(int argc, char *argv[]) { return 0; }
")

SET(_CHECK__Restrict_KEYWORD_SRC "
char f( const char * _Restrict x ) {  return *x; }
int main(int argc, char *argv[]) { return 0; }
")

CHECK_CXX_SOURCE_COMPILES("${_CHECK_restrict_KEYWORD_SRC}"     HAVE_KEYWORD_restrict)
IF(HAVE_KEYWORD_restrict)
  SET(CF_RESTRICT_KEYWORD restrict)
ELSE()
  CHECK_CXX_SOURCE_COMPILES("${_CHECK___restrict_KEYWORD_SRC}"   HAVE_KEYWORD___restrict)
  IF(HAVE_KEYWORD___restrict)
    SET(CF_RESTRICT_KEYWORD __restrict)
  ELSE()
    CHECK_CXX_SOURCE_COMPILES("${_CHECK___restrict___KEYWORD_SRC}" HAVE_KEYWORD___restrict__)
    IF(HAVE_KEYWORD___restrict__)
      SET(CF_RESTRICT_KEYWORD __restrict__)
    ELSE()
      CHECK_CXX_SOURCE_COMPILES("${_CHECK__Restrict_KEYWORD_SRC}"    HAVE_KEYWORD__Restrict)
      IF(HAVE_KEYWORD__Restrict)
        SET(CF_RESTRICT_KEYWORD _Restrict)
      ELSE()
        SET(CF_RESTRICT_KEYWORD) # not supported so keep it empty
      ENDIF()
    ENDIF()
  ENDIF()    
ENDIF()


  