#IF (NOT CF_ENABLE_CUDA) 
SET ( CF_LIBRARY_LINK_FLAGS "" CACHE STRING "Extra link flags for libraries" FORCE )
#ELSE()
#SET ( CF_LIBRARY_LINK_FLAGS "-Xcompiler" CACHE STRING "Extra link flags for libraries" FORCE )
#ENDIF (NOT CF_ENABLE_CUDA)
MARK_AS_ADVANCED ( CF_LIBRARY_LINK_FLAGS )

########################################################################################
# UNIX
########################################################################################

IF(UNIX)
  
  # gnu specific warning flags
  IF( NOT CF_HAVE_CUDA AND CMAKE_COMPILER_IS_GNUCC )

    # use pipe for faster compilation
    CF_ADD_C_FLAGS  ("-pipe")
    CF_ADD_CXX_FLAGS("-pipe")
    # respect c 89 standard (same as -std=c89)
    CF_ADD_C_FLAGS("-ansi")

#IF(CF_ENABLE_MUTATIONPP)
    # use new features from c++0x
#    CF_ADD_CXX_FLAGS("-std=c++0x")
#ELSE()
    # respect c++ 98 standard
    CF_ADD_CXX_FLAGS("-std=c++98")
#ENDIF(CF_ENABLE_MUTATIONPP)
   
    # dont allow gnu extensions
    CF_ADD_CXX_FLAGS("-fno-gnu-keywords")
    
    # dont defined common variables
    CF_ADD_C_FLAGS("-fno-common")
    CF_ADD_CXX_FLAGS("-fno-common")

    if ( CF_ENABLE_WARNINGS )
      # use many warnings
      CF_ADD_CXX_FLAGS("-Wall")
     # this is added for checking 
#      CF_ADD_CXX_FLAGS("-W")
#      CF_ADD_CXX_FLAGS("-Wextra")
#      CF_ADD_CXX_FLAGS("-Woverloaded-virtual")
#      CF_ADD_CXX_FLAGS("-Wsign-promo")
#      CF_ADD_CXX_FLAGS("-Wformat")
      # Warn if an undefined identifier is evaluated in an #if directive.
      # AL: this flag activates many BOOST warnings
      #CF_ADD_CXX_FLAGS("-Wundef" )
      # accept functions that dont use all parameters, due to virtual functions may not need all
      CF_ADD_CXX_FLAGS("-Wno-unused-parameter")
      # remove boost many warnings
      CF_ADD_CXX_FLAGS("-Wno-missing-field-initializers")
      CF_ADD_CXX_FLAGS("-Wno-enum-compare")
      # this is temporary until we all move to using openmpi
      # must turn off non-virtual-dtor because many mpi implementations use it
      # KDE uses -Wnon-virtual-dtor
#      CF_ADD_CXX_FLAGS("-Wno-non-virtual-dtor")
      # must turn long long off because many mpi implementations use it
#      CF_ADD_CXX_FLAGS("-Wno-long-long")
      # be pedantic but issue warnings instead of errors
      # suppressing this option suppresses many BOOST warnings
      # CF_ADD_CXX_FLAGS("-pedantic")
#     CF_ADD_CXX_FLAGS("-fpermissive")

      # could add even these
      #-Wcast-align -Wchar-subscripts -Wpointer-arith -Wformat-security
    endif()

  ENDIF()

ENDIF(UNIX)

########################################################################################
# WINDOWS
########################################################################################

IF(WIN32)

  # stupid VS2005 warning about not using fopen
  ADD_DEFINITIONS( -D_CRT_SECURE_NO_DEPRECATE )
  # for M_PI in cmath
  ADD_DEFINITIONS( -D_USE_MATH_DEFINES )
  # disable auto-linking with boost
  ADD_DEFINITIONS( -DBOOST_ALL_NO_LIB )
  ADD_DEFINITIONS( -DBOOST_ALL_DYN_LINK )
  # Required for auto link not to mess up on vs80.
  # @todo Disable auto link on windows altogether.
  # ADD_DEFINITIONS( -DBOOST_DYN_LINK )

  # compilation flags
  #   /MD use the Multithreaded DLL of runtime library
  CF_ADD_C_FLAGS   ( "/MD" )
  CF_ADD_CXX_FLAGS ( "/MD" )

  # add exception handling
  CF_ADD_C_FLAGS   ( "/EHsc" )
  CF_ADD_CXX_FLAGS ( "/EHsc" )

  # under windows internal deps must be used so force them
  IF( NOT CF_ENABLE_INTERNAL_DEPS )
  	SET ( CF_ENABLE_INTERNAL_DEPS ON CACHE BOOL "Use of internal deps is forced" )
  ENDIF()

  # linker flags:
  #   /OPT:NOREF keeps functions and data that are never referenced ( needed for static libs )
  SET ( CF_LIBRARY_LINK_FLAGS "/OPT:NOREF /OPT:NOICF"  CACHE STRING "Extra link flags for libraries" FORCE )

  #   SET ( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /OPT:NOREF /OPT:NOICF" )
  #   SET ( CMAKE_CXX_CREATE_STATIC_LIBRARY  "lib ${CMAKE_CL_NOLOGO} /OPT:NOREF /OPT:NOICF <LINK_FLAGS> /out:<TARGET> <OBJECTS>" )

ENDIF(WIN32)

########################################################################################
# APPLE
########################################################################################

IF ( APPLE )

	# Apple linker with GCC does not support explicit templates so force them OFF
	IF( CMAKE_COMPILER_IS_GNUCC )
  		SET ( CF_HAVE_CXX_EXPLICIT_TEMPLATES OFF CACHE BOOL "Support for explicit templates deactivated -- Apple with GCC don't support it" )
	ENDIF()
	
	
	#  # improve the linker compiler to avoid unresolved symbols causing errors
        #  SET(CMAKE_CXX_CREATE_SHARED_LIBRARY
        #  "<CMAKE_CXX_COMPILER> -undefined dynamic_lookup <LANGUAGE_COMPILE_FLAGS> <CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS> <LINK_FLAGS> -o <TARGET> -install_name <TARGET_INSTALLNAME_DIR><TARGET_SONAME> <OBJECTS> <LINK_LIBRARIES>")

ENDIF()

########################################################################################
# GENERIC
########################################################################################

# Disable boost pre-1.34 boost::filesystem functions.
# ADD_DEFINITIONS ( -DBOOST_FILESYSTEM_NO_DEPRECATED )

########################################################################################
# FINAL
########################################################################################

# test and add the user defined flags

STRING ( REGEX MATCHALL "[^ ]+" C_FLAGS_LIST "${CF_C_FLAGS}"  )
FOREACH( c_flag ${C_FLAGS_LIST} )
  CF_ADD_C_FLAGS_SIGNAL_ERROR ( ${c_flag} )
ENDFOREACH()
MARK_AS_ADVANCED ( C_FLAGS_LIST   )

STRING ( REGEX MATCHALL "[^ ]+" CXX_FLAGS_LIST "${CF_CXX_FLAGS}"  )
FOREACH( cxx_flag ${CXX_FLAGS_LIST} )
  CF_ADD_CXX_FLAGS_SIGNAL_ERROR ( ${cxx_flag} )
ENDFOREACH()
MARK_AS_ADVANCED ( CXX_FLAGS_LIST  )

IF( NOT CF_SKIP_FORTRAN )
  STRING ( REGEX MATCHALL "[^ ]+" Fortran_FLAGS_LIST "${CF_Fortran_FLAGS}"  )
  # fortran flags currently nont checked
  SET ( CMAKE_Fortran_FLAGS "${CF_Fortran_FLAGS}" )
  # FOREACH( fortran_flag ${Fortran_FLAGS_LIST} )
  #   CF_ADD_Fortran_FLAGS_SIGNAL_ERROR ( ${fortran_flag} )
  # ENDFOREACH()
  MARK_AS_ADVANCED ( Fortran_FLAGS_LISTS )
ENDIF()

IF (CF_ENABLE_CUDA)
STRING ( REGEX MATCHALL "[^ ]+" CUDAC_FLAGS_LIST "${CF_CUDAC_FLAGS}"  )
SET ( CMAKE_CUDA_FLAGS "${CF_CUDAC_FLAGS}" )
MARK_AS_ADVANCED ( CUDAC_FLAGS_LIST )
ENDIF()


