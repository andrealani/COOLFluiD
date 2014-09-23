##############################################################################
# kernel library macro
##############################################################################
MACRO( CF_ADD_KERNEL_LIBRARY LIBNAME )

  # declare this library as part of the kernel
  SET ( ${LIBNAME}_kernellib ON )
  
  CF_ADD_LIBRARY( ${LIBNAME} )

ENDMACRO( CF_ADD_KERNEL_LIBRARY )
##############################################################################


