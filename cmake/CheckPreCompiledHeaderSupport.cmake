# Check support for precompiled headers
IF(WIN32)

  # Assume WIN32 uses MSVC, therefore PCH support is ok
  set(PCHSUPPORT_FOUND TRUE)
  # define macro for using precompiled header
  macro(USE_PRECOMPILED_HEADER _project _hpp _cpp)
    add_executable(${_project}_precomp WIN32 ${_cpp} ${_hpp} )
    add_dependencies(${_project} ${_project}_precomp)
    set_target_properties(${_project}_precomp PROPERTIES COMPILE_FLAGS "/Yc${_hpp} /Fp\"$(IntDir)/${_project}.pch\"")
    # Flags to generate PCH
    set_target_properties(${_project} PROPERTIES COMPILE_FLAGS "/Yu${_hpp} /D${_project}_precomp_EXPORTS")
    # Flags to use PCH
  endmacro(USE_PRECOMPILED_HEADER _project _hdr _cpp)

ELSE()

  # if we have gcc > 3.4 then PCH support is ok
  # but for the moment we keep it out
  set(PCHSUPPORT_FOUND FALSE)

ENDIF()