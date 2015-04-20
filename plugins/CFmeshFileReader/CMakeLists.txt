LIST(APPEND  CFmeshFileReader_files
CFmeshReader.cxx
ReadBase.hh
ReadBase.cxx
ReadCFmesh.cxx
StdSetup.cxx
StdUnSetup.cxx
CFmeshReaderData.cxx
CFmeshFileReader.hh
CFmeshFileReaderAPI.hh
CFmeshReader.hh
CFmeshReaderData.hh
ReadCFmesh.hh
ReadDummy.cxx
ReadDummy.hh
StdSetup.hh
StdUnSetup.hh
)

##############################################################################
# Option on MPI
##############################################################################
# TODO: this dependency on MPI should somehow be removed

LIST ( APPEND OPTIONAL_dirfiles
       ParReadCFmesh.hh
       ParReadCFmesh.ci
       ParReadCFmesh.cxx
       ParCFmeshBinaryFileReader.hh 
       ParCFmeshBinaryFileReader.cxx
       ParCFmeshFileReader.hh 
       ParCFmeshFileReader.cxx
     )

IF ( CF_HAVE_MPI ) 
  LIST ( APPEND CFmeshFileReader_files 
		ParReadCFmesh.hh 
		ParReadCFmesh.ci
		ParReadCFmesh.cxx
		ParCFmeshBinaryFileReader.hh 
		ParCFmeshBinaryFileReader.cxx
		ParCFmeshFileReader.hh 
		ParCFmeshFileReader.cxx
  )
  
  IF ( CF_HAVE_PARMETIS )
     LIST ( APPEND CFmeshFileReader_includedirs ${PARMETIS_INCLUDE_DIR} )
     LIST ( APPEND CFmeshFileReader_libs ${PARMETIS_LIBRARIES} )
  ENDIF ()

ENDIF()

##############################################################################

LIST ( APPEND CFmeshFileReader_cflibs Framework ShapeFunctions )

CF_ADD_PLUGIN_LIBRARY ( CFmeshFileReader )

CF_WARN_ORPHAN_FILES()
