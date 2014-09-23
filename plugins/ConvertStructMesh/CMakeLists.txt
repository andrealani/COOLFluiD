LIST ( APPEND ConvertStructMesh_files
Block.hh
ConvertBlockMesh.hh
ConvertBlockMesh.cxx
ConvertGridProMesh.cxx
ConvertGridProMesh.hh
ConvertQuadMesh.cxx
ConvertQuadMesh.hh
ConvertStructMesh.cxx
ConvertStructMesh.hh
ConvertStructMeshApp.hh
FaceData.cxx
FaceData.hh
NodeData.cxx
NodeData.hh
main.cxx
)

LIST ( APPEND ConvertStructMesh_libs ${CF_KERNEL_LIBS}  ${CF_Boost_LIBRARIES}  )

CF_ADD_PLUGIN_APP ( ConvertStructMesh )
CF_WARN_ORPHAN_FILES()
