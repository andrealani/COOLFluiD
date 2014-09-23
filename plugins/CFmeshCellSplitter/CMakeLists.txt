LIST ( APPEND CFmeshCellSplitter_files
CFmeshCellSplitter.hh
CellSplitter2D.cxx
CellSplitter2D.hh
CellSplitter3D.cxx
CellSplitter3D.hh
CellSplitter2DFVM.cxx
CellSplitter2DFVM.hh
CellSplitter3DFVM.cxx
CellSplitter3DFVM.hh
)

LIST ( APPEND CFmeshCellSplitter_cflibs Framework )


CF_ADD_PLUGIN_LIBRARY ( CFmeshCellSplitter )
CF_WARN_ORPHAN_FILES()
