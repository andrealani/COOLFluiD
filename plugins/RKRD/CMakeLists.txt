LIST ( APPEND RKRD_files
RKRD.cxx
RKRD.hh
RKRDData.cxx
RKRDData.hh
RungeKuttaRD.hh
StdSetup.cxx
StdSetup.hh
StdUnSetup.cxx
StdUnSetup.hh
Backup.hh
Backup.cxx
Update.hpp
Update.cpp
Shift.hpp
Shift.cpp
# strategy
RKRDStrategy.hh
RKRDStrategy.cxx
)

LIST ( APPEND FluctSplitLinEuler_requires_mods FluctSplit )
LIST ( APPEND RKRD_cflibs FluctSplit )

CF_ADD_PLUGIN_LIBRARY ( RKRD )
CF_WARN_ORPHAN_FILES()
