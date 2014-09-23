LIST ( APPEND FluctSplitMeta_files
BaseSplitter.hh
BaseSubElem.hh
ContourIntegral.hh
ElemGeo.hh
FSData.hh
FSMHO.cxx
FSMHO.hh
FluctSplitMeta.hh
IRD.hh
LRD.hh
NschemeT.hh
NullIntegral.hh
Physics.hh
TriagP1.hh
TriagP2.hh
TetraP1.hh
TetraP2.hh
CIntegralTriagP1O1.hh
CIntegralTriagP2O1.hh
CIntegralTriagP2O2.hh
CIntegralTriagP2O3.hh
CIntegralTetraP1O1.hh
CIntegralTetraP2O1.hh
ShapeFunctions.hh
Euler2D.hh
Euler3D.hh
VolumeIntegral.hh
)

LIST ( APPEND FluctSplitMeta_cflibs FluctSplit )

CF_ADD_PLUGIN_LIBRARY ( FluctSplitMeta )

CF_WARN_ORPHAN_FILES()
