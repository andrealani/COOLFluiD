LIST ( APPEND RemeshMeandr_files
AdaptSubSystem.cxx
AdaptSubSystem.hh
CallMeandros.cxx
CallMeandros.hh
RMeshMe.cxx
RMeshMe.hh
RMeshMeData.cxx
RMeshMeData.hh
RemeshMeandr.hh
StdSetup.cxx
StdSetup.hh
StdUnSetup.cxx
StdUnSetup.hh
TAdptMaestro.cxx
TAdptMaestro.hh
WriteMetricCom.cxx
WriteMetricCom.hh
)

LIST ( APPEND RemeshMeandr_cflibs Framework )

CF_ADD_PLUGIN_LIBRARY ( RemeshMeandr )
CF_WARN_ORPHAN_FILES()
