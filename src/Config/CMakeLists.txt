IF ( NOT CF_HAVE_SINGLE_EXEC )
SET ( MYLIBNAME "Config" )
ELSE ()
SET ( MYLIBNAME "coolfluid-solver" ) 
ENDIF()

LIST (APPEND Config_files
BadMatchException.hh
BiggerThanZero.hh
BiggerThanZero.cxx
BuilderParser.cxx
BuilderParser.hh
BuilderParserFrameInfo.cxx
BuilderParserFrameInfo.hh
BuilderParserException.cxx
BuilderParserException.hh
BuilderParserMandType.hh
BuilderParserRules.cxx
BuilderParserRules.hh
Config.hh
ConfigArgs.cxx
ConfigArgs.hh
ConfigFacility.hh
ConfigFacility.cxx
ConfigFileReader.cxx
ConfigFileReader.hh
ConfigOptionException.hh
ConfigObject.cxx
ConfigObject.hh
ConfigOptionPtr.hh
ConverterTools.cxx
ConverterTools.hh
DuplicateNameException.hh
ManagerWorkerFrameType.hh
ManagerWorkerProtocol.hh
ManagerWorkerProtocol.cxx
NestedConfigFileReader.cxx
NestedConfigFileReader.hh
NoEnvVarException.hh
Option.cxx
Option.hh
OptionList.cxx
OptionList.hh
OptionMarkers.hh
OptionT.hh
OptionValidation.hh
OptionValidation.cxx
OptionValidationException.hh
PositiveLessThanOne.cxx
PositiveLessThanOne.hh
XMLConfigFileReader.cxx
XMLConfigFileReader.hh
)

IF ( NOT CF_HAVE_SINGLE_EXEC )
LIST ( APPEND Config_cflibs Common )
CF_ADD_KERNEL_LIBRARY ( Config )
ELSE()
FOREACH (AFILE ${Config_files} )
LIST(APPEND coolfluid-solver_files ../../src/Config/${AFILE} )
ENDFOREACH(AFILE)
CF_CACHE_LIST_APPEND ( ${MYLIBNAME}_files  ${coolfluid-solver_files} )
ENDIF()

CF_WARN_ORPHAN_FILES()
