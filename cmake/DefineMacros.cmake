##############################################################################
# include cmake macros
##############################################################################
INCLUDE(CheckIncludeFile)
INCLUDE(CheckIncludeFileCXX)
INCLUDE(CheckIncludeFiles)
INCLUDE(CheckSymbolExists)
INCLUDE(CheckFunctionExists)
INCLUDE(CheckLibraryExists)
INCLUDE(CheckTypeSize)
INCLUDE(CheckCSourceCompiles)
INCLUDE(CheckCXXSourceCompiles)
INCLUDE(CheckCCompilerFlag)
INCLUDE(CheckCXXCompilerFlag)


##############################################################################
# include coolfluid macros
##############################################################################

INCLUDE(macros/CFVariables)
INCLUDE(macros/CFListOperations)
INCLUDE(macros/CFSearchPaths)
INCLUDE(macros/CFOptionAddSubdirectory)
INCLUDE(macros/CFSetUnion)
INCLUDE(macros/CFLogToFile)
INCLUDE(macros/CFDumpVariables)
INCLUDE(macros/CFBoolTo01)
INCLUDE(macros/CFSeparateSources)
INCLUDE(macros/CFAddLibrary)
INCLUDE(macros/CFAddKernelLibrary)
INCLUDE(macros/CFAddPluginLibrary)
INCLUDE(macros/CFAddPluginApp)
INCLUDE(macros/CFAddTest)
INCLUDE(macros/CFWarnOrphanFiles)
INCLUDE(macros/CFCheckFileLength)
INCLUDE(macros/CFAddCompilationFlags)
INCLUDE(macros/CFAddTestCase)
