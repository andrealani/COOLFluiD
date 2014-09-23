###############################################################################
# Basic files
LIST ( APPEND Common_files
ArrayAllocator.hh
BigAllocator.hh
CFAssert.cxx
CFAssert.hh
CodeLocation.cxx
CodeLocation.hh
Compatibility.hh
Common.hh
CommonAPI.hh
DebugFunctions.hh
ExportAPI.hh
EnumT.hh
COOLFluiD.hh
CppTypeInfo.hh
DemangledTypeID.cxx
DemangledTypeID.hh
DynamicObject.hh
DynamicObject.cxx
Exception.cxx
Exception.hh
FailedAssertionException.hh
FloatingPointException.hh
Fortran.hh
GlobalReduce.hh
MemoryAllocator.hh
MemoryAllocatorNormal.cxx
MemoryAllocatorNormal.hh
NonCopyable.hh
NonInstantiable.hh
NotImplementedException.hh
NullPointerException.hh
NullPointerException.cxx
OSystem.hh
OSystem.cxx
Obj_Helper.hh
ParallelException.hh
ParallelException.cxx
PE.cxx
PE.hh
PEFunctions.hh
PEInterface.cxx
PEInterface.hh
PM.hh
PtrAlloc.hh
ProcessInfo.hh
ProcessInfo.cxx
StlHeaders.hh
ShouldNotBeHereException.hh
ShouldNotBeHereException.cxx
Array2D.hh
BadValueException.hh
CFLog.hh
CFLog.cxx
CFMap2D.ci
CFMap2D.hh
CFMap3D.hh
CFMap.ci
CFMap.hh
CFMultiMap.ci
CFMultiMap.hh
CFPrintContainer.hh
SharedPtr.hh
ConnectivityTable.hh
DynamicFunctionCaller.hh
DynamicFunctionCaller.cxx
EventHandler.cxx
EventHandler.hh
FailedCastException.hh
FastMap.ci
FastMap.hh
FileDownloader.cxx
FileDownloader.hh
FileFormatException.hh
FilesystemException.hh
FilesystemException.cxx
GeneralStorage.hh
GrowArray.hh
HourMinSec.cxx
HourMinSec.hh
IOHelpers.hh
IsFundamental.hh
LibLoader.hh
LibLoader.cxx
LookupTable2D.ci
LookupTable2D.hh
MemFunArg.hh
NamedObject.cxx
NamedObject.hh
NullableObject.hh
NullableObject.cxx
NoSuchStorageException.hh
NoSuchValueException.hh
NoSuchValueException.cxx
OldLookupTable2D.ci
OldLookupTable2D.hh
OldLookupTable.ci
OldLookupTable.hh
OwnedObject.hh
ParserException.hh
ParserException.cxx
ProviderBase.cxx
ProviderBase.hh
Quartet.hh
SafePtr.hh
SelfRegistPtr.hh
SetupException.hh
SetupObject.cxx
SetupObject.hh
SignalHandler.hh
SignalHandler.cxx
WorkerStatus.hh
WorkerStatus.cxx
Stopwatch.ci
Stopwatch.hh
StorageExistsException.hh
SwapEmpty.hh
StringOps.hh
StringOps.cxx
Table.ci
Table.hh
TaggedObject.cxx
TaggedObject.hh
TimePolicies.cxx
TimePolicies.hh
Trio.hh
URLException.hh
VarRegistry.cxx
VarRegistry.hh
xmlParser.h
xmlParser.cpp
)

###############################################################################
# template meta-programming classes
LIST ( APPEND Common_files Meta/Loop.hh Meta/Power.hh )

###############################################################################

# MPI Files
IF (CF_HAVE_MPI)
  LIST(APPEND Common_includedirs ${MPI_INCLUDE_DIR} )
  LIST(APPEND Common_files
    MPI/CommIteratorBase.hh
    MPI/CommPattern.cxx
    MPI/CommPattern.hh
    MPI/GlobalReduceMPI.hh
    MPI/MPICommPattern.hh
    MPI/MPIError.hh
    MPI/MPIDataType.cxx
    MPI/MPIDataTypeHandler.cxx
    MPI/MPIDataTypeHandler.hh
    MPI/MPIDataType.hh
    MPI/MPIDataTypeRegistrar_Helper.cxx
    MPI/MPIDataTypeRegistrar_Helper.hh
    MPI/MPIDataTypeRegistrar.hh
    MPI/MPIException.hh
    MPI/MPIHelper.cxx
    MPI/MPIHelper.hh
    MPI/MPIInitObject.hh
    MPI/MPIStructDef.hh
    MPI/ParVectorException.hh
    MPI/ParVector.hh
    MPI/PatternFullExchange.hh
    MPI/PatternRing.hh
    MPI/PEInterfaceMPI.cxx
    MPI/PEInterfaceMPI.hh
  )
LIST(APPEND Common_libs ${MPI_LIBRARIES})
ENDIF()

###############################################################################
# SERIAL PE
IF(NOT CF_HAVE_MPI)
  LIST(APPEND Common_files
    SERIAL/GlobalReduceSERIAL.hh
    SERIAL/PEInterfaceSERIAL.hh
    SERIAL/PEInterfaceSERIAL.cxx
  )
ENDIF()

###############################################################################

IF (CF_HAVE_CUDA)
  LIST(APPEND Common_includedirs ${CUDA_INCLUDE_DIR} )
  LIST(APPEND Common_files
    CUDA/ArrayT.hh
    CUDA/CFMat.hh
    CUDA/CFMatSlice.hh
    CUDA/CFVec.hh
    CUDA/CFVecSlice.hh
    CUDA/CudaDeviceManager.cu
    CUDA/CudaDeviceManager.hh
    CUDA/CudaEnv.hh
    CUDA/CudaLock.hh
    CUDA/CudaTimer.hh
    CUDA/CudaVector.hh
    CUDA/ExprT.hh
    CUDA/MacrosET.hh
    CUDA/MatExprT.hh
    CUDA/MathFunctions.hh
  )                       
  LIST(APPEND Common_libs ${CUDA_LIBRARIES}) 
ENDIF ()

###############################################################################

# Curl download helper
# LIST ( APPEND OPTIONAL_dirfiles CurlDownloader.hh CurlDownloader.cxx)
IF(CF_HAVE_CURL)
  LIST(APPEND Common_includedirs ${CURL_INCLUDE_DIR} )
  LIST(APPEND Common_files CurlDownloader.hh CurlDownloader.cxx )
  LIST(APPEND Common_libs ${CURL_LIBRARY})
ENDIF()


###############################################################################
# MMap allocation
LIST ( APPEND OPTIONAL_dirfiles MemoryAllocatorMMap.hh MemoryAllocatorMMap.cxx )
IF(CF_HAVE_ALLOC_MMAP AND CF_HAVE_UNISTD_H )
  LIST(APPEND Common_files MemoryAllocatorMMap.hh MemoryAllocatorMMap.cxx )
ENDIF()

###############################################################################
# Operating System dependent files
LIST ( APPEND OPTIONAL_dirfiles PosixDlopenLibLoader.hh PosixDlopenLibLoader.cxx  )
LIST ( APPEND OPTIONAL_dirfiles Win32LibLoader.hh Win32LibLoader.cxx )

LIST ( APPEND OPTIONAL_dirfiles SignalHandlerLinux.hh SignalHandlerLinux.cxx )
LIST ( APPEND OPTIONAL_dirfiles SignalHandlerMacOSX.hh SignalHandlerMacOSX.cxx )
LIST ( APPEND OPTIONAL_dirfiles SignalHandlerWin32.hh SignalHandlerWin32.cxx )

LIST ( APPEND OPTIONAL_dirfiles ProcessInfoLinux.hh  ProcessInfoLinux.cxx )
LIST ( APPEND OPTIONAL_dirfiles ProcessInfoMacOSX.hh ProcessInfoMacOSX.cxx )
LIST ( APPEND OPTIONAL_dirfiles ProcessInfoWin32.hh ProcessInfoWin32.cxx )

### UNIXES ####
IF (UNIX)

  IF ( CF_HAVE_DLOPEN )
    LIST ( APPEND Common_files PosixDlopenLibLoader.hh PosixDlopenLibLoader.cxx )
    LIST ( APPEND Common_libs  dl )
    SET ( CF_HAVE_LIBLOADER ON )
  ENDIF()

  IF    ( CF_OS_MACOSX )
    LIST ( APPEND Common_files SignalHandlerMacOSX.hh SignalHandlerMacOSX.cxx )
    LIST ( APPEND Common_files ProcessInfoMacOSX.hh ProcessInfoMacOSX.cxx )
  ENDIF ()

  IF    ( CF_OS_LINUX )
    LIST ( APPEND Common_files ProcessInfoLinux.hh ProcessInfoLinux.cxx )
    LIST ( APPEND Common_files SignalHandlerLinux.hh SignalHandlerLinux.cxx )
  ENDIF ()

ENDIF(UNIX)

### WINDOWS ####
IF (WIN32)


  LIST ( APPEND Common_files SignalHandlerWin32.hh SignalHandlerWin32.cxx )
  LIST ( APPEND Common_files ProcessInfoWin32.hh ProcessInfoWin32.cxx )

  # should search explicitly for dbghelp.dll
  # should search for psapi.dll if not in windows 7 (already in kernel32.dll)

  # libraries for ProcessInfoWin32
  LIST ( APPEND Common_libs dbghelp psapi )

  IF(CF_HAVE_WINDOWSH)
    LIST ( APPEND Common_files Win32LibLoader.hh Win32LibLoader.cxx)
    SET ( CF_HAVE_LIBLOADER ON )
  ENDIF(CF_HAVE_WINDOWSH)

ENDIF()

IF    ( NOT CF_HAVE_LIBLOADER )
  MESSAGE ( FATAL_ERROR "No dynamic library loader was found!")
ENDIF()

###############################################################################

LIST ( APPEND Common_cflibs logcpp )

CF_ADD_KERNEL_LIBRARY ( Common )

CF_WARN_ORPHAN_FILES()
