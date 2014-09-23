// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <iostream> 

#include "Common/ProcessInfoWin32.hh"
#include "Common/Common.hh"

#include <windows.h> // for CaptureStackBacktrace
#include <dbghelp.h> // for stack trace
#include <psapi.h>   // for memory usage   

// The arraysize(arr) macro returns the # of elements in an array arr.
// The expression is a compile-time constant, and therefore can be
// used in defining new arrays, for example.  If you use arraysize on
// a pointer by mistake, you will get a compile-time error.
//
// One caveat is that arraysize() doesn't accept any array of an
// anonymous type or a type defined inside a function.  In these rare
// cases, you have to use the unsafe ARRAYSIZE_UNSAFE() macro below.  This is
// due to a limitation in C++'s template system.  The limitation might
// eventually be removed, but it hasn't happened yet.

// This template function declaration is used in defining arraysize.
// Note that the function doesn't need an implementation, as we only
// use its type.
template <typename T, size_t N>
char (&ArraySizeHelper(T (&array)[N]))[N];

// That gcc wants both of these prototypes seems mysterious. VC, for
// its part, can't decide which to use (another mystery). Matching of
// template overloads: the final frontier.
#ifndef _MSC_VER
template <typename T, size_t N>
char (&ArraySizeHelper(const T (&array)[N]))[N];
#endif

#define arraysize(array) (sizeof(ArraySizeHelper(array)))

struct _EXCEPTION_POINTERS;

// A stacktrace can be helpful in debugging. For example, you can include a
// stacktrace member in a object (probably around #ifndef NDEBUG) so that you
// can later see where the given object was created from.
class StackTrace {
 public:
  // Creates a stacktrace from the current location
  StackTrace();

  // Creates a stacktrace for an exception.
  // Note: this function will throw an import not found (StackWalk64) exception
  // on system without dbghelp 5.1.
  StackTrace(_EXCEPTION_POINTERS* exception_pointers);

  // Gets an array of instruction pointer values.
  //   count: (output) the number of elements in the returned array
  const void *const *Addresses(size_t* count);
  // Prints a backtrace to stderr
  void PrintBacktrace();

  // Resolves backtrace to symbols and write to stream.
  void OutputToStream(std::ostream* os);

 private:
  // From http://msdn.microsoft.com/en-us/library/bb204633.aspx,
  // the sum of FramesToSkip and FramesToCapture must be less than 63,
  // so set it to 62. Even if on POSIX it could be a larger value, it usually
  // doesn't give much more information.
  static const int MAX_TRACES = 62;
  void* trace_[MAX_TRACES];
  int count_;

};

namespace {

// SymbolContext is a threadsafe singleton that wraps the DbgHelp Sym* family
// of functions.  The Sym* family of functions may only be invoked by one
// thread at a time.
class SymbolContext {
 public:
  static SymbolContext* Get() 
  {
    
    // We use a leaky singleton because code may call this during process termination.
    static SymbolContext aSymContxt;
    return &aSymContxt;
  }

  // Returns the error code of a failed initialization.
  DWORD init_error() const {
    return init_error_;
  }

  // For the given trace, attempts to resolve the symbols, and output a trace
  // to the ostream os.  The format for each line of the backtrace is:
  //
  //    <tab>SymbolName[0xAddress+Offset] (FileName:LineNo)
  //
  // This function should only be called if Init() has been called.  We do not
  // LOG(FATAL) here because this code is called might be triggered by a
  // LOG(FATAL) itself.
  void OutputTraceToStream(const void* const* trace,
                           int count,
                           std::ostream* os) {
    
    /* AutoLock lock(lock_); */ // from chromium

    for (int i = 0; (i < count) && os->good(); ++i) {
      const int kMaxNameLength = 256;
      DWORD_PTR frame = reinterpret_cast<DWORD_PTR>(trace[i]);

      // Code adapted from MSDN example:
      // http://msdn.microsoft.com/en-us/library/ms680578(VS.85).aspx
      ULONG64 buffer[
        (sizeof(SYMBOL_INFO) +
          kMaxNameLength * sizeof(wchar_t) +
          sizeof(ULONG64) - 1) /
        sizeof(ULONG64)];

      // Initialize symbol information retrieval structures.
      DWORD64 sym_displacement = 0;
      PSYMBOL_INFO symbol = reinterpret_cast<PSYMBOL_INFO>(&buffer[0]);
      symbol->SizeOfStruct = sizeof(SYMBOL_INFO);
      symbol->MaxNameLen = kMaxNameLength;
      BOOL has_symbol = SymFromAddr(GetCurrentProcess(), frame,
                                    &sym_displacement, symbol);

      // Attempt to retrieve line number information.
      DWORD line_displacement = 0;
      IMAGEHLP_LINE64 line = {};
      line.SizeOfStruct = sizeof(IMAGEHLP_LINE64);
      BOOL has_line = SymGetLineFromAddr64(GetCurrentProcess(), frame,
                                           &line_displacement, &line);

      // Output the backtrace line.
      (*os) << "\t";
      if (has_symbol) {
        (*os) << symbol->Name << " [0x" << trace[i] << "+"
              << sym_displacement << "]";
      } else {
        // If there is no symbol informtion, add a spacer.
        (*os) << "(No symbol) [0x" << trace[i] << "]";
      }
      if (has_line) {
        (*os) << " (" << line.FileName << ":" << line.LineNumber << ")";
      }
      (*os) << "\n";
    }
  }

 private:

  SymbolContext() : init_error_(ERROR_SUCCESS) {
    // Initializes the symbols for the process.
    // Defer symbol load until they're needed, use undecorated names, and
    // get line numbers.
    SymSetOptions(SYMOPT_DEFERRED_LOADS |
                  SYMOPT_UNDNAME |
                  SYMOPT_LOAD_LINES);
    if (SymInitialize(GetCurrentProcess(), NULL, TRUE)) {
      init_error_ = ERROR_SUCCESS;
    } else {
      init_error_ = GetLastError();
      // TODO(awong): Handle error: SymInitialize can fail with
      // ERROR_INVALID_PARAMETER.
      // When it fails, we should not call debugbreak since it kills the current
      // process (prevents future tests from running or kills the browser
      // process).
      std::cerr << "SymInitialize failed: " << init_error_ << std::endl;
    }
  }

  DWORD init_error_;
};

}  // namespace

StackTrace::StackTrace() {
  // When walking our own stack, use CaptureStackBackTrace().
  count_ = CaptureStackBackTrace(0, arraysize(trace_), trace_, NULL);
}

StackTrace::StackTrace(EXCEPTION_POINTERS* exception_pointers) {
  // When walking an exception stack, we need to use StackWalk64().
  count_ = 0;
  // Initialize stack walking.
  STACKFRAME64 stack_frame;
  memset(&stack_frame, 0, sizeof(stack_frame));
#if defined(_WIN64)
  int machine_type = IMAGE_FILE_MACHINE_AMD64;
  stack_frame.AddrPC.Offset = exception_pointers->ContextRecord->Rip;
  stack_frame.AddrFrame.Offset = exception_pointers->ContextRecord->Rbp;
  stack_frame.AddrStack.Offset = exception_pointers->ContextRecord->Rsp;
#else
  int machine_type = IMAGE_FILE_MACHINE_I386;
  stack_frame.AddrPC.Offset = exception_pointers->ContextRecord->Eip;
  stack_frame.AddrFrame.Offset = exception_pointers->ContextRecord->Ebp;
  stack_frame.AddrStack.Offset = exception_pointers->ContextRecord->Esp;
#endif
  stack_frame.AddrPC.Mode = AddrModeFlat;
  stack_frame.AddrFrame.Mode = AddrModeFlat;
  stack_frame.AddrStack.Mode = AddrModeFlat;
  while (StackWalk64(machine_type,
                     GetCurrentProcess(),
                     GetCurrentThread(),
                     &stack_frame,
                     exception_pointers->ContextRecord,
                     NULL,
                     &SymFunctionTableAccess64,
                     &SymGetModuleBase64,
                     NULL) &&
         count_ < arraysize(trace_)) {
    trace_[count_++] = reinterpret_cast<void*>(stack_frame.AddrPC.Offset);
  }
}

void StackTrace::PrintBacktrace() {
  OutputToStream(&std::cerr);
}

void StackTrace::OutputToStream(std::ostream* os) {
  SymbolContext* context = SymbolContext::Get();
  DWORD error = context->init_error();
  if (error != ERROR_SUCCESS) {
    (*os) << "Error initializing symbols (" << error
          << ").  Dumping unresolved backtrace:\n";
    for (int i = 0; (i < count_) && os->good(); ++i) {
      (*os) << "\t" << trace_[i] << "\n";
    }
  } else {
    (*os) << "Backtrace:\n";
    context->OutputTraceToStream(trace_, count_, os);
  }
}


//////////////////////////////////////////////////////////////////////////////

using namespace std;

namespace COOLFluiD {
  namespace Common {

//////////////////////////////////////////////////////////////////////////////

ProcessInfoWin32::ProcessInfoWin32()
{
}

//////////////////////////////////////////////////////////////////////////////

ProcessInfoWin32::~ProcessInfoWin32()
{
}

//////////////////////////////////////////////////////////////////////////////

std::string ProcessInfoWin32::getBackTrace () const
{
  printf ("\n\nWin32 dumping backtrace ...\n");

  std::ostringstream oss;

  StackTrace trace;
  trace.OutputToStream( &oss );


#if 0

  oss << "No backtace implemented in Win32" << endl;

#endif

#if 0
  const int max_callers = 62;

	void *array[max_callers];
  printf ("Calling  CaptureStackBackTrace ...\n");
  int num = CaptureStackBackTrace(0,max_callers,array, NULL);

  printf ("Returned %d frames ...\n", num);
  for (int i = 0; i < num; i++)
  {
    printf("%s\n",(char*) array[i]);
		oss << (char*) array[i] << "\n";
  }
#endif

#if 0

  // From http://msdn.microsoft.com/en-us/library/bb204633(VS.85).aspx,
  // the sum of FramesToSkip and FramesToCapture must be less than 63,
  // so set it to 62.
  const int kMaxCallers = 62;

  void* callers[kMaxCallers];
  // TODO(ajwong): Migrate this to StackWalk64.
  int count = CaptureStackBackTrace(0, kMaxCallers, callers, NULL);
  if (count > 0) {
    trace_.resize(count);
    memcpy(&trace_[0], callers, sizeof(callers[0]) * count);
  } else {
    trace_.resize(0);
  // When walking our own stack, use CaptureStackBackTrace().
  count_ = CaptureStackBackTrace(0, arraysize(trace_), trace_, NULL);  
  
#endif


  return oss.str();
}

//////////////////////////////////////////////////////////////////////////////

CFuint ProcessInfoWin32::getPID () const
{
  return (CFuint) GetCurrentProcessId();
}

//////////////////////////////////////////////////////////////////////////////

CFdouble ProcessInfoWin32::memoryUsageBytes () const
{
  CFdouble return_value = 0.;  

#if 1
  HANDLE hProcess = GetCurrentProcess();

  PROCESS_MEMORY_COUNTERS pmc;

  if ( hProcess != NULL )
  {
    if ( GetProcessMemoryInfo( hProcess, &pmc, sizeof(pmc)) )
    {
//        printf( "\tPageFaultCount: 0x%08X\n", pmc.PageFaultCount );
//        printf( "\tPeakWorkingSetSize: 0x%08X\n", pmc.PeakWorkingSetSize );
//        printf( "\tWorkingSetSize: 0x%08X\n", pmc.WorkingSetSize );
//        printf( "\tQuotaPeakPagedPoolUsage: 0x%08X\n", pmc.QuotaPeakPagedPoolUsage );
//        printf( "\tQuotaPagedPoolUsage: 0x%08X\n", pmc.QuotaPagedPoolUsage );
//        printf( "\tQuotaPeakNonPagedPoolUsage: 0x%08X\n", pmc.QuotaPeakNonPagedPoolUsage );
//        printf( "\tQuotaNonPagedPoolUsage: 0x%08X\n",  pmc.QuotaNonPagedPoolUsage );
//        printf( "\tPagefileUsage: 0x%08X\n",      pmc.PagefileUsage ); 
//        printf( "\tPeakPagefileUsage: 0x%08X\n",  pmc.PeakPagefileUsage );
    }

    return_value = (CFuint) pmc.WorkingSetSize;
    
  CloseHandle( hProcess );
  }
#endif

  return return_value;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common
} // namespace COOLFluiD

