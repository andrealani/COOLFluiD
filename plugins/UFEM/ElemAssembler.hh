#ifndef COOLFluiD_UFEM_ElemAssembler_hh
#define COOLFluiD_UFEM_ElemAssembler_hh

//////////////////////////////////////////////////////////////////////////////

#include "UFEM/UFEMTerm.hh"
#include "UFEM/ElemProps.hh"
#include "UFEM/AssemblyData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD { namespace Framework { class GeometricEntity; } }

namespace COOLFluiD {
    namespace UFEM {

//////////////////////////////////////////////////////////////////////////////

struct UFEM_API ElemAssembler
{
  typedef std::vector < Common::SelfRegistPtr<UFEMTerm> > VecUFEMTerm;

  explicit ElemAssembler ( const CFuint size ) : terms(), eprops (), assdata(size) {}

  /// element assembler owns the UFEMTerms
  VecUFEMTerm terms;
  /// element assembler owns the ElemProps
  Common::SelfRegistPtr<ElemProps> eprops;
  /// element assembler owns the AssemblyData
  AssemblyData  assdata; // ;)
};

//////////////////////////////////////////////////////////////////////////////

  } // namespace UFEM
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_UFEM_ElemAssembler_hh

