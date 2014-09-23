// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_LinearAdv_ADTerm_hh
#define COOLFluiD_Physics_LinearAdv_ADTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"

#include "LinearAdv/LinearAdv.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
namespace LinearAdv {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the interface for a NavierStokesModel.
/// @author Nadege Villedieu
class LinearAdv_API ADTerm : public Framework::BaseTerm
{
public:
  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Enumerator defining the mapping between
  /// the variable name and its position in the
  /// physical data
  enum {NU=0};

  /// Constructor without arguments
  ADTerm(const std::string& name);

  /// Default destructor
  virtual ~ADTerm();

  /// Physical data size
  CFuint getDataSize() const {  return 1; }

  /// Configures this object by complementing the
  /// implementation in ConfigObject
  virtual void configure ( Config::ConfigArgs& args );

  /// Set physical data
  virtual void setupPhysicalData();

  /// Get the name
  static std::string getName()  {  return "ADTerm"; }

  /// Get the name
  CFreal getDiffusionCoeff() const {   return m_DiffCoef;   }

private:
  /// Diffusion coefficient
  CFreal m_DiffCoef;

}; // end of class ADTerm

//////////////////////////////////////////////////////////////////////////////

} // namespace LinearAdv
  }
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LinearAdv_ADTerm_hh
