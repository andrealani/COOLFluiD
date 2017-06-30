// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_MultiScalarTerm_hh
#define COOLFluiD_Framework_MultiScalarTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include <valarray>

#include <iostream>
#include "MathTools/MathFunctions.hh"
#include "Common/CFLog.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the interface for an Euler convective
/// physical term
/// @author Andrea Lani
template <typename BASE>
class MultiScalarTerm : public BASE {
public:

  /// Constructor without arguments
  MultiScalarTerm(const std::string& name) :
    BASE(name), _nbScalarVars(0), _nbScalarStrides(0), _firstVars(0)
  {
  }

  /// Default destructor
  virtual ~MultiScalarTerm()
  {
  }

  /// Set physical data
  virtual void setupPhysicalData()
  {
    cf_assert(getDataSize() > 0);
    
    this->resizePhysicalData(this->m_physicalData);
    this->resizePhysicalData(this->m_refPhysicalData);
  }
  
  /// Physical data size
  virtual CFuint getDataSize() const
  {
    return BASE::getDataSize() + MathTools::MathFunctions::innerProd(_nbScalarVars, _nbScalarStrides);
  }
  
  /// Get the start of the scalar vars data (mass fractions)
  virtual CFuint getFirstScalarVar(CFuint i) const
  {
    cf_assert(i < _firstVars.size());
    return _firstVars[i];
  }
  
  /// Get the number of scalar vars
  virtual CFuint getNbScalarVars(CFuint i) const
  {
    cf_assert(i < _nbScalarVars.size());
    return _nbScalarVars[i];
  }
  
  /// Get the number of sets of scalar vars
  virtual CFuint getNbScalarVarSets() const
  {
    cf_assert(_nbScalarVars.size() > 0);
    return _nbScalarVars.size();
  }
  
  /// Set the number of scalar vars
  virtual void setNbScalarVars(const std::valarray<CFuint>& nbScalarVars, 
			       std::valarray<CFuint>* nbScalarStrides = CFNULL)
  {
    _nbScalarVars.resize(nbScalarVars.size());
    _nbScalarVars = nbScalarVars;
    
    _nbScalarStrides.resize(nbScalarVars.size());
    _nbScalarStrides = 1;
    if (nbScalarStrides != CFNULL) {
      cf_assert(nbScalarVars.size() == nbScalarStrides->size());
      CFLog(DEBUG_MIN,  _nbScalarStrides.size() << " vs " << nbScalarStrides->size() << "\n");
      cf_assert(_nbScalarStrides.size() == nbScalarStrides->size());
      _nbScalarStrides = *nbScalarStrides; 
    }
    
    _firstVars.resize(nbScalarVars.size());
    
    CFuint counter = BASE::getDataSize();
    for (CFuint i = 0; i < _nbScalarVars.size(); ++i) {
      _firstVars[i] = counter;
      counter += _nbScalarVars[i];
    }
  }

private:

  /// number of scalar vars
  std::valarray<CFuint> _nbScalarVars;

  /// number of scalar strides
  std::valarray<CFuint> _nbScalarStrides;
  
  /// ID of the first var for each subset
  std::valarray<CFuint> _firstVars;
  
}; // end of class MultiScalarTerm

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_MultiScalarTerm_hh
