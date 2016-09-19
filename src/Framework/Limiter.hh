// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_Limiter_hh
#define COOLFluiD_Framework_Limiter_hh

//////////////////////////////////////////////////////////////////////////////

#include "Config/ConfigObject.hh"
#include "BaseMethodStrategyProvider.hh"
#include "Framework/Node.hh"
#include "Common/OwnedObject.hh"
#include "BaseDataSocketSink.hh"
#include "MethodStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Framework {

      class GeometricEntity;

//////////////////////////////////////////////////////////////////////////////

/// This class offers a basic interface for all limiters for FVM
/// @author Andrea Lani
template < typename METHODDATA >
class Limiter : public Framework::MethodStrategy<METHODDATA> {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  typedef Framework::BaseMethodStrategyProvider<METHODDATA,Limiter<METHODDATA> > PROVIDER;

  /// Constructor
  Limiter(const std::string& name);

  /// Default destructor
  virtual ~Limiter();

  /// Set private data that will be used during the computation
  virtual void setup() 
  {
    Framework::MethodStrategy<METHODDATA>::setup(); 
  }
  
  /// Unsetup private data that will be used during the computation
  virtual void unsetup() 
  {
    Framework::MethodStrategy<METHODDATA>::unsetup(); 
  }
  
  /// Compute the flux in the current face
  /// @post in limiterValue you put the value of the limiter for
  ///       each variable
  virtual void limit(const std::vector<std::vector<Node*> >& coord,
      Framework::GeometricEntity* const cell,
      CFreal* limiterValue) = 0;

  /// Apply the face limiter
  virtual void limitOnFace(const RealVector& rLeft,
  		   const RealVector& rRight,
  		   CFreal* limiterValue)
  {
    throw Common::NotImplementedException (FromHere(),"Limiter::limitOnFace()");
  }

  /// Apply the limiter to a scalar quantity
  virtual void limitScalar(CFreal r, CFreal& limiterValue)
  {
    throw Common::NotImplementedException (FromHere(),"Limiter::limitScalar()");
  }

  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args )
  {
    Config::ConfigObject::configure(args);
  }

  /// Returns the DataSocket's that this numerical strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets()
  {
    std::vector<Common::SafePtr<BaseDataSocketSink> > result;

    return result;
  }

  /// Gets the Class name
  static std::string getClassName()
  {
    return "Limiter";
  }
  
  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}
  
protected:

  /// Private Copy Constructor
  Limiter(const Limiter& o);

  /// Private Assignement operator
  const Limiter& operator=(const Limiter& o);
  
  /// parameter to mitigate the control the limiter value <= 1
  CFreal m_alpha;
  
  /// flag telling if to use the full stencil to compute the local extrema
  bool m_useFullStencil;
  
  /// flag telling if to use the stencil used for the nodal extrapolation to compute the local extrema
  bool m_useNodalExtrapolationStencil;
  
}; // end of class Limiter
      
//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "Limiter.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_Limiter_hh
