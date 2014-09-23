#ifndef COOLFluiD_Framework_GET_DM_hh
#define COOLFluiD_Framework_GET_DM_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DomainModel.hh"
#include "Framework/VectorialFunction.hh"




namespace GETSpace {
class GETProxyBase;
}
namespace GET = GETSpace;


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace GETModel {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a model definition
 * described by analytical parametric functions.
 *
 * @author Tiago Quintino
 */
class GET_DM : public Framework::DomainModel {
public: // interface functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  GET_DM(const std::string& name);

  /// Destructor
  virtual ~GET_DM();

  /// Configure the data from the supplied arguments.
  /// @param args configuration arguments
  virtual void configure ( Config::ConfigArgs& args );

  /// gets the number of topological regions that have definition
  /// @return number of topological regions with definition
  virtual TRidx getNbTopoDefs () const ;

  /// computes the parametric coordinates given the real model coordinates
  /// @pre return parameter must be properly resized
  virtual void computeParamCoord (const TRidx idx, const XVector& coord , PVector& pcoord) const ;

  /// computes the coordinates in the model space given the parametric coordinates
  /// @pre return parameter must be properly resized
  virtual void computeCoord (const TRidx idx, const PVector& pcoord, XVector& coord) const ;

  /// computes the first derivatives in the model space given the parametric coordinates
  /// @pre return parameter must be properly resized
  virtual void compute1stDeriv (const TRidx idx, const PVector& pcoord, std::vector< XVector >& deriv1) const	;

  /// computes the second derivatives in the model space given the parametric coordinates
  /// @pre return parameter must be properly resized
  virtual void compute2ndDeriv (const TRidx idx, const PVector& pcoord, std::vector< XVector >& deriv2) const	;

  /// computes the coordinates, first and second derivatives all in one call
  /// @pre return parameters must be properly resized
	virtual void computeAll (const TRidx idx, const PVector& pcoord, XVector& coord, std::vector< XVector >& deriv1, std::vector< XVector >& deriv2) const ;

  /// Gets the class name
  static std::string getClassName() { return "GET_DM"; }

private: // data

	/// dimension of the model space
	CFuint mDim;

	/// name of the coresponding GET file
	std::string mFileName;

	/// pointer to the GET 
	GET::GETProxyBase	*mpGETProxy;

}; // end of class GET_DM

//////////////////////////////////////////////////////////////////////////////

  } // namespace GETModel

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_GET_DM_hh
