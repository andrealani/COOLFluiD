#ifndef COOLFluiD_Physics_MHD_MHD2DProjectionPrimE_hh
#define COOLFluiD_Physics_MHD_MHD2DProjectionPrimE_hh

//////////////////////////////////////////////////////////////////////////////

#include "MHD/MHD2DProjectionPrim.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a MHD physical model 2D for projection scheme
 * for primitive variables
 *
 * @author Andrea Lani
 * @author Radka Keslerova
 */

		class MHD2DProjectionPrimE : public MHD2DProjectionPrim {
		public: //function

			/**
			* Constructor
			* @see MHD2DProjection
			*/
			MHD2DProjectionPrimE(Common::SafePtr<Framework::BaseTerm> term);

			/**
			* Default destructor
			*/
			virtual ~MHD2DProjectionPrimE();

			/**
			* Set up the private data and give the maximum size of states physical
			* data to store
			*/
			virtual void setup();

			/**
			* Get extra variable names
			*/
			virtual std::vector<std::string> getExtraVarNames() const;

			/**
			* Set the jacobians
			*/
			virtual void computeJacobians();

			/**
			* Split the jacobian
			*/
			virtual void splitJacobian(RealMatrix& jacobPlus,
				RealMatrix& jacobMin,
				RealVector& eValues,
				const RealVector& normal);

			/**
			* Set the matrix of the right and left eigenvectors
			* and the matrix of the eigenvalues
			*/
			virtual void computeEigenValuesVectors(RealMatrix& rightEv,
				RealMatrix& leftEv,
				RealVector& eValues,
				const RealVector& normal);

			/**
			* Set the PhysicalData corresponding to the given State
			* @see EulerPhysicalModel
			*/
			virtual void computePhysicalData(const Framework::State& state,
				RealVector& data);

			/**
			* Set the total magnetic field and energy values
			*/
			virtual void setDimensionalValuesPlusExtraValues(const Framework::State& state,
				RealVector& result,
				RealVector& extra);

			/**
			* Set a State starting from the given PhysicalData
			* @see EulerPhysicalModel
			*/
			virtual void computeStateFromPhysicalData(const RealVector& data,
				Framework::State& state);

			/// Set the vector of the eigenValues
			virtual void computeEigenValues(const RealVector& pdata,
				const RealVector& normal,
				RealVector& eValues);

			/// Get the maximum eigenvalue
			virtual CFreal getMaxEigenValue(const RealVector& pdata, const RealVector& normal);

			/// Get the maximum absolute eigenvalue
			virtual CFreal getMaxAbsEigenValue(const RealVector& pdata, const RealVector& normal);

		protected:

			/// Computes the convective flux projected on a normal
			virtual void computeFlux(const RealVector& pdata, const RealVector& normals);

			/// Computes the physical convective flux
			virtual void computeStateFlux(const RealVector& pdata);

}; // end of class MHD2DProjectionPrimE

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif //COOLFluiD_Physics_MHD_MHD2DProjectionPrim_hh
