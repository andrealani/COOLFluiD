// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShapeFunctions/LagrangeShapeFunctionHexaP2.hh"
#include "MathTools/MatrixInverterT.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



   namespace ShapeFunctions {

//////////////////////////////////////////////////////////////////////////////

CFuint LagrangeShapeFunctionHexaP2::_interpolatorID = 0;
RealMatrix LagrangeShapeFunctionHexaP2::m_invJ(3,3);
RealVector LagrangeShapeFunctionHexaP2::_vec1 = RealVector(DIM_3D);
RealVector LagrangeShapeFunctionHexaP2::_vec2 = RealVector(DIM_3D);
RealVector LagrangeShapeFunctionHexaP2::_vec3 = RealVector(DIM_3D);
std::vector<RealVector> LagrangeShapeFunctionHexaP2::_normals = std::vector<RealVector>(6);

//////////////////////////////////////////////////////////////////////////////

void LagrangeShapeFunctionHexaP2::computeGradientStates(
         const std::vector<RealMatrix>& jacob,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& grad)
{
    MathTools::MatrixInverterT<3> inverter;

    //taken from book: "The Finite Element Method Displayed" by G. Dhatt and G. Touzot, Wiley&Sons, 1984
    for (CFuint ip = 0; ip < mappedCoord.size(); ++ip) {
// CFout << "Jacobian : " << jacob[ip] << "\n";
// CFout << "Jacobian Determinant: " <<       jacob[ip].determ3() <<"\n";
      inverter.invert(jacob[ip],m_invJ);
      RealMatrix& lgrad = grad[ip];

      const CFreal xi =  mappedCoord[ip][KSI];
      const CFreal eta = mappedCoord[ip][ETA];
      const CFreal zeta = mappedCoord[ip][ZTA];
      const CFreal xi2  = xi*xi;
      const CFreal eta2 = eta*eta;
      const CFreal zeta2= zeta*zeta;
// CFout << "Mapped Coord: " <<  mappedCoord[ip] <<"\n";

      //Corner nodes
      CFreal xi0 = -1.;
      CFreal eta0 = -1.;
      CFreal zeta0 = -1.;
      const CFreal dN0dxi  = 0.125 * xi0  * (1. + eta * eta0)*(1.+zeta*zeta0)*(-1.+ 2.*xi*xi0 + eta*eta0 + zeta*zeta0);
      const CFreal dN0deta = 0.125 * eta0 * (1. + xi * xi0)*(1.+zeta*zeta0)*(-1.+ xi*xi0 + 2.*eta*eta0 + zeta*zeta0);
      const CFreal dN0dzta = 0.125 * zeta0 * (1. + xi * xi0)*(1.+eta*eta0)*(-1.+ xi*xi0 + eta*eta0 + 2.*zeta*zeta0);

      xi0   = 1.;
      eta0  = -1.;
      zeta0 = -1.;
      const CFreal dN1dxi  = 0.125 * xi0  * (1. + eta * eta0)*(1.+zeta*zeta0)*(-1.+ 2.*xi*xi0 + eta*eta0 + zeta*zeta0);
      const CFreal dN1deta = 0.125 * eta0 * (1. + xi * xi0)*(1.+zeta*zeta0)*(-1.+ xi*xi0 + 2.*eta*eta0 + zeta*zeta0);
      const CFreal dN1dzta = 0.125 * zeta0 * (1. + xi * xi0)*(1.+eta*eta0)*(-1.+ xi*xi0 + eta*eta0 + 2.*zeta*zeta0);

      xi0   = 1.;
      eta0  = 1.;
      zeta0 = -1.;
      const CFreal dN2dxi  = 0.125 * xi0  * (1. + eta * eta0)*(1.+zeta*zeta0)*(-1.+ 2.*xi*xi0 + eta*eta0 + zeta*zeta0);
      const CFreal dN2deta = 0.125 * eta0 * (1. + xi * xi0)*(1.+zeta*zeta0)*(-1.+ xi*xi0 + 2.*eta*eta0 + zeta*zeta0);
      const CFreal dN2dzta = 0.125 * zeta0 * (1. + xi * xi0)*(1.+eta*eta0)*(-1.+ xi*xi0 + eta*eta0 + 2.*zeta*zeta0);

      xi0   = -1.;
      eta0  = 1.;
      zeta0 = -1.;
      const CFreal dN3dxi  = 0.125 * xi0  * (1. + eta * eta0)*(1.+zeta*zeta0)*(-1.+ 2.*xi*xi0 + eta*eta0 + zeta*zeta0);
      const CFreal dN3deta = 0.125 * eta0 * (1. + xi * xi0)*(1.+zeta*zeta0)*(-1.+ xi*xi0 + 2.*eta*eta0 + zeta*zeta0);
      const CFreal dN3dzta = 0.125 * zeta0 * (1. + xi * xi0)*(1.+eta*eta0)*(-1.+ xi*xi0 + eta*eta0 + 2.*zeta*zeta0);

      xi0   = -1.;
      eta0  = -1.;
      zeta0 = 1.;
      const CFreal dN4dxi  = 0.125 * xi0  * (1. + eta * eta0)*(1.+zeta*zeta0)*(-1.+ 2.*xi*xi0 + eta*eta0 + zeta*zeta0);
      const CFreal dN4deta = 0.125 * eta0 * (1. + xi * xi0)*(1.+zeta*zeta0)*(-1.+ xi*xi0 + 2.*eta*eta0 + zeta*zeta0);
      const CFreal dN4dzta = 0.125 * zeta0 * (1. + xi * xi0)*(1.+eta*eta0)*(-1.+ xi*xi0 + eta*eta0 + 2.*zeta*zeta0);

      xi0   = 1.;
      eta0  = -1.;
      zeta0 = 1.;
      const CFreal dN5dxi  = 0.125 * xi0  * (1. + eta * eta0)*(1.+zeta*zeta0)*(-1.+ 2.*xi*xi0 + eta*eta0 + zeta*zeta0);
      const CFreal dN5deta = 0.125 * eta0 * (1. + xi * xi0)*(1.+zeta*zeta0)*(-1.+ xi*xi0 + 2.*eta*eta0 + zeta*zeta0);
      const CFreal dN5dzta = 0.125 * zeta0 * (1. + xi * xi0)*(1.+eta*eta0)*(-1.+ xi*xi0 + eta*eta0 + 2.*zeta*zeta0);

      xi0   = 1.;
      eta0  = 1.;
      zeta0 = 1.;
      const CFreal dN6dxi  = 0.125 * xi0  * (1. + eta * eta0)*(1.+zeta*zeta0)*(-1.+ 2.*xi*xi0 + eta*eta0 + zeta*zeta0);
      const CFreal dN6deta = 0.125 * eta0 * (1. + xi * xi0)*(1.+zeta*zeta0)*(-1.+ xi*xi0 + 2.*eta*eta0 + zeta*zeta0);
      const CFreal dN6dzta = 0.125 * zeta0 * (1. + xi * xi0)*(1.+eta*eta0)*(-1.+ xi*xi0 + eta*eta0 + 2.*zeta*zeta0);

      xi0   = -1.;
      eta0  = 1.;
      zeta0 = 1.;
      const CFreal dN7dxi  = 0.125 * xi0  * (1. + eta * eta0)*(1.+zeta*zeta0)*(-1.+ 2.*xi*xi0 + eta*eta0 + zeta*zeta0);
      const CFreal dN7deta = 0.125 * eta0 * (1. + xi * xi0)*(1.+zeta*zeta0)*(-1.+ xi*xi0 + 2.*eta*eta0 + zeta*zeta0);
      const CFreal dN7dzta = 0.125 * zeta0 * (1. + xi * xi0)*(1.+eta*eta0)*(-1.+ xi*xi0 + eta*eta0 + 2.*zeta*zeta0);

      //Nodes on the side parallel to axis xi
      xi0   = 0.;
      eta0  = -1.;
      zeta0 = -1.;

      const CFreal dN8dxi  = -0.5 * xi * (1. + eta*eta0)*(1. + zeta*zeta0);
      const CFreal dN8deta = 0.25 * eta0 * (1. - xi2)*(1. + zeta*zeta0);
      const CFreal dN8dzta = 0.25 * zeta0 * (1. - xi2)*(1. + eta*eta0);

      eta0  = 1.;
      zeta0 = -1.;
      const CFreal dN10dxi  = -0.5 * xi * (1. + eta*eta0)*(1. + zeta*zeta0);
      const CFreal dN10deta = 0.25 * eta0 * (1. - xi2)*(1. + zeta*zeta0);
      const CFreal dN10dzta = 0.25 * zeta0 * (1. - xi2)*(1. + eta*eta0);

      eta0  = -1.;
      zeta0 = 1.;
      const CFreal dN16dxi  = -0.5 * xi * (1. + eta*eta0)*(1. + zeta*zeta0);
      const CFreal dN16deta = 0.25 * eta0 * (1. - xi2)*(1. + zeta*zeta0);
      const CFreal dN16dzta = 0.25 * zeta0 * (1. - xi2)*(1. + eta*eta0);

      eta0  = 1.;
      zeta0 = 1.;
      const CFreal dN18dxi  = -0.5 * xi * (1. + eta*eta0)*(1. + zeta*zeta0);
      const CFreal dN18deta = 0.25 * eta0 * (1. - xi2)*(1. + zeta*zeta0);
      const CFreal dN18dzta = 0.25 * zeta0 * (1. - xi2)*(1. + eta*eta0);

      //Nodes on the side parallel to axis eta
      eta0  = 0.;

      xi0   = 1.;
      zeta0 = -1.;
      const CFreal dN9dxi  = 0.25 * xi0 * (1. - eta2)*(1. + zeta*zeta0);
      const CFreal dN9deta = -0.5 * eta * (1. + xi*xi0)*(1. + zeta*zeta0);
      const CFreal dN9dzta = 0.25 * zeta0 * (1. - eta2)*(1. + xi*xi0);

      xi0   = -1.;
      zeta0 = -1.;
      const CFreal dN11dxi  = 0.25 * xi0 * (1. - eta2)*(1. + zeta*zeta0);
      const CFreal dN11deta = -0.5 * eta * (1. + xi*xi0)*(1. + zeta*zeta0);
      const CFreal dN11dzta = 0.25 * zeta0 * (1. - eta2)*(1. + xi*xi0);

      xi0   = 1.;
      zeta0 = 1.;
      const CFreal dN17dxi  = 0.25 * xi0 * (1. - eta2)*(1. + zeta*zeta0);
      const CFreal dN17deta = -0.5 * eta * (1. + xi*xi0)*(1. + zeta*zeta0);
      const CFreal dN17dzta = 0.25 * zeta0 * (1. - eta2)*(1. + xi*xi0);

      xi0   = -1.;
      zeta0 = 1.;
      const CFreal dN19dxi  = 0.25 * xi0 * (1. - eta2)*(1. + zeta*zeta0);
      const CFreal dN19deta = -0.5 * eta * (1. + xi*xi0)*(1. + zeta*zeta0);
      const CFreal dN19dzta = 0.25 * zeta0 * (1. - eta2)*(1. + xi*xi0);

      //Nodes on the side parallel to axis eta
      zeta0  = 0.;

      xi0   = -1.;
      eta0  = -1.;
      const CFreal dN12dxi  = 0.25 * xi0 * (1. - zeta2)*(1. + eta*eta0);
      const CFreal dN12deta = 0.25 * eta0 * (1. - zeta2)*(1. + xi*xi0);
      const CFreal dN12dzta = -0.5 * zeta * (1. + xi*xi0)*(1. + eta*eta0);

      xi0   = 1.;
      eta0  = -1.;
      const CFreal dN13dxi  = 0.25 * xi0 * (1. - zeta2)*(1. + eta*eta0);
      const CFreal dN13deta = 0.25 * eta0 * (1. - zeta2)*(1. + xi*xi0);
      const CFreal dN13dzta = -0.5 * zeta * (1. + xi*xi0)*(1. + eta*eta0);

      xi0   = 1.;
      eta0  = 1.;
      const CFreal dN14dxi  = 0.25 * xi0 * (1. - zeta2)*(1. + eta*eta0);
      const CFreal dN14deta = 0.25 * eta0 * (1. - zeta2)*(1. + xi*xi0);
      const CFreal dN14dzta = -0.5 * zeta * (1. + xi*xi0)*(1. + eta*eta0);

      xi0   = -1.;
      eta0  = 1.;
      const CFreal dN15dxi  = 0.25 * xi0 * (1. - zeta2)*(1. + eta*eta0);
      const CFreal dN15deta = 0.25 * eta0 * (1. - zeta2)*(1. + xi*xi0);
      const CFreal dN15dzta = -0.5 * zeta * (1. + xi*xi0)*(1. + eta*eta0);


      const CFreal JXX = m_invJ(0,0);
      const CFreal JXY = m_invJ(0,1);
      const CFreal JXZ = m_invJ(0,2);
      const CFreal JYX = m_invJ(1,0);
      const CFreal JYY = m_invJ(1,1);
      const CFreal JYZ = m_invJ(1,2);
      const CFreal JZX = m_invJ(2,0);
      const CFreal JZY = m_invJ(2,1);
      const CFreal JZZ = m_invJ(2,2);

      lgrad(0,XX) = dN0dxi * JXX + dN0deta * JXY + dN0dzta * JXZ;
      lgrad(0,YY) = dN0dxi * JYX + dN0deta * JYY + dN0dzta * JYZ;
      lgrad(0,ZZ) = dN0dxi * JZX + dN0deta * JZY + dN0dzta * JZZ;

      lgrad(1,XX) = dN1dxi * JXX + dN1deta * JXY + dN1dzta * JXZ;
      lgrad(1,YY) = dN1dxi * JYX + dN1deta * JYY + dN1dzta * JYZ;
      lgrad(1,ZZ) = dN1dxi * JZX + dN1deta * JZY + dN1dzta * JZZ;

      lgrad(2,XX) = dN2dxi * JXX + dN2deta * JXY + dN2dzta * JXZ;
      lgrad(2,YY) = dN2dxi * JYX + dN2deta * JYY + dN2dzta * JYZ;
      lgrad(2,ZZ) = dN2dxi * JZX + dN2deta * JZY + dN2dzta * JZZ;

      lgrad(3,XX) = dN3dxi * JXX + dN3deta * JXY + dN3dzta * JXZ;
      lgrad(3,YY) = dN3dxi * JYX + dN3deta * JYY + dN3dzta * JYZ;
      lgrad(3,ZZ) = dN3dxi * JZX + dN3deta * JZY + dN3dzta * JZZ;

      lgrad(4,XX) = dN4dxi * JXX + dN4deta * JXY + dN4dzta * JXZ;
      lgrad(4,YY) = dN4dxi * JYX + dN4deta * JYY + dN4dzta * JYZ;
      lgrad(4,ZZ) = dN4dxi * JZX + dN4deta * JZY + dN4dzta * JZZ;

      lgrad(5,XX) = dN5dxi * JXX + dN5deta * JXY + dN5dzta * JXZ;
      lgrad(5,YY) = dN5dxi * JYX + dN5deta * JYY + dN5dzta * JYZ;
      lgrad(5,ZZ) = dN5dxi * JZX + dN5deta * JZY + dN5dzta * JZZ;

      lgrad(6,XX) = dN6dxi * JXX + dN6deta * JXY + dN6dzta * JXZ;
      lgrad(6,YY) = dN6dxi * JYX + dN6deta * JYY + dN6dzta * JYZ;
      lgrad(6,ZZ) = dN6dxi * JZX + dN6deta * JZY + dN6dzta * JZZ;

      lgrad(7,XX) = dN7dxi * JXX + dN7deta * JXY + dN7dzta * JXZ;
      lgrad(7,YY) = dN7dxi * JYX + dN7deta * JYY + dN7dzta * JYZ;
      lgrad(7,ZZ) = dN7dxi * JZX + dN7deta * JZY + dN7dzta * JZZ;

      lgrad(8,XX) = dN8dxi * JXX + dN8deta * JXY + dN8dzta * JXZ;
      lgrad(8,YY) = dN8dxi * JYX + dN8deta * JYY + dN8dzta * JYZ;
      lgrad(8,ZZ) = dN8dxi * JZX + dN8deta * JZY + dN8dzta * JZZ;

      lgrad(9,XX) = dN9dxi * JXX + dN9deta * JXY + dN9dzta * JXZ;
      lgrad(9,YY) = dN9dxi * JYX + dN9deta * JYY + dN9dzta * JYZ;
      lgrad(9,ZZ) = dN9dxi * JZX + dN9deta * JZY + dN9dzta * JZZ;

      lgrad(10,XX) = dN10dxi * JXX + dN10deta * JXY + dN10dzta * JXZ;
      lgrad(10,YY) = dN10dxi * JYX + dN10deta * JYY + dN10dzta * JYZ;
      lgrad(10,ZZ) = dN10dxi * JZX + dN10deta * JZY + dN10dzta * JZZ;

      lgrad(11,XX) = dN11dxi * JXX + dN11deta * JXY + dN11dzta * JXZ;
      lgrad(11,YY) = dN11dxi * JYX + dN11deta * JYY + dN11dzta * JYZ;
      lgrad(11,ZZ) = dN11dxi * JZX + dN11deta * JZY + dN11dzta * JZZ;

      lgrad(12,XX) = dN12dxi * JXX + dN12deta * JXY + dN12dzta * JXZ;
      lgrad(12,YY) = dN12dxi * JYX + dN12deta * JYY + dN12dzta * JYZ;
      lgrad(12,ZZ) = dN12dxi * JZX + dN12deta * JZY + dN12dzta * JZZ;

      lgrad(13,XX) = dN13dxi * JXX + dN13deta * JXY + dN13dzta * JXZ;
      lgrad(13,YY) = dN13dxi * JYX + dN13deta * JYY + dN13dzta * JYZ;
      lgrad(13,ZZ) = dN13dxi * JZX + dN13deta * JZY + dN13dzta * JZZ;

      lgrad(14,XX) = dN14dxi * JXX + dN14deta * JXY + dN14dzta * JXZ;
      lgrad(14,YY) = dN14dxi * JYX + dN14deta * JYY + dN14dzta * JYZ;
      lgrad(14,ZZ) = dN14dxi * JZX + dN14deta * JZY + dN14dzta * JZZ;

      lgrad(15,XX) = dN15dxi * JXX + dN15deta * JXY + dN15dzta * JXZ;
      lgrad(15,YY) = dN15dxi * JYX + dN15deta * JYY + dN15dzta * JYZ;
      lgrad(15,ZZ) = dN15dxi * JZX + dN15deta * JZY + dN15dzta * JZZ;

      lgrad(16,XX) = dN16dxi * JXX + dN16deta * JXY + dN16dzta * JXZ;
      lgrad(16,YY) = dN16dxi * JYX + dN16deta * JYY + dN16dzta * JYZ;
      lgrad(16,ZZ) = dN16dxi * JZX + dN16deta * JZY + dN16dzta * JZZ;

      lgrad(17,XX) = dN17dxi * JXX + dN17deta * JXY + dN17dzta * JXZ;
      lgrad(17,YY) = dN17dxi * JYX + dN17deta * JYY + dN17dzta * JYZ;
      lgrad(17,ZZ) = dN17dxi * JZX + dN17deta * JZY + dN17dzta * JZZ;

      lgrad(18,XX) = dN18dxi * JXX + dN18deta * JXY + dN18dzta * JXZ;
      lgrad(18,YY) = dN18dxi * JYX + dN18deta * JYY + dN18dzta * JYZ;
      lgrad(18,ZZ) = dN18dxi * JZX + dN18deta * JZY + dN18dzta * JZZ;

      lgrad(19,XX) = dN19dxi * JXX + dN19deta * JXY + dN19dzta * JXZ;
      lgrad(19,YY) = dN19dxi * JYX + dN19deta * JYY + dN19dzta * JYZ;
      lgrad(19,ZZ) = dN19dxi * JZX + dN19deta * JZY + dN19dzta * JZZ;

// CFout << "grad(0,XX): " << lgrad(0,XX) <<"\n";
// CFout << "grad(1,XX): " << lgrad(1,XX) <<"\n";
// CFout << "grad(2,XX): " << lgrad(2,XX) <<"\n";
// CFout << "grad(3,XX): " << lgrad(3,XX) <<"\n";
// CFout << "grad(4,XX): " << lgrad(4,XX) <<"\n";
// CFout << "grad(5,XX): " << lgrad(5,XX) <<"\n";
// CFout << "grad(6,XX): " << lgrad(6,XX) <<"\n";
// CFout << "grad(7,XX): " << lgrad(7,XX) <<"\n";
// CFout << "grad(8,XX): " << lgrad(8,XX) <<"\n";
// CFout << "grad(9,XX): " << lgrad(9,XX) <<"\n";
// CFout << "grad(10,XX): " << lgrad(10,XX) <<"\n";
// CFout << "grad(11,XX): " << lgrad(11,XX) <<"\n";
// CFout << "grad(12,XX): " << lgrad(12,XX) <<"\n";
// CFout << "grad(13,XX): " << lgrad(13,XX) <<"\n";
// CFout << "grad(14,XX): " << lgrad(14,XX) <<"\n";
// CFout << "grad(15,XX): " << lgrad(15,XX) <<"\n";
// CFout << "grad(16,XX): " << lgrad(16,XX) <<"\n";
// CFout << "grad(17,XX): " << lgrad(17,XX) <<"\n";
// CFout << "grad(18,XX): " << lgrad(18,XX) <<"\n";
// CFout << "grad(19,XX): " << lgrad(19,XX) <<"\n";

    }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ShapeFunctions



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
