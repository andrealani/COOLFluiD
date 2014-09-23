// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "VolumeCalculator.hh"
#include "Common/CFLog.hh"
#include "Framework/Node.hh"
#include "MathTools/MathChecks.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

VolumeCalculator::VolumeCalculator() :
  _detMat(4,4),
  _center(0.0, 3)
{
  _detMat(0,3) = 1.0;
  _detMat(1,3) = 1.0;
  _detMat(2,3) = 1.0;
  _detMat(3,3) = 1.0;
}

//////////////////////////////////////////////////////////////////////////////

VolumeCalculator::~VolumeCalculator()
{
}

//////////////////////////////////////////////////////////////////////////////

CFreal VolumeCalculator::calculateTriagVolume(const vector<Node*>& coord)
{
  RealMatrix matrix(3,3);
  for (CFuint i = 0; i < 3; ++i) {
    for (CFuint j = 0; j < 3; ++j) {
      if (j > 0) {
        matrix(i,j) = (*coord[i])[j-1];
      }
      else {
        matrix(i,j) = 1.0;
      }
    }
  }

  const CFreal volume = 0.5*matrix.determ3();

  return volume;
}

//////////////////////////////////////////////////////////////////////////////

bool VolumeCalculator::checkTriagNumbering(const vector<Node*>& coord)
{

  CFreal volume = calculateTriagVolume(coord);

  bool positive(true);
  if(volume < 0.) positive = false;

  return positive;
}

//////////////////////////////////////////////////////////////////////////////

CFreal VolumeCalculator::calculate3DTriagVolume(const vector<Node*>& coord)
{
  static RealVector vec1(3);
  static RealVector vec2(3);
  static RealVector cross(3);

  for (CFuint iDim = 0; iDim < 3; ++iDim) {
    vec1[iDim] = (*coord[1])[iDim] - (*coord[0])[iDim];
    vec2[iDim] = (*coord[2])[iDim] - (*coord[0])[iDim];
  }

  MathTools::MathFunctions::crossProd(vec1,vec2,cross);

  CFreal volume = sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);

  volume *= 0.5;

  return volume;
}

//////////////////////////////////////////////////////////////////////////////

CFreal VolumeCalculator::calculateQuadVolume(const std::vector<Node*>& nodes)
{
  const CFreal diagonalsProd = ((*nodes[2])[0] - (*nodes[0])[0])*
    ((*nodes[3])[1] - (*nodes[1])[1])-
    ((*nodes[2])[1] - (*nodes[0])[1])*
    ((*nodes[3])[0] - (*nodes[1])[0]);

  const CFreal volume = 0.5*diagonalsProd;

  return volume;
}

//////////////////////////////////////////////////////////////////////////////

bool VolumeCalculator::checkQuadNumbering(const vector<Node*>& coord)
{

  CFreal volume = calculateQuadVolume(coord);

  bool positive = true;
  if(volume < 0.) positive = false;

  return positive;
}

//////////////////////////////////////////////////////////////////////////////

CFreal VolumeCalculator::calculateTetraVolume(const vector<Node*>& coord)
{
  static RealMatrix matrix(4,4);
  for ( CFuint i = 0; i < 4; ++i)  {
    for ( CFuint j = 0; j < 4; ++j) {
      if ( j > 0 )  {
        matrix(i,j) = (*coord[i])[j-1];
      }
      else {
        matrix(i,j) = 1.0;
      }
    }
  }

  const CFreal volume = (1./6.) * matrix.determ4();

  return volume;
}
//////////////////////////////////////////////////////////////////////////////


CFreal VolumeCalculator::calculateTetraVolume(const RealMatrix& coord)
{
  cf_assert(coord.nbRows() == 4);
  cf_assert(coord.nbCols() == 3);

  CFreal volume = 0.;

  setTetraIDs(coord,0,1,2,3);
  volume = calcVolume();

  return volume;
}


//////////////////////////////////////////////////////////////////////////////

bool VolumeCalculator::checkTetraNumbering(const RealMatrix& coord)
{
  CFreal volume = calculateTetraVolume(coord);

  bool positive = true;
  if(volume < 0.) positive = false;

  return positive;
}

//////////////////////////////////////////////////////////////////////////////

CFreal VolumeCalculator::calculatePyramVolume(const RealMatrix& coord)
{
  cf_assert(coord.nbRows() == 5);
  cf_assert(coord.nbCols() == 3);
  
  CFreal volume = 0.;
  CFreal partial1 = 0.;
  CFreal partial2 = 0.;
  
  setTetraIDs(coord,0,1,3,4);
  // setTetraIDs(coord,0,1,2,4);    
  partial1 = calcVolume();
  
  setTetraIDs(coord,1,2,3,4);
  //setTetraIDs(coord,0,2,3,4);
  partial2 = calcVolume();
  
  // cf_assert(partial1*partial2 > 0.);
  if (partial1*partial2 < 0.) {
    CFLog(WARN, "WARNING: VolumeCalculator::calculatePyramVolume() => partial1*partial2 < 0.\n");
  }
  volume = partial1 + partial2;
  
  return volume;
}

//////////////////////////////////////////////////////////////////////////////

CFreal VolumeCalculator::calculatePyramVolume(const std::vector<Node*>& coord)
{
  cf_assert(coord.size() == 5);
  cf_assert(coord[0]->size() == 3);

  CFreal volume = 0.;
  CFreal partial = 0.;

  setTetraIDs(coord,0,1,3,4);
  partial = calcVolume();
  volume += partial;

  setTetraIDs(coord,1,2,3,4);
  partial = calcVolume();
  volume += partial;

  return volume;
}

//////////////////////////////////////////////////////////////////////////////

bool VolumeCalculator::checkPyramNumbering(const RealMatrix& coord)
{
  cf_assert(coord.nbRows() == 5);
  cf_assert(coord.nbCols() == 3);

  CFreal volume = 0.;
  CFreal partial;
  bool positive = true;

  setTetraIDs(coord,0,1,3,4);
  partial = calcVolume();
  positive &= MathChecks::isPositive(partial);
  volume += partial;

  setTetraIDs(coord,1,2,3,4);
  partial = calcVolume();
  positive &= MathChecks::isPositive(partial);
  volume += partial;

  return positive;
}

//////////////////////////////////////////////////////////////////////////////

CFreal VolumeCalculator::calculatePrismVolume(const std::vector<Node*>& coord)
{
  cf_assert(coord.size() == 6);
  cf_assert(coord[0]->size() == 3);
  
  CFreal volume = 0.;
  CFreal partial;

  setTetraIDs(coord,0,4,5,3);
  partial = calcVolume();
  volume += std::abs(partial);

  setTetraIDs(coord,0,1,2,5);
  partial = calcVolume();
  volume += std::abs(partial);

  setTetraIDs(coord,0,4,1,5);
  partial = calcVolume();
  volume += std::abs(partial);

  return volume;
}

//////////////////////////////////////////////////////////////////////////////

CFreal VolumeCalculator::calculatePrismVolume(const RealMatrix& coord)
{
  cf_assert(coord.nbRows() == 6);
  cf_assert(coord.nbCols() == 3);
  
  CFreal volume = 0.;
  CFreal partial = 0;

  setTetraIDs(coord,0,4,5,3);
  partial = calcVolume();
  volume += std::abs(partial);
  
  setTetraIDs(coord,0,1,2,5);
  partial = calcVolume();
  volume += std::abs(partial);

  setTetraIDs(coord,0,4,1,5);
  partial = calcVolume();
  volume += std::abs(partial);

  return volume;
}

//////////////////////////////////////////////////////////////////////////////

bool VolumeCalculator::checkPrismNumbering(const RealMatrix& coord)
{
  cf_assert(coord.nbRows() == 6);
  cf_assert(coord.nbCols() == 3);

  CFreal volume = 0.;
  CFreal partial = 0.;
  bool positive = true;

  setTetraIDs(coord,0,4,5,3);
  partial = calcVolume();
  positive &= MathChecks::isPositive(partial);
  volume += partial;

  setTetraIDs(coord,0,1,2,5);
  partial = calcVolume();
  positive &= MathChecks::isPositive(partial);
  volume += partial;

  setTetraIDs(coord,0,4,1,5);
  partial = calcVolume();
  positive &= MathChecks::isPositive(partial);
  volume += partial;

  return positive;
}

//////////////////////////////////////////////////////////////////////////////

CFreal VolumeCalculator::calculateHexaVolume(const std::vector<Node*>& coord)
{ 
  cf_assert(coord.size() == 8);
  cf_assert(coord[0]->size() == 3);
  
  // cout.precision(12); cout << "6Tet = " << calculateHexaVolumeBy6Tetra(coord) << endl;
  // cout.precision(12); cout << "6Pyr = " << calculateHexaVolumeBy6Pyram(coord) << endl;
  // cout.precision(12); cout << "LD   = " << calculateHexaVolumeByLD(coord) << endl;
  // cout.precision(12); cout << "TH   = " << calculateHexaVolumeByTH(coord) << endl <<endl;
  
  return calculateHexaVolumeByTH(coord);
}

//////////////////////////////////////////////////////////////////////////////

CFreal VolumeCalculator::calculateHexaVolumeBy6Tetra(const std::vector<Node*>& coord)
{
  CFreal volume = 0.;
  CFreal partial = 0.;
    
  setTetraIDs(coord,0,5,7,4);
  partial = calcVolume(); 
  if (partial < 0.0) {CFout << "volume(0,5,7,4) < 0 => " <<  partial  << "\n";}
  volume += std::abs(partial);
  
  setTetraIDs(coord,0,1,3,7);
  partial = calcVolume();
  if (partial < 0.0) {CFout << "volume(0,1,3,7) < 0 => " <<  partial  << "\n";}
  volume += std::abs(partial);
  
  setTetraIDs(coord,0,5,1,7);
  partial = calcVolume();
  if (partial < 0.0) {CFout << "volume(0,5,1,7) < 0 => " <<  partial  << "\n";}
  volume += std::abs(partial);
  
  setTetraIDs(coord,1,6,7,5);
  partial = calcVolume();
  if (partial < 0.0) {CFout << "volume(1,6,7,5) < 0 => " <<  partial  << "\n";}
  volume += std::abs(partial);
  
  setTetraIDs(coord,1,2,3,7);
  partial = calcVolume();
  if (partial < 0.0) {CFout << "volume(1,2,3,7) < 0 => " <<  partial  << "\n";}
  volume += std::abs(partial);
  
  setTetraIDs(coord,1,6,2,7);
  partial = calcVolume();
  if (partial < 0.0) {CFout << "volume(1,6,2,7) < 0 => " <<  partial  << "\n";}
  volume += std::abs(partial);
  
  return volume;
}

//////////////////////////////////////////////////////////////////////////////

CFreal VolumeCalculator::calculateHexaVolumeBy6Pyram(const std::vector<Node*>& coord)
{
  static RealMatrix mat(5,3, 0.0);

  _center = 0.0;
  for (CFuint i = 0; i < coord.size(); ++i) {
    _center += *coord[i];
  }
  _center /= 8;

  CFreal volume = 0.;
  CFreal partial = 0.;
  
  mat.setRow(*coord[0],0);
  mat.setRow(*coord[3],1);
  mat.setRow(*coord[2],2);
  mat.setRow(*coord[1],3);
  mat.setRow(_center,4);
  partial = calculatePyramVolume(mat); 
  volume += std::abs(partial);

//   if (-partial < 0.0) {
//     CFout << "0321 < 0.0" << -partial << "\n";
//   }

  mat.setRow(*coord[4],0);
  mat.setRow(*coord[5],1);
  mat.setRow(*coord[6],2);
  mat.setRow(*coord[7],3);
  mat.setRow(_center,4);
  partial = calculatePyramVolume(mat);
  volume += std::abs(partial);

  // if (-partial < 0.0) {
//     CFout << "4567 < 0.0"<< -partial << "\n";
//   }

  mat.setRow(*coord[0],0);
  mat.setRow(*coord[1],1);
  mat.setRow(*coord[5],2);
  mat.setRow(*coord[4],3);
  mat.setRow(_center,4);
  partial = calculatePyramVolume(mat);
  volume += std::abs(partial);

 //  if (-partial < 0.0) {
//     CFout << "0154 < 0.0"<<-partial << "\n";
//   }

  mat.setRow(*coord[1],0);
  mat.setRow(*coord[2],1);
  mat.setRow(*coord[6],2);
  mat.setRow(*coord[5],3);
  mat.setRow(_center,4);
  partial = calculatePyramVolume(mat);
  volume += std::abs(partial);
  
 //  if (-partial < 0.0) {
//     CFout << "1265 < 0.0"<<-partial << "\n";
//   }

  mat.setRow(*coord[3],0);
  mat.setRow(*coord[7],1);
  mat.setRow(*coord[6],2);
  mat.setRow(*coord[2],3);
  mat.setRow(_center,4);
  partial = calculatePyramVolume(mat);
  volume += std::abs(partial);

//   if (-partial < 0.0) {
//     CFout << "3762 < 0.0"<< -partial << "\n";
//   }

  mat.setRow(*coord[0],0);
  mat.setRow(*coord[4],1);
  mat.setRow(*coord[7],2);
  mat.setRow(*coord[3],3);
  mat.setRow(_center,4);
  partial = calculatePyramVolume(mat);
  volume += std::abs(partial);

  // if (-partial < 0.0) {
//     CFout << "0473 < 0.0"<<-partial << "\n";
//   }
  
  return volume;
}

//////////////////////////////////////////////////////////////////////////////

CFreal VolumeCalculator::calculateHexaVolumeByLD(const std::vector<Node*>& coord)
{
  static RealVector subX6X0(3);
  static RealVector sumX1X7(3);
  static RealVector sumX4X5(3);
  static RealVector subX1X7(3);
  static RealVector subX4X5(3);
  static RealVector sumX2X3(3);
  static RealVector subX2X3(3);
  static RealVector tmp1(3);
  static RealVector tmp2(3);
  static RealMatrix det(3,3, 0.0);
  
  subX6X0 = (*coord[6]) - (*coord[0]);
  sumX1X7 = (*coord[1]) + (*coord[7]);
  sumX4X5 = (*coord[4]) + (*coord[5]);
  tmp1 = sumX1X7 - sumX4X5;
  subX1X7 = (*coord[1]) - (*coord[7]);
  subX4X5 = (*coord[4]) - (*coord[5]);
  tmp2 = subX4X5 - subX1X7;
  
  // set the determinant
  det.setRow(subX6X0,0);
  det.setRow(tmp1,1);
  det.setRow(tmp2,2);
  
  CFreal volume = det.determ3();
  
  sumX2X3 = (*coord[2]) + (*coord[3]);
  tmp1 = sumX1X7 - sumX2X3;
  subX2X3 = (*coord[2]) - (*coord[3]);
  tmp2 = subX1X7 + subX2X3;
  
  // set the determinant
  det.setRow(tmp1,1);
  det.setRow(tmp2,2);
    
  volume += det.determ3(); 
  
  if (volume < 0.0) {
    CFout << "VolumeCalculator::calculateHexaVolumeByLD() => volume < 0 => " <<  volume/12. << "\n";
  }
  
  return volume/12.;
}

//////////////////////////////////////////////////////////////////////////////

CFreal VolumeCalculator::calculateHexaVolumeByTH(const std::vector<Node*>& coord)
{
  static RealVector subX6X1(3);
  static RealVector subX7X0(3);
  static RealVector subX6X3(3);
  static RealVector subX2X0(3);
  static RealVector subX5X0(3);
  static RealVector subX6X4(3);
  static RealVector tmp(3);
  static RealMatrix det(3,3, 0.0);
  
  subX6X1 = (*coord[6]) - (*coord[1]);
  subX7X0 = (*coord[7]) - (*coord[0]);
  subX6X3 = (*coord[6]) - (*coord[3]);
  subX2X0 = (*coord[2]) - (*coord[0]);
  subX5X0 = (*coord[5]) - (*coord[0]);
  subX6X4 = (*coord[6]) - (*coord[4]);
  
  tmp = subX6X1 + subX7X0;
  det.setRow(tmp, 0);
  det.setRow(subX6X3,1);
  det.setRow(subX2X0,2);
  CFreal volume = det.determ3();
  
  tmp = subX6X3 + subX5X0;
  det.setRow(subX7X0, 0);
  det.setRow(tmp,1);
  det.setRow(subX6X4,2);
  volume += det.determ3();
  
  tmp = subX6X4 + subX2X0;
  det.setRow(subX6X1, 0);
  det.setRow(subX5X0,1);
  det.setRow(tmp,2);
  volume += det.determ3();
  
  if (volume < 0.0) {
    CFout << "VolumeCalculator::calculateHexaVolumeByTH() => volume < 0 => " <<  volume/12. << "\n";
  }
  
  return volume/12.;
}

//////////////////////////////////////////////////////////////////////////////

CFreal VolumeCalculator::calculateHexaVolume(const RealMatrix& coord)
{
  cf_assert(coord.nbRows() == 8);
  cf_assert(coord.nbCols() == 3);

  CFreal volume = 0.;
  CFreal partial = 0.;
  
  setTetraIDs(coord,0,5,7,4);
  partial = calcVolume(); 
  if (partial < 0.0) {CFout << "volume(0,5,7,4) < 0 => " <<  partial  << "\n";}
  volume += partial;
  
  setTetraIDs(coord,0,1,3,7);
  partial = calcVolume();
  if (partial < 0.0) {CFout << "volume(0,1,3,7) < 0 => " <<  partial  << "\n";}
  volume += partial;
  
  setTetraIDs(coord,0,5,1,7);
  partial = calcVolume();
  if (partial < 0.0) {CFout << "volume(0,5,1,7) < 0 => " <<  partial  << "\n";}
  volume += partial;

  setTetraIDs(coord,1,6,7,5);
  partial = calcVolume();
  if (partial < 0.0) {CFout << "volume(1,6,7,5) < 0 => " <<  partial  << "\n";}
  volume += partial;

  setTetraIDs(coord,1,2,3,7);
  partial = calcVolume();
  if (partial < 0.0) {CFout << "volume(1,2,3,7) < 0 => " <<  partial  << "\n";}
  volume += partial;
  
  setTetraIDs(coord,1,6,2,7);
  partial = calcVolume();
  if (partial < 0.0) {CFout << "volume(1,6,2,7) < 0 => " <<  partial  << "\n";}
  volume += partial;
  
  return volume;
}

//////////////////////////////////////////////////////////////////////////////

bool VolumeCalculator::checkHexaNumbering(const RealMatrix& coord)
{
  cf_assert(coord.nbRows() == 8);
  cf_assert(coord.nbCols() == 3);

  CFreal volume = 0.;
  CFreal partial = 0.;
  bool positive = true;

  setTetraIDs(coord,0,5,7,4);
  partial = calcVolume();
  positive = MathChecks::isPositive(partial);
  volume += partial;

  setTetraIDs(coord,0,1,3,7);
  partial = calcVolume();
  positive = MathChecks::isPositive(partial);
  volume += partial;

  setTetraIDs(coord,0,5,1,7);
  partial = calcVolume();
  positive = MathChecks::isPositive(partial);
  volume += partial;

  setTetraIDs(coord,1,6,7,5);
  partial = calcVolume();
  positive = MathChecks::isPositive(partial);
  volume += partial;

  setTetraIDs(coord,1,2,3,7);
  partial = calcVolume();
  positive = MathChecks::isPositive(partial);
  volume += partial;

  setTetraIDs(coord,1,6,2,7);
  partial = calcVolume();
  positive = MathChecks::isPositive(partial);
  volume += partial;

  return positive;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

