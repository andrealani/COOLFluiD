// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <iomanip>
#include <cmath>
#include <fstream>
#include "Common/PE.hh"
#include "Common/ParallelException.hh"
#include "Common/BadValueException.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/DirPaths.hh"
#include "Framework/MapGeoEnt.hh"
#include "Framework/CFPolyOrder.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/Stopwatch.hh"
#include "Common/CFMap.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Common/Table.hh"
#include "MeshGenerator1D/MeshGenerator1D.hh"
#include "MeshGenerator1D/MeshGenerator1DImpl.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MeshGenerator1D {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MeshGenerator1DImpl,
			    MeshFormatConverter,
			    MeshGenerator1DModule,
			    1>
meshGenerator1DImplProvider("MeshGenerator1D");

//////////////////////////////////////////////////////////////////////////////

MeshGenerator1DImpl::MeshGenerator1DImpl (const std::string& name)
: MeshFormatConverter(name)
{
  addConfigOptionsTo(this);
  
    m_isPeriodicBoundary = false;
    setParameter("PeriodicBoundary",&m_isPeriodicBoundary);

}

////////////////////////////////////////////////////////////////////////////// 

MeshGenerator1DImpl::~MeshGenerator1DImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

void MeshGenerator1DImpl::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("PeriodicBoundary","Are the boundaries periodic?");

}

//////////////////////////////////////////////////////////////////////////////

void MeshGenerator1DImpl::configure ( Config::ConfigArgs& args )
{
  CFLog(INFO, "Configuring MeshGenerator1DImpl \n");
  MeshFormatConverter::configure(args);
}

////////////////////////////////////////////////////////////////////////////// 

void MeshGenerator1DImpl::checkFormat(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void MeshGenerator1DImpl::readFiles(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;

  using namespace boost::filesystem;
  path meshFile = change_extension(filepath, getOriginExtension());

  Common::SelfRegistPtr<Environment::FileHandlerInput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  ifstream& fin = fhandle->open(meshFile);

  std::string ch = "";

  CFreal inlet = 0.;
  CFreal outlet = 0.;
  CFreal throat = 0.;
  CFreal relax = 0.;
  CFreal p1 = 0.;
  CFreal p2 = 0.;
  CFreal factor1 = 0.;
  CFreal factor2 = 0.;
  CFreal lenCellMin = 0.;

  CFuint flagShock = 0;
  CFuint flagUniMesh = 0; 
  CFuint n = 0;
  CFuint m = n;
  CFuint nodes = 0;
  _nbUpdatableNodes = 1;
  
  fin >> ch;
  fin >> ch;
  fin >> inlet;                     // inlet
  fin >> ch;
  fin >> outlet;                    // outlet
  fin >> ch;
  fin >> throat;                    // throat
  fin >> ch;
  fin >> relax;
  fin >> ch;
  fin >> nodes;                     // initial value of the numbers of nodes
  fin >> ch;
  fin >> p1;                        // stretch parameters upstream the throat
  fin >> ch;
  fin >> p2;                        // stretch parameters downstream the throat
  fin >> ch;
  fin >> factor1;                   // power parameter for the stretch
  fin >> ch;
  fin >> factor2;                   // power parameter for the stretch
  fin >> ch;
  fin >> lenCellMin;                // mininum cell size
  fin >> ch;
  fin >> flagShock;                 // flag telling if the mesh must be generated for a steady shock computation 
  fin >> ch;
  fin >> flagUniMesh;               // flag for uniform mesh (suitable for unsteady shock tube computations)

  fin.close();                      // closing file

  CFreal x = 0.0;
  CFreal xold = 0.0;
  CFreal delta = 0.0;
  CFreal lenCell = (outlet - inlet)/(nodes - 1);

  // Check on the minimum cell size before nodes are generated
  if (lenCellMin >= lenCell) { 
     cout << "MeshGenerator1D:: Error, minimum cell size greater than the initial one\n ";
     cout << "MeshGenerator1D:: Please, modify the 'Parameters.dat' file\n";
     abort(); 
  }

  CFreal x1 = throat - (throat - inlet)*p1;
  CFreal x2 = throat + (outlet - throat)*p2;
  
  _coordinate.push_back(0.0);

  if (flagUniMesh == 0) {
    
    // Generating 1D mesh for steady nozzle or shock  
    while (x <= outlet) {
      
      xold = x;
      
      // Converging part: stretch is not applied
      if (x < x1) {
	xold = x;
	++_nbUpdatableNodes;
	x = x + lenCell / std::pow(factor1, (int) n);
	_coordinate.push_back(x);
	delta = x - xold;
      }
      
      // Converging and diverging part: stretch is applies
      if ((x > x1) && (x < x2)) {
	if ( (lenCell / std::pow(factor1, (int) n)) > lenCellMin) {
	  n++;
	  _nbUpdatableNodes++;
	  x = x + lenCell / std::pow(factor1, (int) n);
	  _coordinate.push_back(x);
	  delta = x -xold;
	} else {
	  _nbUpdatableNodes++;
	  x = x + lenCellMin;
	  _coordinate.push_back(x);
	  delta = x - xold;
	}
      }
      
      // Diverging part: stretch is not applied
      if ((x > x2) && (x <= outlet)) {
	if (delta*(std::pow(factor2, (int) m)) < lenCell) {
	  m++;
	  x = x + delta*(std::pow(factor2, (int) m));
	  if (x <= outlet) { 
	    _coordinate.push_back(x);
	    _nbUpdatableNodes++; 
	    delta = x - xold;
	  }
	} else {
	  x = x + lenCell;
	  if (x <= outlet) { 
	    _coordinate.push_back(x);
	    _nbUpdatableNodes++; 
	    delta = x - xold;
	  }
	}
	
      }
    }
    
    // adding more points for steady shock (relaxation zone must be increased) 
    if (flagShock == 1) { 
      
      const CFuint old = _nbUpdatableNodes;
      x = _coordinate[old-1];
      
      while (x <= relax) { 
	delta = lenCell*(std::pow(factor1, (int)(_nbUpdatableNodes - old)));
	
	if (delta <= 0.1) {
	  x = x + delta;
	  if (x <= relax) { 
	    _nbUpdatableNodes++;
	    _coordinate.push_back(x);  
	  } 
	} else {
	  x = x + 0.1;
	  if (x <= relax) { 
	    _nbUpdatableNodes++;
	    _coordinate.push_back(x);  
	  } 
	}
        
      }    
      
    } 
    
    // Uniform mesh for unsteady shock tube computations
  } else {    
    _nbUpdatableNodes = 1;
    x = 0.0;
    while (x <= outlet+1e-8) { 
      x += lenCell;
      if (x <= outlet+1e-8) {
	_coordinate.push_back(x);
	_nbUpdatableNodes++;
      }  
    } 
  } 
  
  // Number of cells
   _nbCells = _nbUpdatableNodes -1;

  // Radius distribution file ('Radius_avg.dat')
  if (flagShock == 0) {
     nozzleRadius();
  } else { 
     shockRadius();
  }

}
////////////////////////////////////////////////////////////////////////////// 10

void MeshGenerator1DImpl::writeContinuousElements(ofstream& fout)
{
  CFAUTOTRACE;
}

////////////////////////////////////////////////////////////////////////////// 11

void MeshGenerator1DImpl::writeContinuousStates(ofstream& fout)
{
  CFAUTOTRACE;
}

///////////////////////////////////////////////////////////////////////////// 11

void MeshGenerator1DImpl::writeNodes(ofstream& fout)
{
  CFAUTOTRACE;

  fout << "!LIST_NODE " << "\n";
  for (CFuint k = 0; k < _nbUpdatableNodes; ++k)
  {
    fout.precision(14);
    fout.setf(ios::scientific,ios::floatfield);
    fout << _coordinate[k] << " ";
    fout << "\n";
  }
}

//////////////////////////////////////////////////////////////////////////////

void MeshGenerator1DImpl::writeDiscontinuousElements(ofstream& fout)
{
  CFAUTOTRACE;

  const CFuint nbNotUpdatableNodes  = 0;
  const CFuint nbNotUpdatableStates = 0;

  fout << "!NB_NODES "
       << _nbUpdatableNodes
       << " "
       << nbNotUpdatableNodes << "\n";

  /// @todo how to define nbNotUpdatableStates for FVM??
  fout << "!NB_STATES "
       << _nbCells
       << " "
       << nbNotUpdatableStates << "\n";

  fout << "!NB_ELEM "        << _nbCells << "\n";
  fout << "!NB_ELEM_TYPES "  << "1"<< "\n";

  //   /// @todo only first order for now

  fout << "!GEOM_POLYORDER " << (CFuint) CFPolyOrder::ORDER1 << "\n";
  fout << "!SOL_POLYORDER "  << (CFuint) CFPolyOrder::ORDER0 << "\n";

  // this CFPolyOrder::ORDER0 can only be set if here we know that we are dealing
  // with CellCenterFEM

  fout << "!ELEM_TYPES ";
  fout << "Line ";
  fout << "\n";

  fout << "!NB_ELEM_PER_TYPE ";
  fout << _nbCells;
  fout << "\n";

  fout << "!NB_NODES_PER_TYPE ";
  fout << "2";
  fout << "\n";

  fout << "!NB_STATES_PER_TYPE ";
  fout << "1";
  fout << "\n";

  fout << "!LIST_ELEM " << "\n";

  for (CFuint countElem = 0; countElem < (_nbUpdatableNodes-1); countElem++)
  {
    fout << countElem << " " << (countElem + 1) << " "<< countElem << "\n";
  }
}

////////////////////////////////////////////////////////////////////////////// 13

void MeshGenerator1DImpl::writeDiscontinuousStates(ofstream& fout)
{
  CFAUTOTRACE;

  fout << "!LIST_STATE " << "0" << "\n";
}

////////////////////////////////////////////////////////////////////////////// 13

void MeshGenerator1DImpl::writeContinuousTrsData(ofstream& fout)
{
  CFAUTOTRACE;
}

////////////////////////////////////////////////////////////////////////////// 14

void MeshGenerator1DImpl::writeDiscontinuousTrsData(ofstream& fout)
{
  CFAUTOTRACE;
  
  if(m_isPeriodicBoundary) {
    fout << "!NB_TRSs " << "1" << "\n";
    fout << "!TRS_NAME " << "Periodic" << "\n";
    fout << "!NB_TRs " <<  "1"<<  "\n";
    fout << "!NB_GEOM_ENTS " << "2" << "\n";
    fout << "!GEOM_TYPE Face" << "\n";
    fout << "!LIST_GEOM_ENT" << "\n";
    fout << "1 " << "1 " << "0 " << "0" << "\n";
    fout << "1 " << "1 " << (_nbUpdatableNodes -1) << " " << (_nbCells -1) << "\n";
  }
  else {
    fout << "!NB_TRSs " << "2" << "\n";
    fout << "!TRS_NAME " << "Inlet" << "\n";
    fout << "!NB_TRs " <<  "1"<<  "\n";
    fout << "!NB_GEOM_ENTS " << "1" << "\n";
    fout << "!GEOM_TYPE Face" << "\n";
    fout << "!LIST_GEOM_ENT" << "\n";
    fout << "1 " << "1 " << "0 " << "0" << "\n";

    fout << "!TRS_NAME " << "Outlet" << "\n";
    fout << "!NB_TRs " <<  "1"<<  "\n";
    fout << "!NB_GEOM_ENTS " << "1" << "\n";
    fout << "!GEOM_TYPE Face" << "\n";
    fout << "!LIST_GEOM_ENT" << "\n";
    fout << "1 " << "1 " << (_nbUpdatableNodes -1) << " " << (_nbCells -1) << "\n";
  }

}

////////////////////////////////////////////////////////////////////////////// 15

void MeshGenerator1DImpl::convertBack(const boost::filesystem::path& filepath)
{
}

////////////////////////////////////////////////////////////////////////////// 16

void MeshGenerator1DImpl::nozzleRadius()
{
  Common::SelfRegistPtr<Environment::FileHandlerInput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  std::string file = "Radius.dat";
  boost::filesystem::path radfilePath(Environment::DirPaths::getInstance().getWorkingDir()/boost::filesystem::path(file));
  ifstream& indata = fhandle->open(radfilePath);
  
  // std::ifstream indata(file.c_str());
  std::vector<CFreal> Radius;
  std::vector<CFreal> derRadius;
  std::vector<CFreal> axial;
  std::vector<CFreal> Rad;
  std::vector<CFreal> RadPrime;

  CFuint nbPoints = 0;
  CFreal ax = 0.0;
  CFreal R = 0.0;
  CFreal Rprime = 0.0;
  CFreal derR = 0.0;
  CFreal slope = 0.0;
  CFreal slope_der = 0.0;
  indata >> nbPoints;
  
  for (CFuint m = 0; m < nbPoints; m++)
  {
    indata >> ax >> R >> derR;
    axial.push_back(ax);
    Radius.push_back(R);
    derRadius.push_back(derR);
  }
  
  indata.close();
  
  CFuint counter = 0;
  for (CFuint i = 0; i < _nbUpdatableNodes; i++)
  {
    if ( _coordinate[i] < axial[counter+1] )
    {
      slope = (Radius[counter+1] - Radius[counter])/(axial[counter+1]- axial[counter]);
      slope_der = (derRadius[counter+1] - derRadius[counter])/(axial[counter+1]- axial[counter]);
      R = slope*(_coordinate[i] - axial[counter]) + Radius[counter];
      Rprime = slope_der*(_coordinate[i] - axial[counter]) + derRadius[counter];
      Rad.push_back(R);
      RadPrime.push_back(Rprime);
    }
    else
    {
	counter++;
	slope = (Radius[counter+1] - Radius[counter])/(axial[counter+1]- axial[counter]);
	slope_der = (derRadius[counter+1] - derRadius[counter])/(axial[counter+1]- axial[counter]);
	R = slope*(_coordinate[i] - axial[counter]) + Radius[counter];
	Rprime = slope_der*(_coordinate[i] - axial[counter]) + derRadius[counter];
	Rad.push_back(R);
	RadPrime.push_back(Rprime);
    }
  }
  
  
  CFreal rAvg = 0.0;
  CFreal derAvg = 0.0;

  Common::SelfRegistPtr<Environment::FileHandlerOutput> fhandleOut =
    Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  
  std::string geom = "Radius_avg.dat";
  boost::filesystem::path radAvgFilePath(Environment::DirPaths::getInstance().getWorkingDir()/boost::filesystem::path(geom));
  ofstream& outdata = fhandleOut->open(radAvgFilePath);
  
  // std::ofstream outdata(geom.c_str());
  outdata << _nbCells << "\n";
  
  for (CFuint k = 0; k < _nbCells; k++)
  {
    rAvg = (Rad[k+1] + Rad[k])/2.0;
    derAvg = (RadPrime[k+1] + RadPrime[k])/2.0;
    outdata << k << " " << -((2*derAvg)/rAvg) << "\n";
  }
  
  outdata.close();
}

////////////////////////////////////////////////////////////////////////////// 17

void MeshGenerator1DImpl::shockRadius() 
{

  Common::SelfRegistPtr<Environment::FileHandlerOutput> fhandleOut =
    Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();

  std::string geom = "Radius_avg.dat";
  boost::filesystem::path radAvgFilePath(Environment::DirPaths::getInstance().getWorkingDir()/boost::filesystem::path(geom));
  ofstream& outdata = fhandleOut->open(radAvgFilePath);
  
  outdata << _nbCells << "\n";
  
  // File used for the geometrical source term
  for (CFuint k = 0; k < _nbCells; k++)
  {
    outdata << k << " " << 0.0 <<"\n";
  }

  outdata.close();

  // Files with node locations (usable for mesh quality checking)
  ofstream fout("Nodes.dat");

  for (CFuint k = 0; k < _nbUpdatableNodes; k++)
  { 
    fout << _coordinate[k] << " "<< 0.0 << "\n"; 
  }
  
  fout.close();
  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshGenerator1D

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
