#ifndef NDEBUG
#include <fstream>
#endif

#include "ComputeInwardNormalsTriagP3.hh"
#include "InwardNormalsData.hh"
#include "Framework/MeshData.hh"
#include "Environment/ObjectProvider.hh"

#include "Framework/SubSystemStatus.hh"
#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/P3Normal.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ComputeInwardNormalsTriagP3,
           ComputeNormals,
               FluctSplitModule>
computeInwardNormalsTriagP3Provider("InwardTriagP3");

//////////////////////////////////////////////////////////////////////////////

ComputeInwardNormalsTriagP3::ComputeInwardNormalsTriagP3() :
ComputeInwardNormals()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeInwardNormalsTriagP3::operator() (const CFuint& iFirstCell,
                                              const CFuint& iLastCell,
                                              const CFuint& iType)
 {
  //set the IDS mapping a given state with the opposite faceID
  // in the local (within the element) connectivity
  vector<CFuint> ids(6);
  ids[0] = 0; // node 0 opposite to face 0
  ids[1] = 4; // node 1 opposite to face 4
  ids[2] = 8; // node 2 opposite to face 8
  ids[3] = 8; // node 3 opposite to face 8
  ids[4] = 0; // node 4 opposite to face 0
  ids[5] = 4; // node 5 opposite to face 4

  InwardNormalsData::setStateToFaceID(iType, ids);

   Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");

   DataHandle<InwardNormalsData*> normals = socket_normals->getDataHandle();

   DataHandle<CFreal> normalsData = socket_normalsData->getDataHandle();

   DataHandle<CFuint> tempSize = socket_tempSize->getDataHandle();

   cf_assert (cells->getNbNodesInGeo(iFirstCell) == 10);
   cf_assert (PhysicalModelStack::getActive()->getDim() == DIM_2D);

   const CFuint nbFacesInCell = 18;
   const CFuint dim = PhysicalModelStack::getActive()->getDim();

   RealMatrix faceNormals(nbFacesInCell, dim);
   RealVector faceAreas(0.0, nbFacesInCell);
   CFreal* nodalNormals = CFNULL;
   CFreal* nodalAreas = CFNULL;

#ifndef NDEBUG
  std::string filename("TriagP3normals.dat");
  ofstream nfile(filename.c_str());
#endif

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
  geoBuilder = getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = cells;

    CFreal xi, eta, invnorm;
    RealVector facenormal;
    facenormal.resize(DIM_2D);
    P3Normal P3N;

    const CFuint nbnodes = 10;
    const CFreal XI[nbnodes] = { 0.0, 1.0, 0.0, 1.0/3.0, 2.0/3.0, 2.0/3.0, 1.0/3.0, 0.0, 0.0, 1.0/3.0 };
    const CFreal ETA[nbnodes] = { 0.0, 0.0, 1.0, 0.0, 0.0, 1.0/3.0, 2.0/3.0, 2.0/3.0, 1.0/3.0, 1.0/3.0 };

    const CFuint SubFaceStart[nbFacesInCell] = { 3, 8, 0, 1, 5, 4, 6, 2, 7, 5, 6, 9, 9, 7, 8, 4, 9, 3 };
    const CFuint SubFaceEnd[nbFacesInCell] =   { 8, 0, 3, 5, 4, 1, 2, 7, 6, 6, 9, 5, 7, 8, 9, 9, 3, 4 };

  for (CFuint iCell = iFirstCell; iCell < iLastCell; ++iCell) {

    // build the GeometricEntity
    geoData.idx = iCell;

    GeometricEntity& currCell = *geoBuilder->buildGE();

    const vector<Node*>* const nodes = currCell.getNodes();
    cf_assert (nodes != CFNULL);
    cf_assert (nodes->size() == cells->getNbNodesInGeo(iFirstCell));

    //Face 0
    xi = eta = 1.0/6.0;
    P3N.ComputeNormal((*nodes), 0, xi, eta, facenormal);
    invnorm = 1.0/facenormal.norm2();
    faceNormals(0,XX) = facenormal[XX]*invnorm;
    faceNormals(0,YY) = facenormal[YY]*invnorm;

    //Face 1
    xi = 0.0; eta = 1.0/6.0;
    P3N.ComputeNormal((*nodes), 1, xi, eta, facenormal);
    invnorm = 1.0/facenormal.norm2();
    faceNormals(1,XX) = facenormal[XX]*invnorm;
    faceNormals(1,YY) = facenormal[YY]*invnorm;

    //Face 2
    xi = 1.0/6.0; eta = 0.0;
    P3N.ComputeNormal((*nodes), 2, xi, eta, facenormal);
    invnorm = 1.0/facenormal.norm2();
    faceNormals(2,XX) = facenormal[XX]*invnorm;
    faceNormals(2,YY) = facenormal[YY]*invnorm;

    //Face 3
    xi = 5.0/6.0; eta = 1.0/6.0;
    P3N.ComputeNormal((*nodes), 3, xi, eta, facenormal);
    invnorm = 1.0/facenormal.norm2();
    faceNormals(3,XX) = facenormal[XX]*invnorm;
    faceNormals(3,YY) = facenormal[YY]*invnorm;

    //Face 4
    xi = 2.0/3.0; eta = 1.0/6.0;
    P3N.ComputeNormal((*nodes), 4, xi, eta, facenormal);
    invnorm = 1.0/facenormal.norm2();
    faceNormals(4,XX) = facenormal[XX]*invnorm;
    faceNormals(4,YY) = facenormal[YY]*invnorm;

    //Face 5
    xi = 5.0/6.0; eta = 0.0;
    P3N.ComputeNormal((*nodes), 5, xi, eta, facenormal);
    invnorm = 1.0/facenormal.norm2();
    faceNormals(5,XX) = facenormal[XX]*invnorm;
    faceNormals(5,YY) = facenormal[YY]*invnorm;

    //Face 6
    xi = 1.0/6.0; eta = 5.0/6.0;
    P3N.ComputeNormal((*nodes), 6, xi, eta, facenormal);
    invnorm = 1.0/facenormal.norm2();
    faceNormals(6,XX) = facenormal[XX]*invnorm;
    faceNormals(6,YY) = facenormal[YY]*invnorm;

    //Face 7
    xi = 0.0; eta = 5.0/6.0;
    P3N.ComputeNormal((*nodes), 7, xi, eta, facenormal);
    invnorm = 1.0/facenormal.norm2();
    faceNormals(7,XX) = facenormal[XX]*invnorm;
    faceNormals(7,YY) = facenormal[YY]*invnorm;

    //Face 8
    xi = 1.0/6.0; eta = 2.0/3.0;
    P3N.ComputeNormal((*nodes), 8, xi, eta, facenormal);
    invnorm = 1.0/facenormal.norm2();
    faceNormals(8,XX) = facenormal[XX]*invnorm;
    faceNormals(8,YY) = facenormal[YY]*invnorm;

    //Face 9
    xi = eta = 0.5;
    P3N.ComputeNormal((*nodes), 9, xi, eta, facenormal);
    invnorm = 1.0/facenormal.norm2();
    faceNormals(9,XX) = facenormal[XX]*invnorm;
    faceNormals(9,YY) = facenormal[YY]*invnorm;

    //Face 10
    xi = 1.0/3.0; eta = 0.5;
    P3N.ComputeNormal((*nodes), 10, xi, eta, facenormal);
    invnorm = 1.0/facenormal.norm2();
    faceNormals(10,XX) = facenormal[XX]*invnorm;
    faceNormals(10,YY) = facenormal[YY]*invnorm;

    //Face 11
    xi = 0.5; eta = 1.0/3.0;
    P3N.ComputeNormal((*nodes), 11, xi, eta, facenormal);
    invnorm = 1.0/facenormal.norm2();
    faceNormals(11,XX) = facenormal[XX]*invnorm;
    faceNormals(11,YY) = facenormal[YY]*invnorm;

    //Face 12
    xi = 1.0/6.0; eta = 0.5;
    P3N.ComputeNormal((*nodes), 12, xi, eta, facenormal);
    invnorm = 1.0/facenormal.norm2();
    faceNormals(12,XX) = facenormal[XX]*invnorm;
    faceNormals(12,YY) = facenormal[YY]*invnorm;

    //Face 13
    xi = 0.0; eta = 0.5;
    P3N.ComputeNormal((*nodes), 13, xi, eta, facenormal);
    invnorm = 1.0/facenormal.norm2();
    faceNormals(13,XX) = facenormal[XX]*invnorm;
    faceNormals(13,YY) = facenormal[YY]*invnorm;

    //Face 14
    xi = 1.0/6.0; eta = 1.0/3.0;
    P3N.ComputeNormal((*nodes), 14, xi, eta, facenormal);
    invnorm = 1.0/facenormal.norm2();
    faceNormals(14,XX) = facenormal[XX]*invnorm;
    faceNormals(14,YY) = facenormal[YY]*invnorm;

    //Face 15
    xi = 0.5; eta = 1.0/6.0;
    P3N.ComputeNormal((*nodes), 15, xi, eta, facenormal);
    invnorm = 1.0/facenormal.norm2();
    faceNormals(15,XX) = facenormal[XX]*invnorm;
    faceNormals(15,YY) = facenormal[YY]*invnorm;

    //Face 16
    xi = 1.0/3.0; eta = 1.0/6.0;
    P3N.ComputeNormal((*nodes), 16, xi, eta, facenormal);
    invnorm = 1.0/facenormal.norm2();
    faceNormals(16,XX) = facenormal[XX]*invnorm;
    faceNormals(16,YY) = facenormal[YY]*invnorm;

    //Face 17
    xi = 0.5; eta = 0.0;
    P3N.ComputeNormal((*nodes), 17, xi, eta, facenormal);
    invnorm = 1.0/facenormal.norm2();
    faceNormals(17,XX) = facenormal[XX]*invnorm;
    faceNormals(17,YY) = facenormal[YY]*invnorm;


	//*********************************************************************************************************
	//
	//				Compute lengths of edges numerically:
	//
	//*********************************************************************************************************

  const CFreal s  = std::sqrt( 0.6 );
  const CFreal a0 = ( 1.0 - s )*0.5;
  const CFreal a1 = ( 1.0 + s )*0.5;

  CFreal qd0[3], qd1[3], wqd[3];
 
  qd0[0] = a0;  qd1[0] = a1;
  qd0[1] = a1;  qd1[1] = a0;
  qd0[2] = .5;  qd1[2] = .5;

  wqd[0] = 5.0/54.0;
  wqd[1] = 5.0/54.0;
  wqd[2] = 8.0/54.0;

  for(CFuint iface = 0;iface < nbFacesInCell; ++iface) {
      faceAreas[iface] = 0.0;
     
    for(CFuint iqdpt = 0; iqdpt<3; ++iqdpt) {
      xi  = qd0[iqdpt]*XI[SubFaceStart[iface]] + qd1[iqdpt]*XI[SubFaceEnd[iface]];
      eta = qd0[iqdpt]*ETA[SubFaceStart[iface]] + qd1[iqdpt]*ETA[SubFaceEnd[iface]];
      P3N.ComputeNormal((*nodes), iface, xi, eta, facenormal);
      faceAreas[iface] += wqd[iqdpt] * facenormal.norm2();
    }

  }

  ///Scale normals by corresponding face lengths	
	for(CFuint iface=0;iface<nbFacesInCell;++iface)
	{
	faceNormals(iface,XX) *= faceAreas[iface];
	faceNormals(iface,YY) *= faceAreas[iface];	
	}


#ifndef NDEBUG
    nfile << "Cell: "     << "\n" << iCell        << "\n";
    nfile << "Nodes: "    << "\n";
    for(CFuint iNode = 0; iNode < cells->getNbNodesInGeo(iFirstCell); ++iNode){
      nfile << (*(*nodes)[iNode]) << "\n";
    }
    nfile << "NNormals: " << "\n" << nodalNormals ;
    nfile << "FNormals: " << "\n" << faceNormals  ;
    nfile << "NAreas: "   << "\n" << nodalAreas   << "\n";
    nfile << "FAreas: "   << "\n" << faceAreas    << "\n" << "\n";
#endif

    CFreal *const startFaceData = &normalsData[tempSize[0]];
    for (CFuint i = 0; i < faceNormals.size(); ++i) {
      normalsData[tempSize[0]++] = faceNormals[i];
    }

    CFreal *const startFaceArea = &normalsData[tempSize[0]];
    for (CFuint i = 0; i < faceAreas.size(); ++i) {
      normalsData[tempSize[0]++] = faceAreas[i];
    }

    normals[iCell] = new InwardNormalsData(startFaceData,
                                           startFaceArea,
                                           nodalNormals,
                                           nodalAreas,
                                           nbFacesInCell,
                                           iType);
    // release the GeometricEntity
    geoBuilder->releaseGE();
   }
 }

//////////////////////////////////////////////////////////////////////////////

void ComputeInwardNormalsTriagP3::update(const CFuint& iFirstCell,
                                         const CFuint& iLastCell,
                                         const CFuint& iType)
{
  throw Common::NotImplementedException (FromHere(),"ComputeInwardNormalsTriagP3::update()");
}

//////////////////////////////////////////////////////////////////////////////

void ComputeInwardNormalsTriagP3::average(const CFuint& iFirstCell,
                                          const CFuint& iLastCell,
                                          const CFuint& iType)
{
  throw Common::NotImplementedException (FromHere(),"ComputeInwardNormalsTriagP3::average()");
}

//////////////////////////////////////////////////////////////////////////////


    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
