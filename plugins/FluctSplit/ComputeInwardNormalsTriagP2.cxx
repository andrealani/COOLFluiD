#ifndef NDEBUG
#include <fstream>
#endif

#include "ComputeInwardNormalsTriagP2.hh"
#include "InwardNormalsData.hh"
#include "Framework/MeshData.hh"
#include "Environment/ObjectProvider.hh"

#include "Framework/SubSystemStatus.hh"
#include "FluctSplit/FluctSplit.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ComputeInwardNormalsTriagP2,
           ComputeNormals,
               FluctSplitModule>
computeInwardNormalsTriagP2Provider("InwardTriagP2");

//////////////////////////////////////////////////////////////////////////////

ComputeInwardNormalsTriagP2::ComputeInwardNormalsTriagP2() :
ComputeInwardNormals()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeInwardNormalsTriagP2::operator() (const CFuint& iFirstCell,
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

   cf_assert (cells->getNbNodesInGeo(iFirstCell) == 6);
   cf_assert (PhysicalModelStack::getActive()->getDim() == DIM_2D);

   const CFuint nbFacesInCell = 9;
   const CFuint dim = PhysicalModelStack::getActive()->getDim();

   RealMatrix faceNormals(nbFacesInCell, dim);
   RealVector faceAreas(0.0, nbFacesInCell);
   CFreal* nodalNormals = CFNULL;
   CFreal* nodalAreas = CFNULL;

#ifndef NDEBUG
  std::string filename("TriagP2normals.dat");
  ofstream nfile(filename.c_str());
#endif

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
  geoBuilder = getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = cells;

  for (CFuint iCell = iFirstCell; iCell < iLastCell; ++iCell) {

    // build the GeometricEntity
    geoData.idx = iCell;

    GeometricEntity& currCell = *geoBuilder->buildGE();

    const vector<Node*>* const nodes = currCell.getNodes();
    cf_assert (nodes != CFNULL);
    cf_assert (nodes->size() == cells->getNbNodesInGeo(iFirstCell));

    CFreal xi, eta, norm, nx, ny;

	//Face 0
	xi = 0.25;
	//eta = 0.25;
	nx = -(4*xi-1)*(*(*nodes)[1])[1]-(4*xi-1)*(*(*nodes)[2])[1]-2*(*(*nodes)[3])[1]-\
			(2-8*xi)*(*(*nodes)[4])[1]+2*(*(*nodes)[5])[1];
	ny =  (4*xi-1)*(*(*nodes)[1])[0]+(4*xi-1)*(*(*nodes)[2])[0]+2*(*(*nodes)[3])[0]+\
			(2-8*xi)*(*(*nodes)[4])[0]-2*(*(*nodes)[5])[0];
	norm = std::sqrt(nx*nx+ny*ny);
	faceNormals(0,0) = nx/norm;
	faceNormals(0,1) = ny/norm;
	
	//Face 1
	//xi = 0.0;
	eta = 0.25;
	nx = -(-3+4*eta)*(*(*nodes)[0])[1]-(4*eta-1)*(*(*nodes)[2])[1]-(4-8*eta)*(*(*nodes)[5])[1];
	ny =  (-3+4*eta)*(*(*nodes)[0])[0]+(4*eta-1)*(*(*nodes)[2])[0]+(4-8*eta)*(*(*nodes)[5])[0];
	norm = std::sqrt(nx*nx+ny*ny);
	faceNormals(1,0) = nx/norm;
	faceNormals(1,1) = ny/norm;

	//Face 2
	xi = 0.25;
	//eta = 0.0;
	nx =  (-3+4*xi)*(*(*nodes)[0])[1]+(4*xi-1)*(*(*nodes)[1])[1]+(4-8*xi)*(*(*nodes)[3])[1];
	ny = -(-3+4*xi)*(*(*nodes)[0])[0]-(4*xi-1)*(*(*nodes)[1])[0]-(4-8*xi)*(*(*nodes)[3])[0];
	norm = std::sqrt(nx*nx+ny*ny);
	faceNormals(2,0) = nx/norm;
	faceNormals(2,1) = ny/norm;

	//Face 3
	xi = 0.75;
	//eta = 0.0;
	nx =  (-3+4*xi)*(*(*nodes)[0])[1]+(4*xi-1)*(*(*nodes)[1])[1]+(4-8*xi)*(*(*nodes)[3])[1];
	ny = -(-3+4*xi)*(*(*nodes)[0])[0]-(4*xi-1)*(*(*nodes)[1])[0]-(4-8*xi)*(*(*nodes)[3])[0];
	norm = std::sqrt(nx*nx+ny*ny);
	faceNormals(5,0) = nx/norm;
	faceNormals(5,1) = ny/norm;

	//Face 4
	xi = 0.75;
	//eta = 0.25;
	nx = -(4*xi-1)*(*(*nodes)[1])[1]-(-3+4*xi)*(*(*nodes)[2])[1]-(4-8*xi)*(*(*nodes)[4])[1];
	ny =  (4*xi-1)*(*(*nodes)[1])[0]+(-3+4*xi)*(*(*nodes)[2])[0]+(4-8*xi)*(*(*nodes)[4])[0];
	norm = std::sqrt(nx*nx+ny*ny);
	faceNormals(3,0) = nx/norm;
	faceNormals(3,1) = ny/norm;

	//Face 5
	//xi = 0.5;
	eta = 0.25;
	nx = -(4*eta-1)*(*(*nodes)[0])[1]-(4*eta-1)*(*(*nodes)[2])[1]+2*(*(*nodes)[3])[1]-\
						2*(*(*nodes)[4])[1]-(2-8*eta)*(*(*nodes)[5])[1];
	ny =  (4*eta-1)*(*(*nodes)[0])[0]+(4*eta-1)*(*(*nodes)[2])[0]-2*(*(*nodes)[3])[0]+\
						2*(*(*nodes)[4])[0]+(2-8*eta)*(*(*nodes)[5])[0];
	norm = std::sqrt(nx*nx+ny*ny);
	faceNormals(4,0) = nx/norm;
	faceNormals(4,1) = ny/norm;
	
	//Face 6
	xi = 0.25;
	//eta = 0.5;
	nx =  (4*xi-1)*(*(*nodes)[0])[1]+(4*xi-1)*(*(*nodes)[1])[1]+(2-8*xi)*(*(*nodes)[3])[1]+\
				2*(*(*nodes)[4])[1]-2*(*(*nodes)[5])[1];
	ny = -(4*xi-1)*(*(*nodes)[0])[0]-(4*xi-1)*(*(*nodes)[1])[0]-(2-8*xi)*(*(*nodes)[3])[0]-\
				2*(*(*nodes)[4])[0]+2*(*(*nodes)[5])[0];
	norm = std::sqrt(nx*nx+ny*ny);
	faceNormals(8,0) = nx/norm;
	faceNormals(8,1) = ny/norm;

	//Face 7
	xi = 0.25;
	//eta = 0.75;
	nx = -(4*xi-1)*(*(*nodes)[1])[1]-(-3+4*xi)*(*(*nodes)[2])[1]-(4-8*xi)*(*(*nodes)[4])[1];
	ny =  (4*xi-1)*(*(*nodes)[1])[0]+(-3+4*xi)*(*(*nodes)[2])[0]+(4-8*xi)*(*(*nodes)[4])[0];
	norm = std::sqrt(nx*nx+ny*ny);
	faceNormals(6,0) = nx/norm;
	faceNormals(6,1) = ny/norm;

	//Face 8
	//xi = 0.0;
	eta = 0.75;
	nx = -(-3+4*eta)*(*(*nodes)[0])[1]-(4*eta-1)*(*(*nodes)[2])[1]-(4-8*eta)*(*(*nodes)[5])[1];
	ny =  (-3+4*eta)*(*(*nodes)[0])[0]+(4*eta-1)*(*(*nodes)[2])[0]+(4-8*eta)*(*(*nodes)[5])[0];
	norm = std::sqrt(nx*nx+ny*ny);
	faceNormals(7,0) = nx/norm;
	faceNormals(7,1) = ny/norm;

	//Switch the orientation of the normals to make them inward normals:
        //faceNormals *= -1.;


	/* 
     for (CFuint iFace = 0; iFace < nbFacesInCell; ++iFace) {
       CFreal area = 0.;
       for (CFuint iDim = 0; iDim < dim; ++iDim) {
   area += faceNormals(iFace,iDim)*faceNormals(iFace, iDim);
       }
       faceAreas[iFace] = sqrt(area);
     }
	*/
	//*********************************************************************************************************
	//
	//				Compute lengths of edges numerically:
	//
	//*********************************************************************************************************

	CFreal weights[3] = { 5./36., 8./36., 5./36. };
	CFreal s = std::sqrt(0.6);
	CFreal qpPos[3] = { -0.25*s, 0.0, 0.25*s }; //Position of quadrature points 
						    //Integration domain <-0.25;0.25> considered
	
	CFreal dXdxi, dYdxi;
	CFreal shift[2] = { 0.25, 0.75 };	//Shift coordinates to have interval
						// a) <0.0 ; 0.5>
						// b) <0.5 ; 1.0>
						//The reference triangle has edge defined for xi = <0.0; 1.0>
						//We need two normals for this edge for P2 triangle, therefore two 
						//intervals

	CFreal ipXi, ipEta;	//Coordinates of Gauss points in reference space

	vector<CFuint> edgeIdx;
	edgeIdx.resize(2);
	for(CFuint iEdge = 0; iEdge < nbFacesInCell; ++iEdge) faceAreas[iEdge] = 0.0;

	//Length of edges 2 and 5

	edgeIdx[0] = 2;
	edgeIdx[1] = 5;

	for(CFuint iEdge = 0; iEdge < 2; ++iEdge) {
		for(CFuint iq = 0; iq <3; ++iq) {

			ipXi = qpPos[iq]+shift[iEdge];

			dXdxi =(-3+4*ipXi)*(*(*nodes)[0])[0]+\
			(4*ipXi-1)*(*(*nodes)[1])[0]+(4-8*ipXi)*(*(*nodes)[3])[0];
			dYdxi =(-3+4*ipXi)*(*(*nodes)[0])[1]+\
			(4*ipXi-1)*(*(*nodes)[1])[1]+(4-8*ipXi)*(*(*nodes)[3])[1];

			faceAreas[edgeIdx[iEdge]] += weights[iq] * std::sqrt(dXdxi*dXdxi+dYdxi*dYdxi);
		}
	}

	//Length of edges 3 and 6

	edgeIdx[0] = 3;
	edgeIdx[1] = 6;

	for(CFint iEdge = 1; iEdge >= 0; --iEdge) {
		for(CFuint iq = 0; iq <3; ++iq) {

			ipXi = qpPos[iq]+shift[iEdge];

			dXdxi = (4*ipXi-1)*(*(*nodes)[1])[0]+\
			(-3+4*ipXi)*(*(*nodes)[2])[0]+(4-8*ipXi)*(*(*nodes)[4])[0];	
			dYdxi = (4*ipXi-1)*(*(*nodes)[1])[1]+\
			(-3+4*ipXi)*(*(*nodes)[2])[1]+(4-8*ipXi)*(*(*nodes)[4])[1];

			faceAreas[edgeIdx[iEdge]] += weights[iq] * std::sqrt(dXdxi*dXdxi+dYdxi*dYdxi);
		}
	}

	//Length of edges 1 and 7

	edgeIdx[0] = 1;
	edgeIdx[1] = 7;

	for(CFint iEdge = 1; iEdge >= 0; --iEdge) {
		for(CFuint iq = 0; iq <3; ++iq) {

			ipXi = qpPos[iq]+shift[iEdge];

			dXdxi = (-3+4*ipXi)*(*(*nodes)[0])[0]+\
			(4*ipXi-1)*(*(*nodes)[2])[0]+(4-8*ipXi)*(*(*nodes)[5])[0];	
			dYdxi = (-3+4*ipXi)*(*(*nodes)[0])[1]+\
			(4*ipXi-1)*(*(*nodes)[2])[1]+(4-8*ipXi)*(*(*nodes)[5])[1];

			faceAreas[edgeIdx[iEdge]] += weights[iq] * std::sqrt(dXdxi*dXdxi+dYdxi*dYdxi);
		}
	}


	//Inner edge connecting nodes 3 and 4

  	for(CFuint iq = 0; iq <3; ++iq) {

 	ipEta = qpPos[iq]+0.25;

	dXdxi = (4*ipEta-1)*(*(*nodes)[0])[0]+(4*ipEta-1)*(*(*nodes)[2])[0]-2*(*(*nodes)[3])[0]+\
						2*(*(*nodes)[4])[0]+(2-8*ipEta)*(*(*nodes)[5])[0];
	dYdxi = (4*ipEta-1)*(*(*nodes)[0])[1]+(4*ipEta-1)*(*(*nodes)[2])[1]-2*(*(*nodes)[3])[1]+\
						2*(*(*nodes)[4])[1]+(2-8*ipEta)*(*(*nodes)[5])[1];

	faceAreas[4] += weights[iq] * std::sqrt(dXdxi*dXdxi+dYdxi*dYdxi);

 	}

	//Inner edge connecting nodes 4 and 5

  	for(CFuint iq = 0; iq <3; ++iq) {

 	ipXi = qpPos[iq]+0.25;

	dXdxi = (4*ipXi-1)*(*(*nodes)[0])[0]+(4*ipXi-1)*(*(*nodes)[1])[0]+(2-8*ipXi)*(*(*nodes)[3])[0]+\
				2*(*(*nodes)[4])[0]-2*(*(*nodes)[5])[0];
	dYdxi = (4*ipXi-1)*(*(*nodes)[0])[1]+(4*ipXi-1)*(*(*nodes)[1])[1]+(2-8*ipXi)*(*(*nodes)[3])[1]+\
				2*(*(*nodes)[4])[1]-2*(*(*nodes)[5])[1];

	faceAreas[8] += weights[iq] * std::sqrt(dXdxi*dXdxi+dYdxi*dYdxi);

 	}

	//Inner edge connecting nodes 5 and 3

  	for(CFuint iq = 0; iq <3; ++iq) {

 	ipXi = qpPos[iq]+0.25;

	dXdxi = (4*ipXi-1)*(*(*nodes)[1])[0]+(4*ipXi-1)*(*(*nodes)[2])[0]+2*(*(*nodes)[3])[0]+\
			(2-8*ipXi)*(*(*nodes)[4])[0]-2*(*(*nodes)[5])[0];
	dYdxi = (4*ipXi-1)*(*(*nodes)[1])[1]+(4*ipXi-1)*(*(*nodes)[2])[1]+2*(*(*nodes)[3])[1]+\
			(2-8*ipXi)*(*(*nodes)[4])[1]-2*(*(*nodes)[5])[1];

	faceAreas[0] += weights[iq] * std::sqrt(dXdxi*dXdxi+dYdxi*dYdxi);

 	}

	
	for(CFuint iNormal=0;iNormal<9;++iNormal)
	{
	faceNormals(iNormal,XX) *= faceAreas[iNormal];
	faceNormals(iNormal,YY) *= faceAreas[iNormal];	
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

void ComputeInwardNormalsTriagP2::update(const CFuint& iFirstCell,
                                         const CFuint& iLastCell,
                                         const CFuint& iType)
{
  throw Common::NotImplementedException (FromHere(),"ComputeInwardNormalsTriagP2::update()");
}

//////////////////////////////////////////////////////////////////////////////

void ComputeInwardNormalsTriagP2::average(const CFuint& iFirstCell,
                                          const CFuint& iLastCell,
                                          const CFuint& iType)
{
  throw Common::NotImplementedException (FromHere(),"ComputeInwardNormalsTriagP2::average()");
}

//////////////////////////////////////////////////////////////////////////////


    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
