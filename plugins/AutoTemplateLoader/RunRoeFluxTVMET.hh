#include "BaseRunLib.hh"

#include <tvmet/Vector.h>
#include <tvmet/Matrix.h>

template < typename TYPE , int SIZE >
class RunRoeFluxTVMET : public BaseRunLib
{
  public:

  typedef tvmet::Vector<TYPE,SIZE>      DEVector;
  typedef tvmet::Matrix<TYPE,SIZE,SIZE> GEMatrix;

  RunRoeFluxTVMET () : rEV(), lEV(), eV(), reV(), leV(), abseV(), absJ(), rflux(), lflux()
  {
  };

  protected: // helper functions

    void compute()
  {
      for(int i = 0; i < nbelems; ++i)
      {
        DEVector& rstate = *(R[i]);
        DEVector& lstate = *(L[i]);
        DEVector& flux   = *(RHS[i]);

        const TYPE nx    = 0.707106;
        const TYPE ny    = 0.707106;

        const TYPE rinvrho = 1.0 / rstate[0];
        const TYPE linvrho = 1.0 / lstate[0];
        const TYPE avRho = 0.5*( sqrt(rstate[0]) + sqrt(lstate[0]) );
        const TYPE avU   = 0.5*( rstate[1]*rinvrho + lstate[1]*linvrho);
        const TYPE avV   = 0.5*( rstate[2]*rinvrho + lstate[2]*linvrho);
        const TYPE avH   = 326548.18;
        const TYPE avA   = 340.29;

        const TYPE gamma = 1.4;
        const TYPE gammaMinus1 = gamma - 1.;
        const TYPE um = avU*nx + avV*ny;
        const TYPE ra = 0.5*avRho/avA;
        const TYPE avA2 = avA*avA;
        const TYPE coeffM2 = 0.5*gammaMinus1*(avU*avU + avV*avV)/avA2;
        const TYPE invAvRho = 1./avRho;
        const TYPE uDivA = gammaMinus1*avU/avA;
        const TYPE vDivA = gammaMinus1*avV/avA;
        const TYPE rhoA =  avRho*avA;

        rEV(0,0) = 1.;
        rEV(0,1) = 0.;
        rEV(0,2) = ra;
        rEV(0,3) = ra;
        rEV(1,0) = avU;
        rEV(1,1) = avRho*ny;
        rEV(1,2) = ra*(avU + avA*nx);
        rEV(1,3) = ra*(avU - avA*nx);
        rEV(2,0) = avV;
        rEV(2,1) = -avRho*nx;;
        rEV(2,2) = ra*(avV + avA*ny);
        rEV(2,3) = ra*(avV - avA*ny);
        rEV(3,0) = 0.5*(avU*avU +avV*avV);
        rEV(3,1) = avRho*(avU*ny - avV*nx);
        rEV(3,2) = ra*(avH + avA*um);
        rEV(3,3) = ra*(avH - avA*um);

        lEV(0,0) = 1.- coeffM2;
        lEV(0,1) = uDivA/avA;
        lEV(0,2) = vDivA/avA;
        lEV(0,3) = -gammaMinus1/avA2;
        lEV(1,0) = invAvRho*(avV*nx - avU*ny);
        lEV(1,1) = invAvRho*ny;
        lEV(1,2) = -invAvRho*nx;
        lEV(1,3) = 0.0;
        lEV(2,0) = avA*invAvRho*(coeffM2 - um/avA);
        lEV(2,1) = invAvRho*(nx - uDivA);
        lEV(2,2) = invAvRho*(ny - vDivA);
        lEV(2,3) = gammaMinus1/rhoA;
        lEV(3,0) = avA*invAvRho*(coeffM2 + um/avA);
        lEV(3,1) = -invAvRho*(nx + uDivA);
        lEV(3,2) = -invAvRho*(ny + vDivA);
        lEV(3,3) = gammaMinus1/rhoA;

        abseV[0] = um;
        abseV[1] = um;
        abseV[2] = um + avA;
        abseV[3] = um - avA;

//         absJ = rEV * lEV;

        rflux = 10.0;
        lflux = 20.0;

        flux = 0.5 * ( (rflux + lflux) - prod(prod(rEV,lEV),(rstate - lstate)) );

// works
//         std::cout << flux << std::endl;
      }
  }

  void init()
  {
      R = new DEVector*[nbelems];
      L = new DEVector*[nbelems];
      RHS = new DEVector*[nbelems];
      for(int i = 0; i < nbelems; ++i)
      {
          R[i] = new DEVector(SIZE);
          (*(R[i]))[0] = 1.22503;
          (*(R[i]))[1] = 333.412;
          (*(R[i]))[2] = 7.27509;
          (*(R[i]))[3] = 298706.1;

          L[i] = new DEVector(SIZE);
          (*(L[i]))[0] = 1.22503;
          (*(L[i]))[1] = 333.492;
          (*(L[i]))[2] = 0.0;
          (*(L[i]))[3] = 298706.1;

          RHS[i] = new DEVector(SIZE);
      }
  }

  void finalize()
  {
      for(int i = 0; i < nbelems; ++i)
      {
          delete R[i];
          delete L[i];
          delete RHS[i];
      }
      delete [] R;
      delete [] L;
      delete [] RHS;
  }

  private: // data

    DEVector ** R;
    DEVector ** L;
    DEVector ** RHS;

  /// matrix of right eigenvectors
  GEMatrix   rEV;

  /// matrix of left eigenvectors
  GEMatrix   lEV;

  /// vector of eigenvalues
  DEVector   eV;

  /// vector of right state eigenvalues
  DEVector   reV;

  /// vector of left state eigenvalues
  DEVector   leV;

  /// vector of eigenvalues
  DEVector   abseV;

  /// abs of the jacobian matrix
  GEMatrix   absJ;

  /// vector of right state eigenvalues
  DEVector   rflux;

  /// vector of left state eigenvalues
  DEVector   lflux;

};
