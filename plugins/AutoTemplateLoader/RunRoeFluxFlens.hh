#include "BaseRunLib.hh"

#include <flens/flens.h>

template < typename TYPE , int SIZE >
class RunRoeFluxFlens : public BaseRunLib
{
  public:

  typedef flens::GeMatrix<flens::FullStorage<TYPE, flens::ColMajor> >   GEMatrix;
  typedef flens::DenseVector<flens::Array<TYPE> >         DEVector;

  RunRoeFluxFlens () :  tmpv(SIZE), rEV(SIZE,SIZE), lEV(SIZE,SIZE), eV(SIZE), reV(SIZE), leV(SIZE), abseV(SIZE), absJ(SIZE,SIZE), rflux(SIZE), lflux(SIZE)
  {
  };

  protected: // helper functions

    void compute()
  {
      for(unsigned int i = 0; i < nbelems; ++i)
      {

        DEVector& rstate = *(R[i]);
        DEVector& lstate = *(L[i]);
        DEVector& flux   = *(RHS[i]);


        const TYPE nx    = 0.707106;
        const TYPE ny    = 0.707106;

        const TYPE rinvrho = 1.0 / rstate(1);
        const TYPE linvrho = 1.0 / lstate(1);
        const TYPE avRho = 0.5*( sqrt(rstate(1)) + sqrt(lstate(1)) );
        const TYPE avU   = 0.5*( rstate(2)*rinvrho + lstate(2)*linvrho);
        const TYPE avV   = 0.5*( rstate(3)*rinvrho + lstate(3)*linvrho);
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


        rEV(1,1) = 1.;
        rEV(1,2) = 0.;
        rEV(1,3) = ra;
        rEV(1,4) = ra;
        rEV(2,1) = avU;
        rEV(2,2) = avRho*ny;
        rEV(2,3) = ra*(avU + avA*nx);
        rEV(2,4) = ra*(avU - avA*nx);
        rEV(3,1) = avV;
        rEV(3,2) = -avRho*nx;;
        rEV(3,3) = ra*(avV + avA*ny);
        rEV(3,4) = ra*(avV - avA*ny);
        rEV(4,1) = 0.5*(avU*avU +avV*avV);
        rEV(4,2) = avRho*(avU*ny - avV*nx);
        rEV(4,3) = ra*(avH + avA*um);
        rEV(4,4) = ra*(avH - avA*um);


        lEV(1,1) = 1.- coeffM2;
        lEV(1,2) = uDivA/avA;
        lEV(1,3) = vDivA/avA;
        lEV(1,4) = -gammaMinus1/avA2;
        lEV(2,1) = invAvRho*(avV*nx - avU*ny);
        lEV(2,2) = invAvRho*ny;
        lEV(2,3) = -invAvRho*nx;
        lEV(2,4) = 0.0;
        lEV(3,1) = avA*invAvRho*(coeffM2 - um/avA);
        lEV(3,2) = invAvRho*(nx - uDivA);
        lEV(3,3) = invAvRho*(ny - vDivA);
        lEV(3,4) = gammaMinus1/rhoA;
        lEV(4,1) = avA*invAvRho*(coeffM2 + um/avA);
        lEV(4,2) = -invAvRho*(nx + uDivA);
        lEV(4,3) = -invAvRho*(ny + vDivA);
        lEV(4,4) = gammaMinus1/rhoA;

        abseV(1) = um;
        abseV(2) = um;
        abseV(3) = um + avA;
        abseV(4) = um - avA;

        absJ = rEV * lEV;

        rflux = 10.0;
        lflux = 20.0;

        tmpv = (rstate - lstate);

        flux = 0.5 * ( (rflux + lflux)  - absJ * tmpv);

// works
//         std::cout << flux << std::endl;
      }
  }

  void init()
  {
      R = new DEVector*[nbelems];
      L = new DEVector*[nbelems];
      RHS = new DEVector*[nbelems];
      for(unsigned int i = 0; i < nbelems; ++i)
      {
          R[i] = new DEVector(SIZE);
          (*(R[i]))(1) = 1.22503;
          (*(R[i]))(2) = 333.412;
          (*(R[i]))(3) = 7.27509;
          (*(R[i]))(4) = 298706.1;

          L[i] = new DEVector(SIZE);
          (*(L[i]))(1) = 1.22503;
          (*(L[i]))(2) = 333.492;
          (*(L[i]))(3) = 0.0;
          (*(L[i]))(4) = 298706.1;

          RHS[i] = new DEVector(SIZE);
      }
  }

  void finalize()
  {
      for(unsigned int i = 0; i < nbelems; ++i)
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

  DEVector   tmpv;
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
