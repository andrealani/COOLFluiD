#include "BaseRunLib.hh"
#include <flens/flens.h>

template < typename TYPE , int SIZE >
class RunFlens : public BaseRunLib
{
  public:

  typedef flens::GeMatrix<flens::FullStorage<TYPE, flens::ColMajor> >   GEMatrix;
  typedef flens::DenseVector<flens::Array<TYPE> >         DEVector;

  RunFlens () {};

  protected:

  void compute()
  {
using namespace flens;
using namespace std;

      for( unsigned int i = 0; i < nbelems; ++i)
      {
        GEMatrix& a = *A[i];
        DEVector& x = *X[i];
        DEVector& y = *Y[i];
        DEVector& z = *Z[i];

        z = 2*a*x + 1.5*y;

#if 0
        cout << "A = " << a << endl;
        cout << "x = " << x << endl;
        cout << "y = " << y << endl;
        cout << "z = 2*A*x + 1.5*y = " << z << endl;
#endif
      }
  }

  void init()
  {
      A = new GEMatrix*[nbelems];
      for(int i = 0; i < nbelems; ++i)
      {
          A[i] = new GEMatrix(SIZE,SIZE);
          init(*A[i]);
      }

      X = new DEVector*[nbelems];
      Y = new DEVector*[nbelems];
      Z = new DEVector*[nbelems];
      for(int i = 0; i < nbelems; ++i)
      {
          X[i] = new DEVector(SIZE);
          Y[i] = new DEVector(SIZE);
          Z[i] = new DEVector(SIZE);
          init(*X[i]);
          init(*Y[i]);
      }
  }

  void finalize()
  {
      for(int i = 0; i < nbelems; ++i)
          delete A[i];
      delete [] A;

      for(int i = 0; i < nbelems; ++i)
      {
          delete X[i];
          delete Y[i];
          delete Z[i];
      }
      delete [] X;
      delete [] Y;
      delete [] Z;
  }

  private: // helper functions
  void init (GEMatrix& A)
  {
    for (int i=1; i<=SIZE; ++i){
      for (int j=1; j<=SIZE ; ++j){ A(i,j)= i*j; }
    }
  }

  void init (DEVector& v)
  {
    for (int j=1; j<=SIZE ; ++j){ v(j) = j; }
  }

  private: // data

    GEMatrix ** A;
    DEVector ** X;
    DEVector ** Y;
    DEVector ** Z;

};
