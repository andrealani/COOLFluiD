#include <iomanip>

#include "MatrixStabilityMethodWriter/MatrixStabilityMethodWriter.hh"
#include "MatrixStabilityMethodWriter/AddMatrixColumnToFile.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MatrixStabilityMethodWriter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<AddMatrixColumnToFile,
                      MatrixStabilityMethodData,
                      MatrixStabilityMethodWriterModule
                     >
addMatrixColumnToFileProvider("AddMatrixColumnToFile");

//////////////////////////////////////////////////////////////////////////////

AddMatrixColumnToFile::AddMatrixColumnToFile(const std::string& name) :
    MatrixStabilityMethodCom(name),
    socket_rhs("rhs")
{
}

//////////////////////////////////////////////////////////////////////////////

void AddMatrixColumnToFile::execute()
{
  DataHandle< CFreal> rhs = socket_rhs.getDataHandle();

  SafePtr< ofstream > outputFile = getMethodData().getOutputFile();

  const CFuint nbrStates = getMethodData().getNbrStates();
  if (getMethodData().writeBinary())
  {
    for (CFuint i = 0; i < nbrStates; ++i)
    {
      outputFile->write(reinterpret_cast< char* >(&rhs[i]),sizeof(CFreal));
    }
  }
  else
  {
    for (CFuint i = 0; i < nbrStates; ++i)
    {
      *outputFile << setprecision(12) << rhs[i] << " ";
    }
    *outputFile << "\n";
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > AddMatrixColumnToFile::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MatrixStabilityMethodWriter

  } // namespace Numerics

} // namespace COOLFluiD
