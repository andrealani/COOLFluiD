#include <mpi.h>

#include "Config/BuilderParser.hh"
#include "Config/BuilderParserRules.hh"
#include "Config/BuilderParserFrameInfo.hh"

#include "ClientServer/server/MPIReceiver.hh"

using namespace MPI;

using namespace COOLFluiD::Config;
using namespace COOLFluiD::server;

void MPIReceiver::receive(BuilderParserFrameInfo & frame, Intercomm comm,
                          const BuilderParserRules & rules)
{
  char data[65536];
  frame.clear();
  
  comm.Recv(data, sizeof(data), CHAR, ANY_SOURCE, ANY_TAG);
  
  BuilderParser::parseFrame(data, rules, frame);
}
