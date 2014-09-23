#ifndef COOLFluiD_server_RecievingProcess_h
#define COOLFluiD_server_RecievingProcess_h

/////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD
{
  namespace server
  {
    
    /////////////////////////////////////////////////////////////////////////////
    
    //   struct ReceivedFrameInfo
    //   {
    //    COOLFluiD::Config::BuilderParserFrameInfo frame;
    
    //   };
    
    /// @todo this class should be removed
    
    class MPIReceiver
    {
    public:
      
      static void receive(COOLFluiD::Config::BuilderParserFrameInfo & frame, 
                          MPI::Intercomm comm,
                          const COOLFluiD::Config::BuilderParserRules & rules);
    }; // class ReceivingProcess
    
    /////////////////////////////////////////////////////////////////////////////
    
  } // namespace server
} // namespace COOLFluiD

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_server_RecievingProcess_h