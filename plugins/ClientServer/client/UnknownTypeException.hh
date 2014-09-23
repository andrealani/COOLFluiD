#ifndef COOLFluiD_client_UnknownTypeException_h
#define COOLFluiD_client_UnknownTypeException_h

/////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"

namespace COOLFluiD
{
  
  namespace client
  {
    
    /////////////////////////////////////////////////////////////////////////////
    
    class UnknownTypeException : public COOLFluiD::Common::Exception
    {
    public:
      
      /// Constructor
      UnknownTypeException(const COOLFluiD::Common::CodeLocation& where, 
                           const std::string& what);
      
      /// Copy constructor
      UnknownTypeException(const UnknownTypeException& e) throw ();
      
      
    }; // class UnknownTypeException
    
    /////////////////////////////////////////////////////////////////////////////
    
  } // namespace client
} // namespace COOLFluiD

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_client_UnknownTypeException_h
