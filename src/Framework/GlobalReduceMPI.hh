// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_GlobalReduceMPI_hh
#define COOLFluiD_Framework_GlobalReduceMPI_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "Common/PE.hh"
#include "Common/NonCopyable.hh"
#include "Common/MPI/MPIInitObject.hh"
#include "Common/MPI/MPIException.hh"
#include "Common/MPI/MPIStructDef.hh"
#include "Framework/MeshData.hh"
#include "Framework/GlobalReduce.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// The moment you instantiate this class the PE should be initialised
/// (and it should be freed before the PE is stopped)!
/// Warning: not thread safe
/// Warning: RESULTTYPE should not contain virtual methods
/// The user should provide 2 methods and 1 typedef:
/// typedef:
///    GR_RESULTTYPE : this should be the type of the data to be reduce
/// methods:
///    (static): a combine method, taking 2 input and 1 output parameter
///    void GR_Combine (const GR_RESULTTYPE & S1,
///       const GR_RESULTTYPE & S2, GR_RESULTTYPE & Out);
///    (nonstatic): a method for returning the local result
///    GR_RESULTTYPE GR_GetLocalValue () const
template <typename PROVIDER>
class GlobalReduce<PROVIDER, Common::PM_MPI> : public Common::NonCopyable <GlobalReduce<PROVIDER, Common::PM_MPI> >
{
public:
  
  typedef typename PROVIDER::GR_RESULTTYPE RESULTTYPE;
  
  GlobalReduce (const PROVIDER & P);
  ~GlobalReduce ();
  
  /// Return the accumulated global result
  /// WARNING: Global!!
  RESULTTYPE GetGlobalValue (const std::string& nspaceName) const;
  
private: // functions

  /// The actual function called by MPI
  static void CombineMPI (void * A, void * B, int * len, MPI_Datatype *
			  Datatype);
  
  /// Register our function
  void Register ();
  
  /// Unregister the function
  void UnRegister ();
  
private: // variables
  
  /// Forward
  class RegisterHelper;
  
  const PROVIDER & Provider_;
  
  static RegisterHelper * Helper_;
  
  /// Nested helper class
  class RegisterHelper : public Common::NonCopyable<RegisterHelper>, public Common::MPIInitObject
  {
  public:
    RegisterHelper () : UserCount_(0), OperationHandle_(MPI_OP_NULL),  Registered_(false) {}
    
    virtual ~RegisterHelper ()
    {
      if (Registered_)
	throw Common::MPIException (FromHere(), "GlobalReduce: RegisterHelper: "
				    "destructor called before MPI_Done was called on "
				    "this object!\n");
      if (UserCount_)
	throw Common::MPIException (FromHere(), "GlobalReduce: RegisterHelper: "
				    "destructor called on this object while our parent "
				    "GlobalReduce object is still using the function!\n");
    }
    
    virtual void MPI_Init (MPI_Comm  Communicator)
    {
      // For now always assume commutating operations
      MPI_Op_create (CombineMPI, true, &OperationHandle_);
      Registered_=true;
    }
    
    virtual void MPI_Done ()
    {
      if (!Registered_)
	throw Common::MPIException (FromHere(), "MPI_Done called before MPI_Init!!\n");
      if (UserCount_)
	throw Common::MPIException (FromHere(), "GlobalReduce: RegisterHelper: "
				    "destructor called on this object while our parent "
				    "GlobalReduce object is still using the function!\n");
      
      MPI_Op_free (&OperationHandle_);
      Registered_=false;
    }
    
    bool IsRegistered () const  {  return Registered_; }
    
    void RegisterUser () { ++UserCount_; }
    
    void UnRegisterUser ()
    {
      cf_assert (UserCount_>0);
      --UserCount_;
    }
    
    MPI_Op GetOperation ()
    {
      if (!Registered_)
	throw Common::MPIException
	  (FromHere(), "GlobalReduce::RegisterHelper: GetOperation called"
	   " before the PE called MPI_Init!\n");
      return OperationHandle_;
    }
    
    
  private:
    unsigned int UserCount_;
    MPI_Op OperationHandle_;
    bool Registered_;
    
  }; // class RegisterHelper
  
}; // class GlobalReduce

//////////////////////////////////////////////////////////////////////////////

template <typename PROVIDER>
typename GlobalReduce<PROVIDER,Common::PM_MPI>::RegisterHelper *
  GlobalReduce<PROVIDER, Common::PM_MPI>::Helper_;

template <typename PROVIDER>
typename GlobalReduce<PROVIDER, Common::PM_MPI>::RESULTTYPE
GlobalReduce<PROVIDER, Common::PM_MPI>::GetGlobalValue (const std::string& nspaceName) const
{
  RESULTTYPE SendBuf = Provider_.GR_GetLocalValue ();
  RESULTTYPE ReceiveBuf;
  MPI_Allreduce (&SendBuf, &ReceiveBuf, 1,
		 Common::MPIStructDef::getMPIType(&SendBuf),
		 Helper_->GetOperation (), 
		 Common::PE::GetPE().GetCommunicator(nspaceName));
  return ReceiveBuf;
}
   
//////////////////////////////////////////////////////////////////////////////
   
template <typename PROVIDER>
void GlobalReduce<PROVIDER, Common::PM_MPI>::CombineMPI (void * a, void * b, int * len, MPI_Datatype * datatype)
{
    cf_assert (len != CFNULL);
    cf_assert (*len==1);

    const RESULTTYPE & S1 = *(static_cast<RESULTTYPE *>(a));
    RESULTTYPE & S2 = * (static_cast<RESULTTYPE *>(b));

    RESULTTYPE R;

    PROVIDER::GR_Combine(S1,S2,R);
    S2 = R;
}

//////////////////////////////////////////////////////////////////////////////
   
template <typename PROVIDER>
GlobalReduce<PROVIDER, Common::PM_MPI>::GlobalReduce (const PROVIDER & P) : Provider_(P)
{
    Register ();
}

template <typename PROVIDER>
GlobalReduce<PROVIDER, Common::PM_MPI>::~GlobalReduce ()
{
    UnRegister ();
}

//////////////////////////////////////////////////////////////////////////////
   
template <typename PROVIDER>
void GlobalReduce<PROVIDER, Common::PM_MPI>::Register ()
{
    // We never free the helper. It outlives it's parent object
    if (!Helper_)
    {
      Helper_ = new RegisterHelper ();
      Common::PE::GetPE().RegisterInitObject (Helper_);
    }

    // We already have a helper
    Helper_->RegisterUser ();
}

//////////////////////////////////////////////////////////////////////////////
   
template <typename PROVIDER>
void GlobalReduce<PROVIDER, Common::PM_MPI>::UnRegister ()
{
    cf_assert (Helper_ != CFNULL );
    Helper_->UnRegisterUser ();
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAGCLASS, typename BASETYPE>
class GR_CombineHelper;

//////////////////////////////////////////////////////////////////////////////

template <typename BASETYPE>
class GR_CombineHelper<GRO_SUM,BASETYPE>
{
public:
    static void DoOperation (const BASETYPE & S1, const BASETYPE & S2, BASETYPE & Out)
    {
      Out = S1+S2;
    }
};

//////////////////////////////////////////////////////////////////////////////

template <typename BASETYPE>
class GR_CombineHelper<GRO_MAX,BASETYPE>
{
public:
    static void DoOperation (const BASETYPE & S1, const BASETYPE & S2, BASETYPE & Out)
    {
       Out = std::max(S1,S2);
    }
};

//////////////////////////////////////////////////////////////////////////////

template <typename BASETYPE>
class GR_CombineHelper<GRO_MIN,BASETYPE>
{
public:
    static void DoOperation (const BASETYPE & S1, const BASETYPE & S2, BASETYPE & Out)
    {
       Out = std::min(S1,S2);
    }
};

//////////////////////////////////////////////////////////////////////////////

template <typename TAGCLASS, typename BASETYPE>
class GlobalReduceOperationMPIHelper
{
public:
    typedef BASETYPE GR_RESULTTYPE;

    static void GR_Combine (const GR_RESULTTYPE & S1, const GR_RESULTTYPE & S2, GR_RESULTTYPE & Out)
    {
      GR_CombineHelper<TAGCLASS, GR_RESULTTYPE>::DoOperation(S1, S2, Out);
    }

    GR_RESULTTYPE GR_GetLocalValue () const  {  return *Source_; }

    GlobalReduceOperationMPIHelper (BASETYPE * S, BASETYPE * D, unsigned int Count)
      : Source_(S), Dest_(D), Count_(Count), Reduce_(*this) {}
  
  void GetGlobalValue (BASETYPE * Dest) 
  {
    const std::string nsp = Framework::MeshDataStack::getActive()->getPrimaryNamespace();
    *Dest = Reduce_.GetGlobalValue (nsp); 
  }
  
private:
    const BASETYPE * Source_;
    BASETYPE * Dest_;
    unsigned int Count_;

    GlobalReduce < GlobalReduceOperationMPIHelper <TAGCLASS, BASETYPE> , Common::PM_MPI > Reduce_;

};

//////////////////////////////////////////////////////////////////////////////

/// Default (non MPI specialized)
/// GlobalReduceOperation for PE PM_MPI
template <typename TAGCLASS, typename BASETYPE>
class GlobalReduceOperationHelper<TAGCLASS, BASETYPE, Common::PM_MPI>
{
public:
    static inline void DoReduce (BASETYPE * S, BASETYPE * D, unsigned int Count)
    {
      GlobalReduceOperationMPIHelper <TAGCLASS, BASETYPE> Helper (S,D,Count);
      Helper.GetGlobalValue (D);
    }
};

//////////////////////////////////////////////////////////////////////////////

    }
}

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_GlobalReduceMPI_hh
