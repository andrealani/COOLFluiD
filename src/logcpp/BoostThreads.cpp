#include "logcpp/BoostThreads.hh"

namespace logcpp {
    namespace threading {

        std::string getThreadId()
        {
           std::ostringstream oss; 
           oss << boost::this_thread::get_id();
           return oss.str();
        }

    }
}
