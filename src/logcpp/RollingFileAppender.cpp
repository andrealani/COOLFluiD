/*
 * RollingFileAppender.cpp
 *
 * See the COPYING file for the terms of usage and distribution.
 */

#include "PortabilityImpl.hh"
#ifdef LOGCPP_HAVE_IO_H
#    include <io.h>
#endif
#ifdef CF_HAVE_UNISTD_H
#    include <unistd.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <logcpp/RollingFileAppender.hh>
#include <logcpp/Category.hh>
#include <sstream>

namespace logcpp {

    RollingFileAppender::RollingFileAppender(const std::string& name,
                                             const std::string& fileName,
                                             size_t maxFileSize,
                                             unsigned int maxBackupIndex,
                                             bool append,
                                             mode_t mode) :
        FileAppender(name, fileName, append, mode),
        _maxBackupIndex(maxBackupIndex),
        _maxFileSize(maxFileSize) {
    }

    void RollingFileAppender::setMaxBackupIndex(unsigned int maxBackups) {
        _maxBackupIndex = maxBackups;
    }

    unsigned int RollingFileAppender::getMaxBackupIndex() const {
        return _maxBackupIndex;
    }

    void RollingFileAppender::setMaximumFileSize(size_t maxFileSize) {
        _maxFileSize = maxFileSize;
    }

    size_t RollingFileAppender::getMaxFileSize() const {
        return _maxFileSize;
    }

    void RollingFileAppender::rollOver() {
        ::close(_fd);
        if (_maxBackupIndex > 0) {
            std::ostringstream oldName;
            oldName << _fileName << "." << _maxBackupIndex << std::ends;
            ::remove(oldName.str().c_str());
            size_t n = _fileName.length() + 1;
            for(unsigned int i = _maxBackupIndex; i > 1; i--) {
            	std::string newName = oldName.str();
                oldName.seekp(n, std::ios::beg);
                oldName << i-1 << std::ends;
                ::rename(oldName.str().c_str(), newName.c_str());
            }
            ::rename(_fileName.c_str(), oldName.str().c_str());
        }
        _fd = ::open(_fileName.c_str(), _flags, _mode);
    }

    void RollingFileAppender::_append(const LoggingEvent& event) {
        FileAppender::_append(event);
        off_t offset = ::lseek(_fd, 0, SEEK_END);
        if (offset < 0) {
            // XXX we got an error, ignore for now
        } else {
            if(static_cast<size_t>(offset) >= _maxFileSize) {
                rollOver();
            }
        }
    }
}
