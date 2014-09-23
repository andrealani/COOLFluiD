/*
 * RemoteSyslogAppender.cpp
 *
 * Copyright 2001, LifeLine Networks BV (www.lifeline.nl). All rights reserved.
 * Copyright 2001, Walter Stroebel. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */

#include "PortabilityImpl.hh"

#ifdef CF_HAVE_UNISTD_H
#    include <unistd.h>
#endif
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <logcpp/RemoteSyslogAppender.hh>
#ifdef WIN32
#include <winsock2.h>
#else
#include <netdb.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#endif

namespace logcpp {

    int RemoteSyslogAppender::toSyslogPriority(Priority::Value priority) {
        static int priorities[8] = { LOG_EMERG, LOG_ALERT, LOG_CRIT, LOG_ERR,
                                     LOG_WARNING, LOG_NOTICE, LOG_INFO, 
                                     LOG_DEBUG };
        int result;

        priority++;
        priority /= 100;

        if (priority < 0) {
            result = LOG_EMERG;
        } else if (priority > 7) {
            result = LOG_DEBUG;
        } else {
            result = priorities[priority];
        }

        return result;
    }
        

    RemoteSyslogAppender::RemoteSyslogAppender(const std::string& name, 
                                   const std::string& syslogName, 
				   const std::string& relayer,
                                   int facility,
				   int portNumber) : 
        LayoutAppender(name),
        _syslogName(syslogName),
	_relayer(relayer),
        _facility((facility == -1) ? LOG_USER : facility),
	_portNumber((portNumber == -1) ? 514 : portNumber),
	_socket (0),
	_ipAddr (0),
	_cludge (0)
    {
        open();
    }
    
    RemoteSyslogAppender::~RemoteSyslogAppender() {
        close();
#ifdef WIN32
	if (_cludge) {
	    // we started it, we end it.
	    WSACleanup ();
	}
#endif
    }

    void RemoteSyslogAppender::open() {
	if (!_ipAddr) {
	    struct hostent *pent = gethostbyname (_relayer.c_str ());
	    if (pent == NULL) {
#ifdef WIN32
		if (WSAGetLastError () == WSANOTINITIALISED) {
		    WSADATA wsaData;
		    int err;
 
		    err = WSAStartup (0x101, &wsaData );
		    if (err) {
                        // loglog("RemoteSyslogAppender: WSAStartup returned %d", err);
                        return; // fail silently
                    }
		    pent = gethostbyname (_relayer.c_str ());
		    _cludge = 1;
		} else {
		    // loglog("RemoteSyslogAppender: gethostbyname returned error");
                    return; // fail silently
		}
#endif
	    }
	    if (pent == NULL) {
		unsigned long ip = (unsigned long) inet_addr (_relayer.c_str ());
		pent = gethostbyaddr ((const char *) &ip, 4, AF_INET);
                if (pent == NULL) {
                    // loglog("RemoteSyslogAppender: failed to resolve host %s", _relayer.c_str());
                    return; // fail silently                    
                }
            }
	    _ipAddr = *((unsigned long *) pent->h_addr);
	}
	// Get a datagram socket.
	
	if ((_socket = socket(AF_INET, SOCK_DGRAM, 0)) < 0) {
            // loglog("RemoteSyslogAppender: failed to open socket");
            return; // fail silently                    
	}
    }

    void RemoteSyslogAppender::close() {
	if (_socket) {
#ifdef WIN32
	    closesocket (_socket);
#else
	    ::close (_socket);
#endif
	    _socket = 0;
	}
    }

    void RemoteSyslogAppender::_append(const LoggingEvent& event) {
	std::string message(_getLayout().format(event));
	int len = message.length() + 16;
	char *buf = new char [len];
        int priority = _facility + toSyslogPriority(event.priority);
	int len2 = sprintf (buf, "<%d>", priority);
	memcpy (buf + len2, message.data(), len - 16);
	sockaddr_in sain;
	sain.sin_family = AF_INET;
	sain.sin_port   = htons (_portNumber);
	// NO, do NOT use htonl on _ipAddr. Is already in network order.
        sain.sin_addr.s_addr = _ipAddr;
	len = len - 16 + len2;
	while (len > len2) {
	    if (len > 900) {
    		sendto (_socket, buf, 900, 0, (struct sockaddr *) &sain, sizeof (sain));
		std::memmove (buf + len2, buf + 900, len - 900 - len2);
		len -= (900 - len2);
		// note: we might need to sleep a bit here
	    } else {
		sendto (_socket, buf, len, 0, (struct sockaddr *) &sain, sizeof (sain));
		break;
	    }
	}
	delete buf;
    }

    bool RemoteSyslogAppender::reopen() {
        close();
        open();
        return true;
    }      
}
