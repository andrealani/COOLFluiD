/*
 * PatternLayout.cpp
 *
 * Copyright 2002, Bastiaan Bakker. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */

#include "PortabilityImpl.hh"

#include <logcpp/PatternLayout.hh>
#include <logcpp/Priority.hh>
#include <logcpp/NDC.hh>
#include <logcpp/TimeStamp.hh>

#include <sstream>

#include <iomanip>
#include <ctime>
#include <cmath>

namespace logcpp {

    struct StringLiteralComponent : public PatternLayout::PatternComponent {
        StringLiteralComponent(const std::string& literal) :
            _literal(literal) {
        }

        virtual void append(std::ostringstream& out, const LoggingEvent& event) {
            out << _literal;
        }

        private:
        std::string _literal;
    };

    struct CategoryNameComponent : public PatternLayout::PatternComponent {
        CategoryNameComponent(std::string specifier) {
            if (specifier == "") {
                _precision = -1;
            } else {
                std::istringstream s(specifier);
                s >> _precision;
            }
        }

        virtual void append(std::ostringstream& out, const LoggingEvent& event) {
            if (_precision == -1) {
                out << event.categoryName;
            } else {
                std::string::size_type begin = std::string::npos;
                for(int i = 0; i < _precision; i++) {
                    begin = event.categoryName.rfind('.', begin - 2);
                    if (begin == std::string::npos) {
                        begin = 0;
                        break;
                    }
                    begin++;
                }
                out << event.categoryName.substr(begin);
            }
        }

        private:
        int _precision;
    };

    struct MessageComponent : public PatternLayout::PatternComponent {
        virtual void append(std::ostringstream& out, const LoggingEvent& event) {
            out << event.message;
        }
    };

    struct NDCComponent : public PatternLayout::PatternComponent {
        virtual void append(std::ostringstream& out, const LoggingEvent& event) {
            out << event.ndc;
        }
    };

    struct PriorityComponent : public PatternLayout::PatternComponent {
        virtual void append(std::ostringstream& out, const LoggingEvent& event) {
            out << Priority::getPriorityName(event.priority);
        }
    };

    struct ThreadNameComponent : public PatternLayout::PatternComponent {
        virtual void append(std::ostringstream& out, const LoggingEvent& event) {
            out << event.threadName;
        }
    };

    struct ProcessorTimeComponent : public PatternLayout::PatternComponent {
        virtual void append(std::ostringstream& out, const LoggingEvent& event) {
            out << ::clock();
        }
    };

    struct TimeStampComponent : public PatternLayout::PatternComponent {
        static const char* const FORMAT_ISO8601;
        static const char* const FORMAT_ABSOLUTE;
        static const char* const FORMAT_DATE;

        TimeStampComponent(std::string timeFormat) {
            if ((timeFormat == "") || (timeFormat == "ISO8601")) {
                timeFormat = FORMAT_ISO8601;
            } else if (timeFormat == "ABSOLUTE") {
                timeFormat = FORMAT_ABSOLUTE;
            } else if (timeFormat == "DATE") {
                timeFormat = FORMAT_DATE;
            }
            std::string::size_type pos = timeFormat.find("%l");
            if (pos == std::string::npos) {
                _printMillis = false;
                _timeFormat1 = timeFormat;
            } else {
                _printMillis = true;
                _timeFormat1 = timeFormat.substr(0, pos);
                _timeFormat2 = timeFormat.substr(pos + 2);
            }
        }

        virtual void append(std::ostringstream& out, const LoggingEvent& event) {
            struct tm *currentTime;
            time_t t = event.timeStamp.getSeconds();
            currentTime = std::localtime(&t);
            char formatted[100];
            std::string timeFormat;
            if (_printMillis) {
                std::ostringstream formatStream;
                formatStream << _timeFormat1
                             << std::setw(3) << std::setfill('0')
                             << event.timeStamp.getMilliSeconds()
                             << _timeFormat2;
                timeFormat = formatStream.str();
            } else {
                timeFormat = _timeFormat1;
            }
            std::strftime(formatted, sizeof(formatted), timeFormat.c_str(), currentTime);
            out << formatted;
        }

        private:
        std::string _timeFormat1;
        std::string _timeFormat2;
        bool _printMillis;
    };

    const char* const TimeStampComponent::FORMAT_ISO8601 = "%Y-%m-%d %H:%M:%S,%l";
    const char* const TimeStampComponent::FORMAT_ABSOLUTE = "%H:%M:%S,%l";
    const char* const TimeStampComponent::FORMAT_DATE = "%d %b %Y %H:%M:%S,%l";

    struct SecondsSinceEpochComponent : public PatternLayout::PatternComponent {
        virtual void append(std::ostringstream& out, const LoggingEvent& event) {
            out << event.timeStamp.getSeconds();
        }
    };

    struct MillisSinceEpochComponent : public PatternLayout::PatternComponent {
        virtual void append(std::ostringstream& out, const LoggingEvent& event) {
            double t = event.timeStamp.getSeconds() -
                TimeStamp::getStartTime().getSeconds();
            t *= 1000;
            t += event.timeStamp.getMilliSeconds() -
                TimeStamp::getStartTime().getMilliSeconds();

            out << std::setiosflags(std::ios::fixed) << std::setprecision(0) << t;
        }
    };

    struct FormatModifierComponent : public PatternLayout::PatternComponent {
        FormatModifierComponent(PatternLayout::PatternComponent* component,
                                int minWidth, int maxWidth) :
            _component(component) ,
            _minWidth(std::abs( (float) minWidth)),
            _maxWidth(maxWidth),
            _alignLeft(minWidth < 0) {
        }

        virtual ~FormatModifierComponent() {
            delete _component;
        }

        virtual void append(std::ostringstream& out, const LoggingEvent& event) {
            std::ostringstream s;
            _component->append(s, event);
            std::string msg = s.str();
            if (_maxWidth > 0) {
                msg.erase(_maxWidth);
            }
            int fillCount = _minWidth - (int) msg.length();
            if (fillCount > 0) {
                if (_alignLeft) {
                    out << msg << std::string(fillCount, ' ');
                } else {
                    out << std::string(fillCount, ' ') << msg;
                }
            } else {
                out << msg;
            }
        }

        private:
        PatternLayout::PatternComponent* _component;
        int _minWidth;
        int _maxWidth;
        bool _alignLeft;
    };

    const char* PatternLayout::DEFAULT_CONVERSION_PATTERN = "%m%n";
    const char* PatternLayout::SIMPLE_CONVERSION_PATTERN = "%p - %m%n";
    const char* PatternLayout::BASIC_CONVERSION_PATTERN = "%R %p %c %x: %m%n";
    const char* PatternLayout::TTCC_CONVERSION_PATTERN = "%r [%t] %p %c %x - %m%n";

    PatternLayout::PatternLayout() {
        try {
            setConversionPattern(DEFAULT_CONVERSION_PATTERN);
        } catch(ConfigureFailure& e) {
		}
    }

    PatternLayout::~PatternLayout() {
        clearConversionPattern();
    }

    void PatternLayout::clearConversionPattern() {
        for(ComponentVector::const_iterator i = _components.begin();
            i != _components.end(); ++i) {
            delete (*i);
        }
        _components.clear();
        _conversionPattern = "";
    }

    void PatternLayout::setConversionPattern(const std::string& conversionPattern) throw(ConfigureFailure) {
        std::istringstream conversionStream(conversionPattern);
        std::string literal;

        char ch;
        PatternLayout::PatternComponent* component = NULL;
        int minWidth = 0;
        int maxWidth = 0;
        clearConversionPattern();
        while (conversionStream.get(ch)) {
            if (ch == '%') {
                // readPrefix;
                {
                    char ch2;
                    conversionStream.get(ch2);
                    if ((ch2 == '-') || ((ch2 >= '0') && (ch2 <= '9'))) {
                        conversionStream.putback(ch2);
                        conversionStream >> minWidth;
                        conversionStream.get(ch2);
                    }
                    if (ch2 == '.') {
                        conversionStream >> maxWidth;
                    } else {
                        conversionStream.putback(ch2);
                    }
                }
                if (!conversionStream.get(ch)) {
                    std::ostringstream msg;
                    msg << "unterminated conversion specifier in '" << conversionPattern << "' at index " << conversionStream.tellg();
                    throw ConfigureFailure(msg.str());
                }
                std::string specPostfix = "";
                // read postfix
                {
                    char ch2;
                    if (conversionStream.get(ch2)) {
                        if (ch2 == '{') {
                            while(conversionStream.get(ch2) && (ch2 != '}'))
                                specPostfix += ch2;
                        } else {
                            conversionStream.putback(ch2);
                        }
                    }
                }
                switch (ch) {
                case '%':
                    literal += ch;
                    break;
                case 'm':
                    component = new MessageComponent();
                    break;
                case 'n':
                    {
                        std::ostringstream endline;
                        endline << std::endl;
                        literal += endline.str();
                    }
                    break;
                case 'c':
                    component = new CategoryNameComponent(specPostfix);
                    break;
                case 'd':
                    component = new TimeStampComponent(specPostfix);
                    break;
                case 'p':
                    component = new PriorityComponent();
                    break;
                case 'r':
                    component = new MillisSinceEpochComponent();
                    break;
                case 'R':
                    component = new SecondsSinceEpochComponent();
                    break;
                case 'u':
                    component = new ProcessorTimeComponent();
                    break;
                case 'x':
                    component = new NDCComponent();
                    break;
                default:
                    std::ostringstream msg;
                    msg << "unknown conversion specifier '" << ch << "' in '" << conversionPattern << "' at index " << conversionStream.tellg();
                    throw ConfigureFailure(msg.str());
                }
                if (component) {
                    if (!literal.empty()) {
                        _components.push_back(new StringLiteralComponent(literal));
                        literal = "";
                    }
                    if ((minWidth != 0) || (maxWidth != 0)) {
                        component = new FormatModifierComponent(component, minWidth, maxWidth);
                        minWidth = maxWidth = 0;
                    }
                    _components.push_back(component);
                    component = NULL;
                }
            } else {
                literal += ch;
            }
        }
        if (!literal.empty()) {
            _components.push_back(new StringLiteralComponent(literal));
        }

        _conversionPattern = conversionPattern;
    }

    std::string PatternLayout::getConversionPattern() const {
        return _conversionPattern;
    }

    std::string PatternLayout::format(const LoggingEvent& event) {
        std::ostringstream message;

        for(ComponentVector::const_iterator i = _components.begin();
            i != _components.end(); ++i) {
            (*i)->append(message, event);
        }

        return message.str();
    }
}
