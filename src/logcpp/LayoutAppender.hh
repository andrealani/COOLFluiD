// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/*
 * LayoutAppender.hh
 *
 * Copyright 2000, LifeLine Networks BV (www.lifeline.nl). All rights reserved.
 * Copyright 2000, Bastiaan Bakker. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */

#ifndef CF_LOGCPP_LAYOUTAPPENDER_HH
#define CF_LOGCPP_LAYOUTAPPENDER_HH

#include <string>
#include <logcpp/Portability.hh>
#include <logcpp/AppenderSkeleton.hh>
#include <logcpp/BasicLayout.hh>

namespace logcpp {

    /**
     * LayoutAppender is a common superclass for all Appenders that require
     * a Layout. 
     **/
    class logcpp_API LayoutAppender : public AppenderSkeleton {
        public:

        typedef BasicLayout DefaultLayoutType;

        LayoutAppender(const std::string& name);
        virtual ~LayoutAppender();
        
        /**
         * Check if the appender requires a layout. All LayoutAppenders do,
         * therefore this method returns true for all subclasses.
         * 
         * @returns true.
         **/
        virtual bool requiresLayout() const;
        virtual void setLayout(Layout* layout = NULL);

        protected:
        /**
         * Return the layout of the appender.
         * This method is the Layout accessor for subclasses of LayoutAppender.
         * @returns the Layout.
         **/
        Layout& _getLayout();

        private:
        Layout* _layout;
    };
}

#endif // CF_LOGCPP_LAYOUTAPPENDER_HH

