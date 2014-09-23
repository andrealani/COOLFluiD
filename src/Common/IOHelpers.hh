// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef UTILS_IOHELPERS_HH
#define UTILS_IOHELPERS_HH

#include "Array2D.hh"

namespace std
{

/*template <typename T, template <class> class CONT>
std::ostream & operator << (std::ostream & O, const CONT<T> & A)
{
    O << '(';
    unsigned int Count = 0;
    const unsigned int Size = A.size();
    for (typename CONT<T>::const_iterator I = A.begin(); I!=A.end(); ++I)
    {
        ++Count;
        O << *I;
        if (Count!=Size)
            O << ',';
    }
    O << ')';
    return O;
}*/

template <typename T>
std::ostream & operator << (std::ostream & O, const COOLFluiD::Common::Array2D<T> & A)
{
    O << '[';
    for (unsigned int i=0; i<A.GetSize(); ++i)
    {
        O << '(';
        for (unsigned int j=0; j<A.GetElementSize (); ++j)
        {
            O << A.Get(i,j);
            if (j!=(A.GetElementSize()-1))
                O << ',';
        }
        O << ')';
    }
    O << ']';
    return O;
}

template <typename T>
std::ostream & operator << (std::ostream & O, const std::vector<T> & Array)
{
    O << '(';
    for (unsigned int i=0; i<Array.size(); ++i)
    {
        O << Array[i];
        if (i!=(Array.size()-1))
            O << ',';
    }
    O << ')';
    return O;
}

template <typename T, unsigned int N>
std::ostream & operator << (std::ostream & O, const T * Array)
{
    O << '(';
    for (unsigned int i=0; i<N; ++i)
    {
        O << Array[i];
        if (i!=(N-1))
            O << ',';
    }
    O << ')';
    return O;
}

}

#endif
