#include <iostream>

#include "LineOfSight.h"

using namespace std;

LineOfSight::LineOfSight(const FieldData& field, 
                         const int& index, const double& mu)
    : m_field(field)
{

    switch (m_field.geom()) {
        case CARTESIAN: 
            buildCartesianPath(index, mu);
            break;
        case SPHERICAL: 
            buildSphericalPath(index, mu);
            break;
        default: 
            cout << "ERROR " << endl;
    }


}

LineOfSight::LineOfSight(const FieldData& field, int i1, int i2)
    : m_field(field)
{
    for (int i=i2-1; i >= i1; i--) {
       double length = (m_field.loc(i+1)-m_field.loc(i));
       m_path.push_back(make_pair(length, i));
   }
}

void LineOfSight::buildCartesianPath(const int& index, const double& mu)
{

    double length;

    if (mu <=  1. & mu > 0.) {
       m_wall = 0;
       for (int i=index-1; i>=0; i--) {
           length = (m_field.loc(i+1)-m_field.loc(i))/mu;
           m_path.push_back(make_pair(length, i));
       }
    } else if (mu >= -1. & mu < 0.) {
       m_wall = 1;
       for (int i=index; i<m_field.nPoints()-1; i++) {
           length = (m_field.loc(i)-m_field.loc(i+1))/mu;
           m_path.push_back(make_pair(length, i));
       }
    } else {
        cout << "CANNOT BUILD THE PATH !" << endl;
    }

}

void LineOfSight::buildSphericalPath(const int& index, const double& mu)
{
//    cout << "BUILD SPHERICAL PATH " << index << " " << mu << endl;

    const double R0 = m_field.loc(0);
    const double r = m_field.loc(index);
    const double y = r * sqrt(1-mu*mu);
//    cout << "R0 = " << R0
//         << " R = " << r  
//         << " Y = " << y << endl;

    double length, xi, xj, ri, rj;
    int imin;

    if (mu <=  1. & mu > 0.) {
        ri = r; 
        xi = sqrt(ri*ri - y*y);

        if (y <= R0) {
            m_wall = 0;

            for (int i=index-1; i>=0; i--) {
                rj = m_field.loc(i); 
                xj = sqrt(rj*rj - y*y);
                length = xi - xj;
                m_path.push_back(make_pair(length, i));
//                cout << "I = " << i << ", LEN = " << length << endl;
                ri = rj;
                xi = xj;
            }
        } else {
            m_wall = 1;

            imin = index -1;
            while (m_field.loc(imin) > y) imin--;
 
            for (int i=index-1; i>imin; i--) {
                rj = m_field.loc(i); 
                xj = sqrt(rj*rj - y*y);
                length = xi - xj;
                m_path.push_back(make_pair(length, i));
//                cout << "I = " << i << ", LEN = " << length << endl;
                ri = rj;
                xi = xj;
            }

                length = 2*xi;
                m_path.push_back(make_pair(length, imin));
//                cout << "I = " << imin << ", LEN = " << length << endl;

            for (int i=imin+1; i<m_field.nPoints()-1; i++) {
                rj = m_field.loc(i+1); 
                xj = sqrt(rj*rj - y*y);
                length = xj - xi;
                m_path.push_back(make_pair(length, i));
//                cout << "I = " << i << ", LEN = " << length << endl;
                ri = rj;
                xi = xj;
            }

        }

    } else if (mu >= -1. & mu < 0.) {
        ri = r; 
        xi = sqrt(ri*ri - y*y);
        m_wall = 1;
        for (int i=index; i<m_field.nPoints()-1; i++) {
            rj = m_field.loc(i+1); 
            xj = sqrt(rj*rj - y*y);
            length = xj - xi;
            m_path.push_back(make_pair(length, i));
//            cout << "I = " << i << ", LEN = " << length << endl;
            ri = rj;
            xi = xj;
        }

    } else if (mu == 0.0) {
        ri = r;
        xi =0.0;
        m_wall = 1;
        for (int i=index; i<m_field.nPoints()-1; i++) {
            rj = m_field.loc(i+1); 
            xj = sqrt(rj*rj - ri*ri);
            length = xj - xi;
            m_path.push_back(make_pair(length, i));
//            cout << "I = " << i << ", LEN = " << length << endl;
            ri = rj;
            xi = xj;
        }

    } else {
        cout << "CANNOT BUILD THE PATH !" << endl;
    }

}
