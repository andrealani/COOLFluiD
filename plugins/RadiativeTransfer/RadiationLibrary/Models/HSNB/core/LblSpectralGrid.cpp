#include <cassert>
#include <cmath>
#include <list>
#include <fstream>
#include "LblSpectralGrid.h"
#include "LblAtomicSystem.h"

using namespace std;

void swap(LblSpectralGrid& g1, LblSpectralGrid& g2)
{
    std::swap(g1.m_size,  g2.m_size);
    std::swap(g1.mp_grid, g2.mp_grid);
}

////==============================================================================
//
LblSpectralGrid::LblSpectralGrid(double min, double max)
{
    const double a = (0.12-0.01)/max;
    const double b = 0.01;

    m_size  = int(std::log((max+b/a)/(min+b/a))/std::log(a+1.0))+1;
    mp_grid = new double [m_size];

    for (size_t i = 0; i < m_size; ++i) {
        mp_grid[i] = std::pow(a+1.0,(double)i);
        mp_grid[i] = mp_grid[i]*min+(b/a)*(mp_grid[i]-1.0);
    }
}

//==============================================================================

LblSpectralGrid::LblSpectralGrid(double min, double max, int n)
{
    // Check data
    assert(n > 1);
    assert(min < max);

    m_size = n;
    mp_grid = new double [m_size];

    const double delta = (max-min)/double(m_size-1);
    for (int i = 0; i < m_size; ++i)
        mp_grid[i] = min + i*delta;
}

//==============================================================================

LblSpectralGrid::LblSpectralGrid(std::string directory)
{

    ifstream gf( (directory + "/spectral-grid.dat").c_str(), ios::binary);

    gf.read( (char*) &m_size, sizeof(size_t));
    mp_grid = new double [m_size];

    gf.read( (char*) mp_grid, m_size*sizeof(double));
    gf.close();

}

//==============================================================================

LblSpectralGrid::~LblSpectralGrid()
{
    if (mp_grid != 0) delete [] mp_grid;
}

void LblSpectralGrid::save(const std::string& file)
{
    ofstream f(file.c_str());

    f << m_size << "\n";
    f.setf(ios::scientific);
    f.precision(10);

    for (int i = 0; i < m_size; ++i)
        f << mp_grid[i] << "\n";
    f.close();
}

void LblSpectralGrid::load(const std::string& file)
{
    ifstream f(file.c_str());

    f >> m_size;

    if (mp_grid != NULL) delete [] mp_grid;
    mp_grid = new double [m_size];

    for (int i = 0; i < m_size; ++i)
        f >> mp_grid[i];
    f.close();
}

//==============================================================================
//
//size_t LblSpectralGrid::index(double sigma) const
//{
//    const double a = (0.12-0.01)/m_sigma_max;
//    const double b = 0.01;
//
//    if (sigma <= m_sigma_min)
//        return 0;
//    else if (sigma > m_sigma_max)
//        return m_size-1;
//    else
//        return (size_t) (std::log((sigma+b/a)/(m_sigma_min+b/a))/
//                std::log(a+1.0));
//}
//
////==============================================================================



int LblSpectralGrid::index(double sigma) const
{
    return (std::lower_bound(mp_grid, mp_grid+m_size, sigma) - mp_grid);
}


void LblSpectralGrid::fillPointsForward(
    std::list<double>& points, std::list<double>::iterator iter)
{
    const double ratio = 1.0;

    if (iter   == points.begin()) return;
    if (--iter == points.begin()) return;
    if (--iter == points.begin()) return;

    double v1 = *iter++;
    double v2 = *iter++;
    double v3 = *iter;

//    cout << "forward" << endl;
//    cout << "v1: " << v1 << endl;
//    cout << "v2: " << v2 << endl;
//    cout << "v3: " << v3 << endl;

    double goal = 2.0*ratio*(v2-v1) + v2;

    int i = 0;
    while (v3 > goal+1.0e-6){// && i++ < 5) {
        v3 = 0.5*(v2+v3);
        iter = points.insert(iter, v3);
    }

}

void LblSpectralGrid::fillPointsReverse(
    std::list<double>& points, std::list<double>::iterator iter)
{
    const double ratio = 1.0;

    if (iter == points.begin()) return;
    --iter;
    double v1 = *iter; iter++;
    double v2 = *iter; iter++;
    if (iter == points.end()) return;
    double v3 = *iter; iter--;

//    cout << "reverse" << endl;
//    cout << "v1: " << v1 << endl;
//    cout << "v2: " << v2 << endl;
//    cout << "v3: " << v3 << endl;

    double goal = v2 - 2.0*ratio*(v3-v2);

    while (v1 < goal-1.0e-6) {
        v1 = 0.5*(v1+v2);
        points.insert(iter, v1);
    }
}


