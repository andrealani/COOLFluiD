#ifndef COOLFluiD_RadiativeTransfer_LBL_SPECTRAL_GRID_H
#define COOLFluiD_RadiativeTransfer_LBL_SPECTRAL_GRID_H

#include <algorithm>
#include <string>
#include <cstdlib>
#include <vector>
#include <list>
#include <iostream>
#include <cmath>

enum AdaptType {
	SMART,
	DA_SILVA
};

class LblSpectralGrid
{
public:

    LblSpectralGrid() : m_size(0), mp_grid(NULL) {}

    /**
     * Reads a spectral grid from a file.
     */
    LblSpectralGrid(std::string directory);

    /**
     * Generates a spectral grid with the given minimum and maximum wavenumbers.
     */
    LblSpectralGrid(double min, double max);

    /**
     * Generates a uniform spectral grid with between min and max with n nodes.
     */
    LblSpectralGrid(double min, double max, int n);

    /**
     * Constructs a LBL spectral grid using an adaptive technique based on the
     * supplied line data and maximum and minimum allowed wavenumbers.
     */
    template <typename T>
    LblSpectralGrid(
        std::vector<T>& lines, double min, double max, AdaptType type = SMART);

    /**
     * Copy constructor.
     */
    LblSpectralGrid(const LblSpectralGrid& grid)
        : m_size(grid.m_size)
    {
        mp_grid = new double [m_size];
        std::copy(grid.mp_grid, grid.mp_grid+m_size, mp_grid);
    }

    /**
     * Destructor.
     */
    ~LblSpectralGrid();

    /**
     * Assignment operator.
     */
    LblSpectralGrid& operator= (LblSpectralGrid grid)
    {
        swap(*this, grid);
        return *this;
    }

    void save(const std::string& file);
    void load(const std::string& file);

    /**
     * Returns number of points in the grid.
     */
    int size() const { return m_size; }

    /**
     * Grid point access.
     */
    double operator[] (int i) const { return mp_grid[i]; }

    /**
     * Returns the maximum wavenumber of the grid in cm-1.
     */
    double max() const { return mp_grid[m_size-1]; }

    /**
     * Returns the minimum wavenumber of the grid in cm-1.
     */
    double min() const { return mp_grid[0]; }

    /**
     * Returns the index to the grid point with the given wavenumber.
     */
    int index(double sigma) const;

    /**
     * Returns the address of the wavenumber grid.
     */
    const double* const ptr() const { return mp_grid; }

    friend void swap(LblSpectralGrid&, LblSpectralGrid&);

private:

    template <typename T>
    void addPointsDaSilva(
        std::vector<T>& lines, double min, double max,
        std::list<double>& list);

    void fillPointsForward(
        std::list<double>& points, std::list<double>::iterator iter);

    void fillPointsReverse(
        std::list<double>& points, std::list<double>::iterator iter);

    template <typename T>
    void createSmartGrid(
        std::vector<T>& lines, double min, double max,
        std::list<double>& points);

private:

    size_t  m_size;
    double* mp_grid;
};

void swap(LblSpectralGrid& g1, LblSpectralGrid& g2);


struct ClosePredicate {
    ClosePredicate(double c)
        : close(c) { }
    bool operator() (double a, double b) {
        return (b-a < close);
    }
    double close;
};

template <typename T>
LblSpectralGrid::LblSpectralGrid(
    std::vector<T>& lines, double min, double max, AdaptType type)
{
    std::list<double> list;

    switch (type) {
    case SMART:
        createSmartGrid(lines, min, max, list);
        break;
    case DA_SILVA:
        addPointsDaSilva(lines, min, max, list);
        break;
    }


    std::list<double>::iterator iter = list.begin();
    double min_dist = 1.0e-6;
    while ((iter = std::adjacent_find(iter, list.end(), ClosePredicate(min_dist)))
            != list.end())
        iter = list.erase(++iter);

    m_size  = list.size();
    mp_grid = new double [m_size];
    //std::cout << m_size << std::endl;
    std::copy(list.begin(), list.end(), mp_grid);
}

template <typename T>
void LblSpectralGrid::addPointsDaSilva(
        std::vector<T>& lines, double min, double max, std::list<double>& list)
{
    for (int i = 0; i < lines.size(); ++i) {
        T& line = lines[i];
        const double V  =
            std::sqrt(line.gaml*line.gaml + line.gamd*line.gamd);
        const double W  = line.gaml/std::atan(1.0) + 1.8*line.gamd;
        const double FW = 1.8/std::atan(1.0)*line.gaml + 5.8*line.gamd;
        const double AV = std::sqrt(W*FW);

        // Create the points for this line
        std::list<double> points;
        points.push_back(line.sigc);

        //points.push_back(line.sigc+0.0625*V);
        points.push_back(line.sigc+0.125*V);
        points.push_back(line.sigc+0.25*V);
        //points.push_back(line.sigc+0.5*V);
        points.push_back(line.sigc+W);
        points.push_back(line.sigc+FW);
        points.push_back(line.sigc+12.5*V);

        //points.push_front(line.sigc-0.0625*V);
        points.push_front(line.sigc-0.125*V);
        points.push_front(line.sigc-0.25*V);
        //points.push_front(line.sigc-0.5*V);
        points.push_front(line.sigc-W);
        points.push_front(line.sigc-FW);
        points.push_front(line.sigc-12.5*V);

        // Merge with the full list
        list.merge(points);
    }

    while (list.front() < min)
        list.pop_front();

    while (list.back() > max)
        list.pop_back();
}

template <typename T>
struct LineDataSorter {
    bool operator() (const T& d1, const T& d2) {
        return d1.sigc < d2.sigc;
    }
};

template <typename T>
void LblSpectralGrid::createSmartGrid(
    std::vector<T>& lines, double min, double max,
    std::list<double>& points)
{
    // Start by sorting the line data by line center position
    std::sort(lines.begin(), lines.end(), LineDataSorter<T>());

    double v1  =
        std::sqrt(lines[0].gaml*lines[0].gaml + lines[0].gamd*lines[0].gamd);
    double w1  = lines[0].gaml/std::atan(1.0) + 1.8*lines[0].gamd;
    double f1 = 1.8/std::atan(1.0)*lines[0].gaml + 5.8*lines[0].gamd;
    double a1 = std::sqrt(w1*f1);
    double c1 = lines[0].sigc;

    double v2, w2, f2, a2, x, y, c2;

    // Left side of first line
    points.push_back(c1-f1);
    points.push_back(c1-a1);
    points.push_back(c1-w1);
    points.push_back(c1-0.5*v1);
    points.push_back(c1-0.125*v1);
    points.push_back(c1);

    for (int i = 0; i < lines.size()-1; ++i) {
        // Compute sizing parameters
        v2 = std::sqrt(
            lines[i+1].gaml*lines[i+1].gaml + lines[i+1].gamd*lines[i+1].gamd);
        w2 = lines[i+1].gaml/std::atan(1.0) + 1.8*lines[i+1].gamd;
        f2 = 1.8/std::atan(1.0)*lines[i+1].gaml + 5.8*lines[i+1].gamd;
        a2 = std::sqrt(w2*f2);
        c2 = lines[i+1].sigc;

        //cout << v2 << endl;

        // Distances to "center" from first and second line centers
        x  = (c2-c1)/(1.0+v2/v1);
        y  = (c2-c1)-x;

        // Line 1 points
        if (x > 0.125*v1) {
            points.push_back(c1+0.125*v1);
            if (x > 0.5*v1) {
                points.push_back(c1+0.5*v1);
                if (x > w1) {
                    points.push_back(c1+w1);
                    if (x > a1) {
                        points.push_back(c1+a1);
                        if (x > f1) {
                            points.push_back(c1+f1);
//                            if (x > 12.5*v1) {
//                                points.push_back(c1+12.5*v1);
//                                if (x > 100.0*v1) {
//                                    points.push_back(c1+100.0*v1);
                                    //points.push_back(c1+x);
                                } } } } } //} }

        // Line 2 points (note they must be inserted to maintain order)
        std::list<double>::iterator iter = points.insert(points.end(), c2);
        if (y > 0.125*v2) {
            iter = points.insert(iter, c2-0.125*v2);
            if (y > 0.5*v2) {
                iter = points.insert(iter, c2-0.5*v2);
                if (y > w2) {
                    iter = points.insert(iter, c2-w2);
                    if (y > a2) {
                        iter = points.insert(iter, c2-a2);
                        if (y > f2) {
                            iter = points.insert(iter, c2-f2);
//                            if (y > 12.5*v2) {
//                                iter = points.insert(iter, c2-12.5*v2);
//                                if (y > 100.0*v2) {
//                                    iter = points.insert(iter, c2-100.0*v2);
                                } } } } } //} }


        fillPointsForward(points, iter);
        fillPointsReverse(points, iter);

        // Copy the sizing parameters
        v1 = v2;
        w1 = w2;
        f1 = f2;
        a1 = a2;
        c1 = c2;
    }

    // Points on right side of last line
    points.push_back(c1+0.125*v1);
    points.push_back(c1+0.5*v1);
    points.push_back(c1+w1);
    points.push_back(c1+a1);
    points.push_back(c1+f1);

    // Remove any points below min or above max
    while (points.front() < min)
        points.pop_front();

    while (points.back() > max)
        points.pop_back();

    // Finally fill in points on either side between min and max
    points.push_front(min);
    fillPointsReverse(points, ++(points.begin()));
    points.push_back(max);
    fillPointsForward(points, --(points.end()));
}


#endif // LBL_SPECTRAL_GRID_H
