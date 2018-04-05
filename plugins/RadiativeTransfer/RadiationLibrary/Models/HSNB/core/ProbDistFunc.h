#ifndef COOLFluiD_RadiativeTransfer_PROB_DIST_FUNC_H
#define COOLFluiD_RadiativeTransfer_PROB_DIST_FUNC_H

#include <algorithm>
#include <cassert>
#include <vector>
#include <fstream>

#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/Random.h"

using namespace std;

class ProbDistFunc
{
public:

    /// Empty constructor.
    ProbDistFunc() { }

    /// Constructs a new probability distribution function.
    template <typename Vector>
    ProbDistFunc(const Vector& pdf) :
        m_cumulative_pdf(pdf.size())
    {
        // Make sure all the probabilities are positive
        for (int i = 0; i < pdf.size(); ++i)
            assert(pdf[i] >= 0.0);


        // Create a cumulative PDF
        m_cumulative_pdf[0] = pdf[0];
        for (int i = 1; i < m_cumulative_pdf.size(); ++i)
            m_cumulative_pdf[i] = pdf[i] + m_cumulative_pdf[i-1];

        // Normalize
        const double sum = m_cumulative_pdf.back();
        for (int i = 0; i < m_cumulative_pdf.size(); ++i)
            m_cumulative_pdf[i] /= sum;
    }

    /// Copy constructor.
    ProbDistFunc(const ProbDistFunc& pdf) :
        m_cumulative_pdf(pdf.m_cumulative_pdf)
    { }

    /// Assignment operator.
    ProbDistFunc& operator= (ProbDistFunc pdf) {
        std::swap(m_cumulative_pdf, pdf.m_cumulative_pdf);
        return *this;
    }

    /// Returns the size of the PDF.
    std::size_t size() const {
        return m_cumulative_pdf.size();
    }

    /// Returns a random index following the PDF.
    std::size_t index() const
    {
        return std::distance(
            m_cumulative_pdf.begin(),
            std::lower_bound(
                m_cumulative_pdf.begin(), m_cumulative_pdf.end(),
                Random::uniform()));

    }

    /// Samples the PDF n times and returns a histogram.
    std::vector<unsigned int> sample(int n) {
        std::vector<unsigned int> samples(size(), 0);
        for (int i = 0; i < n; ++i)
            samples[index()]++;
        return samples;
    }

    /// Returns the probability that index i is returned by index().
    double probability(std::size_t i) const
    {
        if (i > 0)
            return m_cumulative_pdf[i] - m_cumulative_pdf[i-1];
        return m_cumulative_pdf[0];
    }

    /// Saves the PDF to a file.
    void save(const std::string& file_name) const {
        // Open the file
        std::ofstream f(file_name.c_str());
        f << std::scientific;
        for (int i = 0; i < size(); ++i) {
            f << i << " " << probability(i) << " " << m_cumulative_pdf[i] << "\n";
        }

        f.close();
    }

private:

    std::vector<double> m_cumulative_pdf;
};

#endif // PROB_DIST_FUNC_H

