#ifndef PHOTON_PATH_H
#define PHOTON_PATH_H

#include <cassert>

/// Represents a path that a photon takes.
class PhotonPath
{
public:

    /// Helper struct for storing cell data along a photon path
    struct CellData {
        CellData(int i, double d) : id(i), distance(d) { }
        int id;
        double distance;
    };

    /// Helper struct for storing wall data along a photon path
    struct WallData {
        WallData(int p, int i, double e) :
            pos(p), id(i), emissivity(e) { }
        int pos;
        int id;
        double emissivity;
    };

    PhotonPath() :
        m_energy(0), m_spec_loc(0), m_mechanism(-1)
    { }

    int nCells() const { return m_cells.size(); }
    int nWalls() const { return m_walls.size(); }

    void setEnergy(double energy) {
        assert(energy > 0);
        m_energy = energy;
    }
    double energy() const { return m_energy; }

    void setSpectralLocation(double loc) {
        assert(loc > 0);
        m_spec_loc = loc;
    }
    double spectralLocation() const { return m_spec_loc; }

    void setMechanism(int mechanism) {
        assert(mechanism >= 0);
        m_mechanism = mechanism;
    }
    int mechanism() const { return m_mechanism; }

    void addCell(int id, double distance) {
        assert(id >= 0);  assert(distance > 0);
        m_cells.push_back(CellData(id, distance));
    }

    int cellID(int ic) const {
        assert(ic >= 0); assert(ic < nCells());
        return m_cells[ic].id;
    }

    double cellDistance(int ic) const {
        assert(ic >= 0); assert(ic < nCells());
        return m_cells[ic].distance;
    }

    void addWall(int id, double emissivity) {
        assert(id >= 0); assert(emissivity >= 0);
        m_walls.push_back(WallData(nCells(), id, emissivity));
    }

    int wallPosition(int iw) const {
        assert(iw >= 0); assert(iw < nWalls());
        return m_walls[iw].pos;
    }

    int wallID(int iw) const {
        assert(iw >= 0); assert(iw < nWalls());
        return m_walls[iw].id;
    }

    double wallEmissivity(int iw) const {
        assert(iw >= 0); assert(iw < nWalls());
        return m_walls[iw].emissivity;
    }

private:

    double m_energy;
    double m_spec_loc;
    int m_mechanism;

    std::vector<CellData> m_cells;
    std::vector<WallData> m_walls;
};

#endif // PHOTON_PATH_H
