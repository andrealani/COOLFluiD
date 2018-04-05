#ifndef COOLFluiD_RadiativeTransfer_OPERATORS_H
#define COOLFluiD_RadiativeTransfer_OPERATORS_H

typedef struct {
    void operator () (double& x, const double& y) const { x += y; }
} SumEq;

typedef struct {
    void operator () (double& x, const double& y) const { x *= y; }
} TimesEq;

typedef struct {
    void operator () (double& x, const double& y) const { x = y; }
} Eq;

#endif
