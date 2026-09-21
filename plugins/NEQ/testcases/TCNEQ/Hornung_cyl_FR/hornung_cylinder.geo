// Hornung cylinder, front half, Q2 quads, for the FR TCNEQ cases in this
// directory. Body radius 12.75 mm, shock adapted outer spline.
//
// Two blocks:
//  - an inner boundary layer block between the wall and a concentric circle at
//    R + DBL. Its first cell is H1 = 1e-6 R (12.75 nm, cell Reynolds number of
//    order one) and NBL = 44 cells grow by PROGBL = 1.2, so the block is
//    DBL = 0.19 mm thick and its last cell (32 um) matches the first cell of
//    the outer block.
//  - the outer block, NRAD x NCIRC = 56 x 41 cells growing by PROG = 1.05.
//
// Both circles are discretised uniformly in theta, so the inner block is purely
// radial and cannot shear: the shear between the wall circle and the arc length
// parametrised outer spline stays in the outer block.
//
// NCIRC is odd, so the stagnation line (theta = 0) cuts a cell in half instead
// of running along a cell edge. The two quarter arcs are split at the node just
// below theta = 0 (THS = -90 + NLOW*180/NCIRC deg) to keep the spacing uniform
// in theta; with an even NCIRC this falls back to a split at theta = 0.
//
// Build:
//   gmsh -2 -order 2 -format msh2 hornung_cylinder.geo -o hornung_cylinder_raw.msh
// followed by a projection of the Q2 mid nodes of the inner block onto their
// arcs: gmsh leaves them on the chords, which folds the 13 nm wall cells.

// resolution
NCIRC = DefineNumber[ 41, Name "Cells around the body" ];
NRAD  = DefineNumber[ 56, Name "Cells in the radial direction" ];
PROG  = DefineNumber[ 1.05, Name "Radial growth ratio, from the wall out" ];
NBL    = DefineNumber[ 44,  Name "Cells in the boundary layer block" ];
PROGBL = DefineNumber[ 1.2, Name "Growth ratio in the boundary layer block" ];

R  = 0.01275;   // body radius [m]
H1 = 1e-6*R;    // first cell height at the wall
DBL = H1*(PROGBL^NBL - 1)/(PROGBL - 1);   // boundary layer block thickness
RB  = R + DBL;

// split of the half circle into two arcs (gmsh arcs must be below 180 deg)
NLOW = Floor(NCIRC/2);           // cells from theta = -90 deg to the split
NUP  = NCIRC - NLOW;             // cells from the split to theta = +90 deg
THS  = -Pi/2 + NLOW*Pi/NCIRC;    // split angle, a mesh node

Point(1) = {0, 0, 0};            // centre, for the wall arcs
Point(2) = {0, -R, 0};
Point(3) = {R*Cos(THS), R*Sin(THS), 0};
Point(4) = {0,  R, 0};
Point(5) = {0, -RB, 0};          // boundary layer block outer circle
Point(6) = {RB*Cos(THS), RB*Sin(THS), 0};
Point(7) = {0,  RB, 0};

// outer boundary control points, from theta = -90 to +90 deg
Point(10) = { 1.69990624e-17, -3.36400899e-02, 0};   // theta =  -90.00
Point(11) = { 4.80432510e-03, -3.06043865e-02, 0};   // theta =  -81.08
Point(12) = { 8.92981601e-03, -2.74831477e-02, 0};   // theta =  -72.00
Point(13) = { 1.23541456e-02, -2.42973577e-02, 0};   // theta =  -63.05
Point(14) = { 1.53519435e-02, -2.11301375e-02, 0};   // theta =  -54.00
Point(15) = { 1.77873301e-02, -1.78052737e-02, 0};   // theta =  -45.03
Point(16) = { 1.98659309e-02, -1.44334436e-02, 0};   // theta =  -36.00
Point(17) = { 2.13882472e-02, -1.09054938e-02, 0};   // theta =  -27.02
Point(18) = { 2.25103071e-02, -7.31404213e-03, 0};   // theta =  -18.00
Point(19) = { 2.30692727e-02, -3.65640552e-03, 0};   // theta =   -9.01
Point(20) = { 2.33083834e-02,  5.17550079e-18, 0};   // theta =    0.00
Point(21) = { 2.30858474e-02,  3.65949736e-03, 0};   // theta =    9.01
Point(22) = { 2.24996339e-02,  7.31057420e-03, 0};   // theta =   18.00
Point(23) = { 2.13824757e-02,  1.09024670e-02, 0};   // theta =   27.02
Point(24) = { 1.98394403e-02,  1.44141971e-02, 0};   // theta =   36.00
Point(25) = { 1.77981508e-02,  1.78179340e-02, 0};   // theta =   45.03
Point(26) = { 1.53503429e-02,  2.11279345e-02, 0};   // theta =   54.00
Point(27) = { 1.23703898e-02,  2.43320222e-02, 0};   // theta =   63.05
Point(28) = { 8.93315013e-03,  2.74934091e-02, 0};   // theta =   72.00
Point(29) = { 4.81051402e-03,  3.06467479e-02, 0};   // theta =   81.08
Point(30) = { 2.05986142e-18,  3.36400899e-02, 0};   // theta =   90.00

Circle(1) = {2, 1, 3};           // wall, from -90 deg to the split
Circle(2) = {3, 1, 4};           // wall, from the split to +90 deg
Spline(3) = {10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};   // outer boundary
Circle(6) = {5, 1, 6};           // BL block outer circle, from -90 deg to the split
Circle(7) = {6, 1, 7};           // BL block outer circle, from the split to +90 deg
Line(4) = {10, 5};        // outflow, lower, outer block (x = 0)
Line(5) = {7, 30};        // outflow, upper, outer block (x = 0)
Line(8) = {5, 2};         // outflow, lower, BL block
Line(9) = {4, 7};         // outflow, upper, BL block

// counter-clockwise, so the quads come out with a positive Jacobian. Walking
// the other way round gives every cell a negative determinant, which COOLFluiD
// flags cell by cell at setup.
Line Loop(1) = {3, -5, -7, -6, -4};
Plane Surface(1) = {1};
Line Loop(2) = {6, 7, -9, -2, -1, -8};
Plane Surface(2) = {2};

// structured quads: NCIRC around the body, NRAD outwards, clustered at the wall
Transfinite Line{1, 6} = NLOW + 1;
Transfinite Line{2, 7} = NUP + 1;
Transfinite Line{3}    = NCIRC + 1;
Transfinite Line{8} = NBL + 1 Using Progression 1/PROGBL;
Transfinite Line{9} = NBL + 1 Using Progression PROGBL;
Transfinite Line{4} = NRAD + 1 Using Progression 1/PROG;
Transfinite Line{5} = NRAD + 1 Using Progression PROG;
Transfinite Surface{1} = {5, 7, 30, 10};
Transfinite Surface{2} = {2, 4, 7, 5};
Recombine Surface{1};
Recombine Surface{2};

// complete (9 node) quads, and put the new mid nodes on the real geometry
Mesh.SecondOrderIncomplete = 0;
Mesh.SecondOrderLinear     = 0;
Mesh.MshFileVersion        = 2.2;

Physical Line("Wall")   = {1, 2};
Physical Line("Inlet")  = {3};
Physical Line("Outlet") = {4, 5, 8, 9};
Physical Surface("InField") = {1, 2};
