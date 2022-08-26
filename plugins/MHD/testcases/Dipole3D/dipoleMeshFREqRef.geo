

Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0;
Mesh.MshFileVersion = 2.2;

cl2 = .0490873852;

r = 1.0;// inner radius
R = 20.795;// outer radius
t = 17;//17;//8.0;//4.0;//4.0;// distance in xy-plane from origin to corner points of polar square patch on outer sphere

nbSq = 11;//5;//6;// number of points on one side of square patch
nbNS = 13;//11;//21;// number of points on north-south line
nbN = 26;//15;//26;// number of points normal to inner sphere
prog = 1.3;//1.9;// point progression factor normal to inner sphere

NSEqBumpOuter = 4.0;
NSEqBumpInner = 4.0;

xyR = Sqrt(2.0)/2.0*t;
zR = Sqrt(R^2-t^2);

xyr = r/R*xyR;
zr = r/R*zR;

Point(1) = {0, 0, 0, cl2};
Point(2) = {xyr, xyr, zr, cl2};
Point(3) = {-xyr, xyr, zr, cl2};
Point(4) = {-xyr, -xyr, zr, cl2};
Point(5) = {xyr, -xyr, zr, cl2};
Point(6) = {xyr, xyr, -zr, cl2};
Point(7) = {-xyr, xyr, -zr, cl2};
Point(8) = {-xyr, -xyr, -zr, cl2};
Point(9) = {xyr, -xyr, -zr, cl2};
Point(10) = {xyR, xyR, zR, cl2};
Point(11) = {-xyR, xyR, zR, cl2};
Point(12) = {-xyR, -xyR, zR, cl2};
Point(13) = {xyR, -xyR, zR, cl2};
Point(14) = {xyR, xyR, -zR, cl2};
Point(15) = {-xyR, xyR, -zR, cl2};
Point(16) = {-xyR, -xyR, -zR, cl2};
Point(17) = {xyR, -xyR, -zR, cl2};
//+
Circle(1) = {4, 1, 5};
//+
Circle(2) = {5, 1, 2};
//+
Circle(3) = {2, 1, 3};
//+
Circle(4) = {3, 1, 4};
//+
Circle(5) = {7, 1, 8};
//+
Circle(6) = {8, 1, 9};
//+
Circle(7) = {9, 1, 6};
//+
Circle(8) = {6, 1, 7};
//+
Circle(9) = {5, 1, 9};
//+
Circle(10) = {4, 1, 8};
//+
Circle(11) = {3, 1, 7};
//+
Circle(12) = {2, 1, 6};
//+
Circle(13) = {12, 1, 11};
//+
Circle(14) = {11, 1, 10};
//+
Circle(15) = {10, 1, 13};
//+
Circle(16) = {13, 1, 12};
//+
Circle(17) = {15, 1, 16};
//+
Circle(18) = {16, 1, 17};
//+
Circle(19) = {17, 1, 14};
//+
Circle(20) = {14, 1, 15};
//+
Circle(21) = {12, 1, 16};
//+
Circle(22) = {11, 1, 15};
//+
Circle(23) = {10, 1, 14};
//+
Circle(24) = {13, 1, 17};
//+
Line(25) = {4, 12};
//+
Line(26) = {3, 11};
//+
Line(27) = {2, 10};
//+
Line(28) = {5, 13};
//+
Line(29) = {7, 15};
//+
Line(30) = {8, 16};
//+
Line(31) = {9, 17};
//+
Line(32) = {6, 14};
//+
Curve Loop(1) = {18, -31, -6, 30};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {31, 19, -32, -7};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {32, 20, -29, -8};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {30, -17, -29, 5};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {25, -16, -28, -1};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {28, -15, -27, -2};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {27, -14, -26, -3};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {25, 13, -26, 4};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {1, 2, 3, 4};
//+
Surface(9) = {9} In Sphere {1};
//+
Curve Loop(10) = {10, -5, -11, 4};
//+
Surface(10) = {10} In Sphere {1};
//+
Curve Loop(11) = {10, 6, -9, -1};
//+
Surface(11) = {11} In Sphere {1};
//+
Curve Loop(12) = {9, 7, -12, -2};
//+
Surface(12) = {12} In Sphere {1};
//+
Curve Loop(13) = {12, 8, -11, -3};
//+
Surface(13) = {13} In Sphere {1};
//+
Curve Loop(14) = {5, 6, 7, 8};
//+
Surface(14) = {14} In Sphere {1};
//+
Curve Loop(15) = {13, 14, 15, 16};
//+
Surface(15) = {15} In Sphere {1};
//+
Curve Loop(16) = {14, 23, 20, -22};
//+
Surface(16) = {16} In Sphere {1};
//+
Curve Loop(17) = {15, 24, 19, -23};
//+
Surface(17) = {17} In Sphere {1};
//+
Curve Loop(18) = {16, 21, 18, -24};
//+
Surface(18) = {18} In Sphere {1};
//+
Curve Loop(19) = {13, 22, 17, -21};
//+
Surface(19) = {19} In Sphere {1};
//+
Curve Loop(20) = {18, 19, 20, 17};
//+
Surface(20) = {20} In Sphere {1};
//+
Transfinite Curve {13, 14, 15, 16, 20, 19, 18, 17, 3, 4, 1, 2, 5, 6, 7, 8} = nbSq Using Progression 1.0;
//+
Transfinite Curve {10, 11, 12, 9} = nbNS Using Bump NSEqBumpInner;
//+
Transfinite Curve {23, 24, 21, 22} = nbNS Using Bump NSEqBumpOuter;
//+
Transfinite Curve {25, 26, 27, 28, 30, 29, 32, 31} = nbN Using Progression prog;
//+
Transfinite Surface {15};
//+
Transfinite Surface {16};
//+
Transfinite Surface {17};
//+
Transfinite Surface {18};
//+
Transfinite Surface {19};
//+
Transfinite Surface {20};
//+
Transfinite Surface {5};
//+
Transfinite Surface {7};
//+
Transfinite Surface {6};
//+
Transfinite Surface {8};
//+
Transfinite Surface {1};
//+
Transfinite Surface {3};
//+
Transfinite Surface {2};
//+
Transfinite Surface {4};
//+
Transfinite Surface {9};
//+
Transfinite Surface {14};
//+
Transfinite Surface {13};
//+
Transfinite Surface {12};
//+
Transfinite Surface {11};
//+
Transfinite Surface {10};
//+
Recombine Surface {9, 14, 13, 12, 11, 10, 5, 7, 8, 6, 1, 2, 3, 4, 15, 16, 17, 18, 19, 20};
//+
Curve Loop(21) = {21, -30, -10, 25};
//+
Plane Surface(21) = {21};
//+
Curve Loop(22) = {28, 24, -31, -9};
//+
Plane Surface(22) = {22};
//+
Curve Loop(23) = {23, -32, -12, 27};
//+
Plane Surface(23) = {23};
//+
Curve Loop(24) = {22, -29, -11, 26};
//+
Plane Surface(24) = {24};
//+
Transfinite Surface {24};
//+
Transfinite Surface {21};
//+
Transfinite Surface {22};
//+
Transfinite Surface {23};
//+
Recombine Surface {24, 21, 22, 23};
//+
Surface Loop(1) = {12, 22, 23, 6, 2, 17};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {13, 7, 3, 23, 24, 16};
//+
Volume(2) = {2};
//+
Surface Loop(3) = {10, 8, 4, 21, 24, 19};
//+
Volume(3) = {3};
//+
Surface Loop(4) = {11, 22, 21, 5, 1, 18};
//+
Volume(4) = {4};
//+
Transfinite Volume{3} = {7, 3, 4, 8, 15, 11, 12, 16};
//+
Transfinite Volume{4} = {8, 4, 5, 9, 16, 12, 13, 17};
//+
Transfinite Volume{1} = {9, 5, 2, 6, 17, 13, 10, 14};
//+
Transfinite Volume{2} = {6, 2, 3, 7, 14, 10, 11, 15};
//+
Surface Loop(5) = {15, 6, 5, 7, 8, 9};
//+
Volume(5) = {5};
//+
Surface Loop(6) = {14, 1, 4, 3, 2, 20};
//+
Volume(6) = {6};
//+
Transfinite Volume{5} = {11, 10, 13, 12, 3, 2, 5, 4};
//+
Transfinite Volume{6} = {7, 6, 9, 8, 15, 14, 17, 16};
//+
Physical Surface("Outlet") = {16, 19, 18, 17, 15, 20};
//+
Physical Surface("Inlet") = {11, 10, 13, 12, 9, 14};
//+
Physical Volume("InField") = {5, 3, 4, 1, 2, 6};
