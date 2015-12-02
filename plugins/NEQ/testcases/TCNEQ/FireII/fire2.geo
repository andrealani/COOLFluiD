r = .805;
t = .05;
rn = .036;
angle = .372988;

bl1 = 0.0067372;

bl2 = .03;

sh1 = .038;
sh2 = .044;

dr    = .005;
dr_bl = .05  * dr;
dr_sh = .1   * dr ;


Point(1) = {-0.  , 0. , 0. , dr_bl};
Point(2) = {-bl1 , 0. , 0. , dr   };
Point(3) = {-bl2 , 0. , 0. , dr   };
Point(4) = {-sh1 , 0. , 0. , dr_sh};
Point(5) = {-sh2 , 0. , 0. , dr_sh};
Point(6) = {-t   , 0. , 0. , dr   };


Point(7 ) = { -(r + 0  ) * Cos(angle) + r, (r + 0  ) * Sin(angle) ,0., dr };
Point(8 ) = { -(r + bl1) * Cos(angle) + r, (r + bl1) * Sin(angle) ,0., dr };
Point(9 ) = { -(r + bl2) * Cos(angle) + r, (r + bl2) * Sin(angle) ,0., dr };
Point(10) = { -(r + sh1) * Cos(angle) + r, (r + sh1) * Sin(angle) ,0., dr };
Point(11) = { -(r + sh2) * Cos(angle) + r, (r + sh2) * Sin(angle) ,0., dr };
Point(12) = { -(r + t  ) * Cos(angle) + r, (r + t  ) * Sin(angle) ,0., dr };
                
Point(13) = {r, 0. ,0., dr};

//first section

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};

Circle(6 ) = {1, 13,  7 };
Circle(7 ) = {2, 13,  8 };
Circle(8 ) = {3, 13,  9 };
Circle(9 ) = {4, 13, 10 };
Circle(10) = {5, 13, 11 };
Circle(11) = {6, 13, 12 };

Line(12) = {7 ,  8};
Line(13) = {8 ,  9};
Line(14) = {9 , 10};
Line(15) = {10, 11};
Line(16) = {11, 12};

Line Loop(1) = {-1, -7,  12 , 6};
Line Loop(2) = {-2, -8,  13 , 7};
Line Loop(3) = {-3, -9,  14 , 8};
Line Loop(4) = {-4, -10, 15 , 9};
Line Loop(5) = {-5, -11, 16 , 10};

Plane Surface(1) = {1}; 
Plane Surface(2) = {2}; 
Plane Surface(3) = {3}; 
Plane Surface(4) = {4}; 
Plane Surface(5) = {5}; 

//second section

x00 = -(r-rn) * Cos(angle) + r;
y00 = (r-rn)*Sin(angle);

Point(14) = { x00, y00          ,0. ,dr  };
Point(15) = { x00, y00 +rn      ,0. ,dr  };
Point(16) = { x00, y00 +rn +bl1 ,0. ,dr  };

xpos =  (r - rn ) * Cos(angle) - r ;  
angle1 = Acos(  ( r + xpos ) / ( r + bl2 )  );
angle2 = Acos(  ( r + xpos ) / ( r + sh1 )  );
angle3 = Acos(  ( r + xpos ) / ( r + sh2 )  );
angle4 = Acos(  ( r + xpos ) / ( r + t   )  );

//Point(17) = {  -xpos , (r+bl2) * Sin(angle1), 0. ,dr  };
//Point(18) = {  -xpos , (r+sh1) * Sin(angle2), 0. ,dr  };
//Point(19) = {  -xpos , (r+sh2) * Sin(angle3), 0. ,dr  };

sh22 = .414;
sh12 = sh22-(sh2-sh1)*1.5;
bl22 = sh22-(sh2-bl2)*1.6; 

Point(19) = {  -xpos , sh22, 0. ,dr  };

Point(18) = {  -xpos , sh12, 0. ,dr  };
Point(17) = {  -xpos , bl22, 0. ,dr  };

Point(20) = {  -xpos , (r+t  ) * Sin(angle4), 0. ,dr  };


//Add a B-spline between the inner points and its known derivatives
//In the case of a B-spline that means adding a control point that 
//is the intersection of the lines defined by the points and its 
// derivatives

//POINT3
//Point(11) = { -(r + sh2) * Cos(angle) + r, (r + sh2) * Sin(angle) ,0., dr };
x01= -(r + sh2) * Cos(angle) + r;
y01= (r + sh2) * Sin(angle);
d11= Sin(angle);
d12= Cos(angle);

//Point(19) = {  -xpos , sh22, 0. ,dr  };
x02= -xpos;
y02= sh22;
d21= Cos(-132.15*Pi/180.);
d22= Sin(-132.14*Pi/180.);

t1 = (d21*(y01 - y02))/(d11*d22 - d12*d21) - (d22*(x01 - x02))/(d11*d22 - d12*d21);
t2 = (d11*(y01 - y02))/(d11*d22 - d12*d21) - (d12*(x01 - x02))/(d11*d22 - d12*d21);
x_intersection = x01 + d11*t1;
y_intersection = y01 + d12*t1;

Point(103) = {x_intersection, y_intersection, 0, dr};

//POINT2
//Point(10) = { -(r + sh1) * Cos(angle) + r, (r + sh1) * Sin(angle) ,0., dr };
x01= -(r + sh1) * Cos(angle) + r ;
y01= (r + sh1) * Sin(angle) ;
d11= Sin(angle);
d12= Cos(angle);

//Point(18) = {  -xpos , sh12, 0. ,dr  };
x02= -xpos;
y02= sh12;
d21= d21;
d22= d22;

t1 = (d21*(y01 - y02))/(d11*d22 - d12*d21) - (d22*(x01 - x02))/(d11*d22 - d12*d21);
t2 = (d11*(y01 - y02))/(d11*d22 - d12*d21) - (d12*(x01 - x02))/(d11*d22 - d12*d21);
x_intersection = x01 + d11*t1;
y_intersection = y01 + d12*t1;

Point(102) = {x_intersection, y_intersection, 0, dr};

//POINT3
//Point(9 ) = { -(r + bl2) * Cos(angle) + r, (r + bl2) * Sin(angle) ,0., dr };
x01= -(r + bl2) * Cos(angle) + r;
y01= (r + bl2) * Sin(angle);
d11= Sin(angle);
d12= Cos(angle);

//Point(17) = {  -xpos , bl22, 0. ,dr  };
x02= -xpos;
y02= bl22;
d21= d21;
d22= d22;

t1 = (d21*(y01 - y02))/(d11*d22 - d12*d21) - (d22*(x01 - x02))/(d11*d22 - d12*d21);
t2 = (d11*(y01 - y02))/(d11*d22 - d12*d21) - (d12*(x01 - x02))/(d11*d22 - d12*d21);
x_intersection = x01 + d11*t1;
y_intersection = y01 + d12*t1;

Point(101) = {x_intersection, y_intersection, 0, dr};

Circle(17)  = {7,  14,  15 };
Circle(18)  = {8,  14,  16 };
BSpline(19) = {9,  101, 17 };
BSpline(20) = {10, 102, 18 };
BSpline(21) = {11, 103, 19 };
Circle(22)  = {12, 13,  20 };

Line(23) = {15, 16};
Line(24) = {16, 17};
Line(25) = {17, 18};
Line(26) = {18, 19};
Line(27) = {19, 20};

Line Loop(6 ) = {-12, -18, 23, 17};
Line Loop(7 ) = {-13, -19, 24, 18};
Line Loop(8 ) = {-14, -20, 25, 19};
Line Loop(9 ) = {-15, -21, 26, 20};
Line Loop(10) = {-16, -22, 27, 21};

Plane Surface(6 ) = {6 } ;
Plane Surface(7 ) = {7 };
Plane Surface(8 ) = {8 };
Plane Surface(9 ) = {9 };
Plane Surface(10) = {10};

//Transfinite Definitions

// [wall -> boundary layer] region
// First cell needs to be a= 1.23646059e-6
// obtained solving  "a * (1-r^m) /(1-r) = bl1" for r
m=30;

Transfinite Line{1, 12, 23} = m Using Progression 1.1894799;
//Transfinite Line{1, 12, 23} = m Using Progression 1.2;

// [boundary layer -> inner layer] region
Transfinite Line{2, 13, 24} = 18;

//inner layer -> shock layer
Transfinite Line{-3, -14, -25} = 10 Using Progression 1.2;

//shock region
Transfinite Line{4, 15, 26} = 18;

// [shock -> inlet] region
Transfinite Line{5, 16, 27} = 8 Using Progression 1.25;


Transfinite Line{6, 7, 8, 9, 10,11} = 40;
Transfinite Line{17,18,19,20,21,22} = 15;


Transfinite Surface "*";
Recombine Surface "*";

Smoother Surface{7} = 100;

Physical Line("Wall") = {17, 6};
Physical Line("Outlet") = {23, 24, 25, 26, 27};
Physical Line("Inlet") = {22, 11};
Physical Line("Symmetry") = {1, 2, 3, 4, 5};
Physical Surface("InnerFaces") = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
