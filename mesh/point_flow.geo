Point(1) = {0, 0, 0, 1.0};
Point(2) = {0, 1, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {1, 0, 0, 1.0};

Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Curve {1, 3} = 11 Using Progression 1;
Transfinite Curve {2, 4} = 11 Using Progression 1;

Transfinite Surface {1};
Recombine Surface {1};

Extrude {0, 0, 1} {
  Surface{1};
  Layers {1};
  Recombine;
};

Physical Volume("interior", 10) = {1};
Physical Surface("wall", 13) = {25, 21, 17, 13};
Physical Surface("sym", 11) = {26, 1};
