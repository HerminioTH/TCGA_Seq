// Gmsh project created on Thu May 13 14:28:40 2021
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(5) = {0, 1, 0, 1.0};
Point(6) = {0, 1, 1.0, 1.0};
Point(7) = {1, 1, 1.0, 1.0};
Point(8) = {1, 0, 1.0, 1.0};
Point(9) = {0, 0, 1.0, 1.0};

Line(1) = {5, 3};
Line(2) = {3, 2};
Line(3) = {2, 1};
Line(4) = {1, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 3};
Line(8) = {7, 8};
Line(9) = {8, 2};
Line(10) = {8, 9};
Line(11) = {9, 1};
Line(12) = {9, 6};

Line Loop(1) = {-2, 9, 8, -7};
Plane Surface(1) = {1};
Line Loop(2) = {-3, 11, 10, -9};
Plane Surface(2) = {2};
Line Loop(3) = {-4, -5, 12, -11};
Plane Surface(3) = {3};
Line Loop(4) = {-1, 7, 6, 5};
Plane Surface(4) = {4};
Line Loop(5) = {3, 4, 1, 2};
Plane Surface(5) = {5};
Line Loop(6) = {-10, -12, -6, -8};
Plane Surface(6) = {6};
Surface Loop(1) = {2, 5, 3, 4, 1, 6};
Volume(1) = {1};

Physical Volume("Body") = {1};
Physical Surface("East") = {1};
Physical Surface("West") = {3};
Physical Surface("South") = {2};
Physical Surface("North") = {4};
Physical Surface("Top") = {6};
Physical Surface("Bottom") = {5};


Transfinite Line {3, 2, 1, 4, 10, 8, 6, 12} = 5 Using Progression 1;
Transfinite Line {9, 7, 5, 11} = 10 Using Progression 1;

Transfinite Surface {2};
Transfinite Surface {1};
Transfinite Surface {4};
Transfinite Surface {3};
Transfinite Surface {5};
Transfinite Surface {6};
//+
Transfinite Volume{1} = {1, 2, 3, 5, 9, 8, 7, 6};
//+
Recombine Surface {2, 1, 4, 3, 5, 6};
