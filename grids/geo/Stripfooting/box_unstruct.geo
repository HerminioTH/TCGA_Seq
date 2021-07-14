gridSize = 0.044;
L = 1.0;
a = 0.25*L;
nz = 10;
nh = 5;
nt = 10;


Point(1) = {0, 0, 0, gridSize};
Point(2) = {0, L, 0, gridSize};
Point(3) = {L, 0, 0, gridSize};
Point(4) = {L, L, 0, gridSize};
Point(5) = {L, L, L, gridSize};
Point(6) = {L, 0, L, gridSize};
Point(7) = {0, 0, L, gridSize};
Point(8) = {0, L, L, gridSize};
Point(9) = {0, a, L, gridSize};
Point(10) = {a, a, L, gridSize};
Point(11) = {a, 0, L, gridSize};

Line(1) = {1, 3};
Line(2) = {3, 4};
Line(3) = {4, 2};
Line(4) = {2, 1};
Line(5) = {1, 7};
Line(6) = {7, 11};
Line(7) = {11, 6};
Line(8) = {6, 5};
Line(9) = {5, 8};
Line(10) = {8, 8};
Line(11) = {9, 8};
Line(12) = {7, 9};
Line(13) = {9, 10};
Line(14) = {10, 11};
Line(15) = {3, 6};
Line(16) = {4, 5};
Line(17) = {2, 8};

Line Loop(1) = {1, 15, -7, -6, -5};
Plane Surface(1) = {1};
Line Loop(2) = {4, 5, 12, 11, -17};
Plane Surface(2) = {2};
Line Loop(3) = {3, 17, -9, -16};
Plane Surface(3) = {3};
Line Loop(4) = {2, 16, -8, -15};
Plane Surface(4) = {4};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(5) = {5};
Line Loop(6) = {7, 8, 9, -11, 13, 14};
Plane Surface(6) = {6};
Line Loop(7) = {6, -14, -13, -12};
Plane Surface(7) = {7};
Surface Loop(1) = {5, 1, 4, 3, 2, 7, 6};
Volume(1) = {1};

Physical Volume("BODY") = {1};

Physical Surface("NORTH") = {3};
Physical Surface("SOUTH") = {1};
Physical Surface("WEST") = {2};
Physical Surface("EAST") = {4};
Physical Surface("BOTTOM") = {5};
Physical Surface("TOP") = {6};
Physical Surface("TOP_LOAD") = {7};
