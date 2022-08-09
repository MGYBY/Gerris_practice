xb = DefineNumber[ 40.008, Name "loc" ];
width = DefineNumber[ 0.60, Name "r" ];
nd = DefineNumber[ 0.00798, Name "H" ];

leftCoord = xb-width/2.0;
rightCoord = leftCoord+width;
midCoord = leftCoord+width/2.0;

Point(1) = {leftCoord, 1.50, 0.0, 0.020};
Point(2) = {midCoord, 1.80, 0.0, 0.020};
Point(3) = {rightCoord, 1.50, 0.0, 0.020};
Point(4) = {midCoord, 1.20, 0.0, 0.020};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};

Extrude {0,0,nd*9.6} {
  Surface{6};
}

// Transfinite Line{1} = lr+1 Using Bump p;
