# Numerical Green's Theorem
This repository is for a project I made to numerically calculuate the area enclosed by complicated two-dimensional curves using a variation of Green's Theorem. There are currently 4 different versions of this project, which Matt and I have named 
Project Los Gatos, after the town in which our high school is located.

The General Algorithm is the same for all versions:
1. Find points on the curve f(x,y) = 0 that forms the boundary of the shape
2. Order the points correctly so as to generate a parametrization for the boundary
3. Plud the correctly ordered points into the polygonal area formula(implemeneted in the prototype-pa3.c) to numerically
calculate the area enclosed by the curve)

Pleides was the first version that worked successfully. Pleides is a special case where the points are given and the shape they form is convex

Hyades is the general case where only the function(or functions-if the boundary is piecewise defined) that forms the boudnary of the shape is given. It uses a traversal algortihm called BlackBird to find points on the boundary curve in the correct order.

MidnightOil is a variation of Pleides designed to perform thermodynamic calculations, specifically with regards to engine power and efficiency.

Cygnus(currently in development) is a new version of Hyades that is designed to be more robust as well as to integrate some 
functionality from MidnightOil into the general algorithm. It also continues the trend of naming things after celestial objects.

The powerpoint contained in this repository is for a presentation I made to UCSD's Math Department's undergraduate student colloqium about the project. I was the first undergarduate in several years to present his or her own research at the colloqium.

