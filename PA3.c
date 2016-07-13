//Eric Gelphman
//May 16, 2016
//PA3 v1.0
/*This program calculates the area of any polygon in the plane (R^2) given
its verticies using Green's Theorem. The greens function actually does the
computation*/

#include <stdio.h>

double greens(double[]);//Function declaration, coordinates can be decimals
int numVerticies;//Number of verticies polygin has
double area;//Area of polygon

int main()
{
  printf("Enter the number of verticies the polygon has\n");
  scanf("%d", &numVerticies);
  /*Each element in coordinate (x or y) gets 1 cell in array, so array needs
  to be of size numVerticies * 2 */
  double verticies[numVerticies * 2];
  int idx = 0;
  while(idx < (numVerticies * 2))//Getting inputs
  {
    if(idx % 2 == 0)
      printf("Enter x-coordinate of vertex %d \n", (idx / 2) + 1);
    else
      printf("Enter y-coordinate of vertex %d \n", (idx  / 2) + 1);
    scanf("%lf", &verticies[idx]);
    idx++;
  }
  area = greens(verticies);//Calculating area
  printf("The area of the polygon with %d verticies is %lf \n", numVerticies, area);
  return 0;
}

//Function that does actual computation using Green's Theorem
double greens(double v[])
{
  double area = 0.0;
  int idx = 0;//Index in array
  int numItr = 0;//Number of iterations(additions). Actual formula is defined using series
  while(numItr < numVerticies)
  {
    if(numItr < numVerticies - 1)
      area += ((v[idx] * v[idx + 3]) - (v[idx + 2] * v[idx + 1]));
    else
      area += ((v[idx] * v[1]) - (v[0] * v[idx + 1]));
    idx += 2;
    numItr = numItr + 1;
  }
  if(area < 0)//Area cannot be negative
    area *= -1;
  return (area / 2);
}
