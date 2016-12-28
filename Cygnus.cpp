/*Eric Gelphman
  University of California, San Diego Department of Physics
  Matthew Uffenheimer
  University of California, Santa Barbara College of Creative Studies(CCS)
December 27, 2016
Cygnus -The seventh major iteration of our Algorithm to numerically do Green's Theroem*/

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <ctime>
#include "exprtk.hpp"
using namespace std;

double STEP_SIZE = 0.1;//Step Size
double h = 0.0000001;//Needed for numerical differentiation
double EPSILON = 0.00001;//Epsilon needed for operations with doubles
/*****************************Structure Definitions**************************/
//Point structure for a point (x,y)
struct point{
    double x;//x-coordinate
    double y;//y-coordinate
};

//Overload == operator for use with points
inline bool operator == (point const& p, point const& q)
{
    return (abs(p.x - q.x) <= EPSILON && abs(p.y - q.y) <= EPSILON);
}

/*Function structure contains function f(x,y) as well as the minimum and maximum x-
and y- coordinates for which the function is defined*/
struct functionStruct {
    string function;//f(x,y) that forms boundary of shape
    point start;//Starting point of traversal
    point end;//End point of traversal. For a closed loop, start = end
    double xmin;//Min. x-coordinate of function
    double xmax;//Max. x-coordinate of function
    double ymin;//Min. y-coordinate of function
    double ymax;//Max. y-coordinate of function
};

vector<point> orderedPoints;//Stores points in order
//Type definitions from ExprTk library
typedef exprtk::symbol_table<double> symbol_table_t;
typedef exprtk::expression<double> expression_t;
typedef exprtk::parser<double> parser_t;
expression_t expression;//Setting up evaluation infrastructure-no need to do this multiple times
parser_t parser;

/**********************Function Declarations**********************************/
void printPoint();//Function to print a strung representation of a point
void traversal(functionStruct);//Function to obtain points along boundary of shape
point sfd2(point, functionStruct);//Function to determine the next point in the search
int inBounds(point, functionStruct);//Function to determine if a point is within the overall boundaries of the search
double eval(string, double, double);//Function to evaluate the function at a point (x,y)
double *numericalGrad(string, double, double);//Function to numerically calculate the partial derivatives a two-variable function f(x,y)
point newton(string, double, double, double);//Newton's Method
double calcArea(vector<point>);//Function to calculate area

//Function to print a string representation of a point
void printPoint(point point1)
{
  printf("(%lf, %lf) ", point1.x, point1.y);
}

/*
Given function f(x,y) = 0 for boundary(or segment of boundary of shape), obtain
points (x,y) along boundary of shape to then calculate area.
*/
void traversal(functionStruct fs1)
{
    point curPoint, end, next;
    curPoint = fs1.start;
    end = fs1.end;
    orderedPoints.push_back(curPoint);//Adding starting and first point to storage
    curPoint = sfd2(curPoint, fs1);
    orderedPoints.push_back(curPoint);
    while(!(curPoint == end))//Doing traversal
    {
      next = sfd2(curPoint, fs1);//Obtain next search point
      orderedPoints.push_back(next);//Add to storage
      curPoint = next;
    }
}

/*Function that determines the next point on the boundary to traverse by finging the tangent line
to the boundary curve at that point, moving a fixed step distance along the tangent line,
comuting the normal to the line at that point, and then finding where the normal line and
the boundary curve intersect. This is called the SirFrancisDrake Algorithm*/
point sfd2(point curPoint, functionStruct f1)
{
    point returned;
    double *gradPtr = numericalGrad(f1.function, curPoint.x, curPoint.y);//Obtain partial derivatives
    double x_prime = curPoint.x + (*(gradPtr + 1) * STEP_SIZE);//Use pointer arithmetic
    double y_prime = curPoint.y + ((-1 * (*gradPtr)) * STEP_SIZE);//Move a step distance down the tangent line to point (x_prime, y_prime)
    double *gradPtr_prime = numericalGrad(f1.function, x_prime, y_prime);//Compute gradient at this point to determine line H
    double M = (*(gradPtr_prime + 1) - y_prime) / (*gradPtr_prime - x_prime);
    double B = y_prime - M * x_prime;
    return newton(f1.function, M, B, x_prime);//Find intersection between line H and function F, return result
}

//Function to determine if a point is within the overall boundaries of the search
int inBounds(point p, functionStruct grid1)
{
    if (p.x >= grid1.xmin && p.x <= grid1.xmax && p.y >= grid1.ymin && p.y <= grid1.ymax)
        return 1;
    return 0;
}

//Function to numerically calculate the gradient of a function f(x,y) at point (a,b)
double* numericalGrad(string function, double a, double b) {
    double partials[2];
    double *returnPtr;
    returnPtr = partials;
    partials[0] = (eval(function, a + h, b) - eval(function, a - h, b)) / (2 * h);//df/dx
    partials[1] = (eval(function, a, b + h) - eval(function, a, b - h)) / (2 * h);//df/dy
    return returnPtr;//Return as array
}

//Evaluate a function f(x,y) at a point (x,y), return value of function at point (x,y)
double eval(string function, double a, double b)
{
  symbol_table_t symbol_table;
  symbol_table.add_constants();
  symbol_table.add_variable("x", a);
  symbol_table.add_variable("y", b);
  expression.register_symbol_table(symbol_table);
  if(!(parser.compile(function, expression)))//If f(x,y) is not a valid expression that can be evaluated by ExprTk
  {
    printf("Error: %s\tExpression: %s\n", parser.error().c_str(), function.c_str());
    exit(1);
  }
  double result = expression.value();
  return result;
}

//Function to find the intersection of Boundary curve F and line H using Newton's Method
point newton(string function, double M, double B, double x_init)
{
    double x_i1, x_i, H, fxh, dfdx;
    int n = 0;
    x_i1 = x_init;
    while(abs(x_i1 - x_i) >= EPSILON && n <= 8)//2 Constraints: epsilon and number of iterations
    {
        H = M * x_i1 + B;
        fxh = eval(function, x_i1, H);//y is function of x, sub in for y, reducing to a nonlinear equation in 1 variable
        dfdx = (eval(function, x_i1 + h, H) - eval(function, x_i1 - h, H)) / (2 * h);//Perfom Newton's Method
        x_i = x_i1;
        x_i1 -= fxh/dfdx;
        n += 1;
    }
    point returned;
    returned.x = x_i;
    returned.y = H;
    return returned;//Return point
}


/*Function that actually calculates area, given points along boundary of shape
using variation of Green's Theorem*/
double calcArea(vector<point> orderedPoints)
{
  double area = 0.0;
  int idx = 0;//Index in vector
  int numItr = 0;//Number of iterations(additions). Actual formula is defined using series
  int size = orderedPoints.size();
  while(numItr < size)
  {
      if(numItr < size - 1)
          area += ((orderedPoints[idx].x * orderedPoints[idx + 1].y) - (orderedPoints[idx].y * orderedPoints[idx + 1].x));
      else
          area += ((orderedPoints[idx].x * orderedPoints[0].y) - (orderedPoints[0].x * orderedPoints[idx].y));
      idx++;
      numItr++;
  }
  if(area < 0)//Area cannot be negative
      area *= -1;
  return (area / 2);
}

int main() {
  functionStruct fs0;
  fs0.function = "(x*x)+(x*y)+(y*y)-4";//Shape(Ellipse)
  fs0.xmin = -2.31;
  fs0.xmax = 2.31;
  fs0.ymin = -2.31;
  fs0.ymax = 2.31;
  point start1;
  start1.x = 0.0;
  start1.y = 2.0;
  fs0.start = start1;
  fs0.end = start1;
  vector<functionStruct> functionVector;
  functionVector.push_back(fs0);
  clock_t clk;
  clk = clock();
  int i;
  for(i = 0; i < functionVector.size(); i++)//Do search
  {
    functionStruct fs1 = functionVector[i];
    traversal(fs1);
  }
  double area = calcArea(orderedPoints);//Calculate area
  clk = clock() - clk;
  printf("The area of the shape is: %lf\n", area);
  printf("Size of vector: %lu\n", orderedPoints.size());
  printf("Runtime: %lf\n", ((double)clk) / CLOCKS_PER_SEC);
  return 0;
}
