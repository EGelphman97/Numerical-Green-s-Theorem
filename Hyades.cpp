/*Eric Gelphman
  University of California, San Diego Jacobs School of Engineering
  Department of Mechanical and Aerospace Engineering(MAE)*/
/*Matthew Uffenheimer
  University of California, Santa Barbara College of Creative Studies(CCS)*/
//October 19, 2016
/*Hyades v1.0-C++ program that numerically calculates the area of complicated shapes
in R^2 using a variation of Green's Theorem*/

#include <iostream>
#include <vector>
#include <cmath>
#include <unordered_set>
#include <string>
#include <ctime>
#include "exprtk.hpp"
using namespace std;

double DELTA = 0.05;//Step Size
double h = 0.0000001;//Needed for numerical differentiation
double EPSILON = 0.00001;//Epsilon needed for operations with doubles
//Search order: up, down, left, right, upper left, lower right, upper right, lower left
double xc[] = {0, 0, -1 * DELTA, DELTA, -1 * DELTA, DELTA, DELTA, -1 * DELTA};//Search grid
double yc[] = {DELTA, -1 * DELTA, 0, 0, DELTA, -1 * DELTA, DELTA, -1 * DELTA};

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

namespace std
{
    //Implements R17L23XOR Hash Algorithm
    template <>
    struct hash<point>
    {
        size_t operator()(point const &point1) const{
          unsigned int id1, id2;
          if(point1.x > 0)
            id1 = point1.x * 10000010;
          else
            id1 = point1.x * 10000000;
          if(point1.y > 0)
            id2 = point1.y * 10000010;
          else
            id2 = point1.y * 10000000;
          id1 = id1 >> 17;
          id1 = id1 << 23;
          id2 = id2 >> 17;
          id2 = id2 << 23;
          return ((id1 ^ id2));
        }
    };
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
unordered_set<point, hash<point> > visitedPoints;//Hash table to store visited points
//Type definitions from ExprTk library
typedef exprtk::symbol_table<double> symbol_table_t;
typedef exprtk::expression<double> expression_t;
typedef exprtk::parser<double> parser_t;
expression_t expression;//Setting up evaluation infrastructure-no need to do this multiple times
parser_t parser;

/**********************Function Declarations**********************************/
void printPoint();//Function to print a strung representation of a point
void dfs(functionStruct);//Function to obtain points along boundary of shape
point blackBirdN(point, functionStruct);//Function to determine the next point in the search
int inBounds(point, functionStruct);//Function to determine if a point is within the overall boundaries of the search
double eval(string, double, double);//Function to evaluate the function at a point (x,y)
double *numericalPartialDiff(string, double, double);//Function to numerically calculate the partial derivatives a two-variable function f(x,y)
double calcArea(vector<point>);//Function to calculate area

//Function to cprint a string representation of a point
void printPoint(point point1)
{
  printf("(%lf, %lf) ", point1.x, point1.y);
}

/*
Given function f(x,y) = 0 for boundary(or segment of boundary of shape), obtain
points (x,y) along boundary of shape to give to chenLai to calculate area.
*/
void dfs(functionStruct fs1)
{
    point curPoint, end, next;
    curPoint = fs1.start;
    end = fs1.end;
    orderedPoints.push_back(curPoint);//Adding starting and first points to storage
    visitedPoints.insert(curPoint);
    curPoint = blackBirdN(curPoint, fs1);
    orderedPoints.push_back(curPoint);
    visitedPoints.insert(curPoint);
    while(!(curPoint == end))//Doing traversal
    {
      next = blackBirdN(curPoint, fs1);//Obtain next search point
      printPoint(next);
      orderedPoints.push_back(next);//Add to storage
      visitedPoints.insert(next);
      curPoint = next;
      if(visitedPoints.load_factor() >= 0.67)//Rehash if necessary
          visitedPoints.rehash(visitedPoints.bucket_count() * 2);
      if(orderedPoints.size() == 25)
          visitedPoints.erase(fs1.start);
    }
}

/*Function that determines the next point on the boundary to traverse by identifying the point closest to the linearization
 of the boundary function at prevPoint. This is called the Blackbird Algorithm*/
point blackBirdN(point curPoint, functionStruct grid1)
{
    point returned;
    double *partialPtr = numericalPartialDiff(grid1.function, curPoint.x, curPoint.y);//Obtain partial derivatives
    double A = *partialPtr;//Use pointer arithmetic
    double B = *(partialPtr + 1);
    double C = (*(partialPtr + 2)) - (curPoint.x * A + curPoint.y * B);
    double minDist = HUGE_VAL;
    int i;
    for (i = 0; i < 8; i++)//Iterate through search grid
    {
        point tPoint;
        tPoint.x = curPoint.x + xc[i];
        tPoint.y = curPoint.y + yc[i];
        if(visitedPoints.count(tPoint) == 0 && inBounds(tPoint, grid1) == 1)//If point is in bounds and has not been visited before
        {
              double tDist = abs(A * tPoint.x + B * tPoint.y + C) / sqrt((A * A) + (B * B));//Using formula to find distance between point and line
              if (tDist < minDist)//Closest point to line becomes next point in traversal
              {
                  minDist = tDist;
                  returned = tPoint;
              }
        }
    }
    return returned;
}

//Function to determine if a point is within the overall boundaries of the search
int inBounds(point p, functionStruct grid1)
{
    if (p.x >= grid1.xmin && p.x <= grid1.xmax && p.y >= grid1.ymin && p.y <= grid1.ymax)
        return 1;
    return 0;
}

//Function to numerically calculate the first-order partial derivatives of a function f(x,y) at point (a,b)
double* numericalPartialDiff(string function, double a, double b) {
    double partials[3];
    double *returnPtr;
    returnPtr = partials;
    partials[0] = (eval(function, a + h, b) - eval(function, a - h, b)) / (2 * h);//df/dx
    partials[1] = (eval(function, a, b + h) - eval(function, a, b - h)) / (2 * h);//df/dy
    partials[2] = eval(function, a, b);//Value of f at (a,b)
    return returnPtr;
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
    exit(0);
  }
  double result = expression.value();
  return result;
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
  {
      area *= -1;
  }
  return (area / 2);
}

int main() {
    functionStruct fs0, fs1;
    fs0.function = "x*x";
    fs1.function = "2*x";
    fs0.xmin = 0.0;
    fs0.xmax = 2.0;
    fs0.ymin = 0.0;
    fs0.ymax = 4.0;
    fs1.xmin = 0.0;
    fs1.xmax = 2.0;
    fs1.ymin = 0.0;
    fs1.ymax = 4.0;
    point start1, start2;
    start1.x = 0.0;
    start1.y = 0.0;
    start2.x = 2.0;
    start2.y = 4.0;
    fs0.start = start1;
    fs0.end = start2;
    fs1.start = start2;
    fs1.end = start1;
    vector<functionStruct> functionVector;
    functionVector.push_back(fs0);
    functionVector.push_back(fs1);
    clock_t clk;
    clk = clock();
    int i;
    for(i = 0; i < functionVector.size(); i++)//Do search
    {
      functionStruct fs1 = functionVector[i];
      dfs(fs1);
    }
    double area = calcArea(orderedPoints);//Calculate area
    clk = clock() - clk;
    printf("The area of the shape is: %lf\n", area);
    printf("Size of vector: %lu\n", orderedPoints.size());
    printf("Runtime: %lf\n", ((double)clk) / CLOCKS_PER_SEC);
    return 0;
}
