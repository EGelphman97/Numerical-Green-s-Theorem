/*Eric Gelphman
  University of California, San Diego Department of Physics
  Matthew Uffenheimer
  University of California, Santa Barbara College of Creative Studies(CCS)
December 27, 2016
Cygnus -The seventh major iteration of our Algorithm to numerically do Green's Theroem*/

#include <iostream>
#include <vector>
#include <util>
#include <cmath>
#include <string>
#include <ctime>
#include "exprtk.hpp"
using namespace std;

double STEP_SIZE = 0.1;//Step Size
double h = 0.0000001;//Needed for numerical differentiation
double EPSILON = 0.00001;//Epsilon needed for operations with doubles

//Overload == operator for use with pairs of doubles
inline bool operator == (pair<double,double> const& p, pair<double,double> const& q)
{
    return (abs(p.first - q.first) <= EPSILON && abs(p.second - q.second) <= EPSILON);
}

/*Function structure contains function f(x,y) as well as the minimum and maximum x-
and y- coordinates for which the function is defined*/
struct functionStruct {
    string function;//f(x,y) that forms boundary of shape
    pair<double,double> start;//Starting point of traversal
    pair<double,double> end;//End point of traversal. For a closed loop, start = end
};

vector<pair<double,double>> orderedPoints;//Stores points in order
//Type definitions from ExprTk library
typedef exprtk::symbol_table<double> symbol_table_t;
typedef exprtk::expression<double> expression_t;
typedef exprtk::parser<double> parser_t;
expression_t expression;//Setting up evaluation infrastructure-no need to do this multiple times
parser_t parser;
symbol_table_t symbol_table;

/**********************Function Declarations**********************************/
void printPoint();//Function to print a strung representation of a point
void traversal(functionStruct);//Function to obtain points along boundary of shape
pair<double,double> getPoint(pair<double,double>, functionStruct);//Function to determine the next point in the search
double eval(pair<double,double>, double, double);//Function to evaluate the function at a point (x,y)
double *numericalGrad(string, double, double);//Function to numerically calculate the partial derivatives a two-variable function f(x,y)
double calcArea(vector<pair<double,double>>);//Function to calculate area

//Evaluate a function f(x,y) at a point (x,y), return value of function at point (x,y)
double eval(string function, double a, double b)
{
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

/*
Given function f(x,y) = 0 for boundary(or segment of boundary of shape), obtain
points (x,y) along boundary of shape to then calculate area.
*/
void traversal(functionStruct fs1)
{
    pair<double,double> curPoint, end, next;
    curPoint = fs1.start;
    end = fs1.end;
    orderedPoints.push_back(curPoint);//Adding starting and first point to storage
    curPoint = getPoint(curPoint, fs1);
    orderedPoints.push_back(curPoint);
    while(!(curPoint == end))//Doing traversal
    {
      next = getPoints(curPoint, fs1);//Obtain next search point
      orderedPoints.push_back(next);//Add to storage
      curPoint = next;
    }
}

/*Function that determines the next point on the boundary to traverse by finging the tangent line
to the boundary curve at that point, moving a fixed step distance along the tangent line,
comuting the normal to the line at that point, and then finding where the normal line and
the boundary curve intersect. This is called the SirFrancisDrake Algorithm*/
pair<double,double> getPoint(pair<double,double> curPoint, functionStruct f1)
{
    /*REDO*/
}

//Function to numerically calculate the gradient of a function f(x,y) at point (a,b)
double* numericalGrad(string function, double a, double b)
{
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

/*Function that actually calculates area, given points along boundary of shape
using variation of Green's Theorem*/
double calcArea(vector<pair<double,double>> orderedPoints)
{
  double area = 0.0;
  int idx = 0;//Index in vector
  int numItr = 0;//Number of iterations(additions). Actual formula is defined using series
  int size = orderedPoints.size();
  while(numItr < size)
  {
      if(numItr < size - 1)
          area += ((orderedPoints[idx].first * orderedPoints[idx + 1].second) - (orderedPoints[idx].second * orderedPoints[idx + 1].first));
      else
          area += ((orderedPoints[idx].first * orderedPoints[0].second) - (orderedPoints[0].first * orderedPoints[idx].second));
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
  pair<double,double> start1;
  start1.first = 0.0;
  start1.second = 2.0;
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
