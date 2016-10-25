/*Eric Gelphman
  University of California, San Diego Jacobs School of Engineering
  Department of Mechanical and Aerospace Engineering(MAE)*/
/*Matthew Uffenheimer
  University of California, Santa Barbara College of Creative Studies(CCS)*/
//October 19, 2016
/*MidnightOil v1.0-C++ program that numerically calculates the area of complicated shapes
in R^2 using a variation of Green's Theorem*/

#include <iostream>
#include <vector>
#include <string>
#include <ctime>
#include "exprtk.hpp"
using namespace std;

double DELTA = 0.05;//Step Size
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

/*Function structure contains function f(x) as well as the starting/ending
points of the traversal*/
struct functionStruct {
    string function;//f(x) that forms one side of shape
    point start;//Starting point of traversal
    point end;//End point of traversal
};

//Type definitions from ExprTk library
typedef exprtk::symbol_table<double> symbol_table_t;
typedef exprtk::expression<double> expression_t;
typedef exprtk::parser<double> parser_t;
expression_t expression;//Setting up evaluation infrastructure-no need to do this multiple times
parser_t parser;

/**********************Function Declarations**********************************/
void printPoint();//Function to print a strung representation of a point
vector<point> dfs(functionStruct, vector<point>);//Function to obtain points along boundary of shape
double eval(string, double);//Function to evaluate the function at a point (x,y)
double calcArea(vector<point>);//Function to calculate area

//Function to cprint a string representation of a point
void printPoint(point point1)
{
  printf("(%lf, %lf) ", point1.x, point1.y);
}

/*
Given function f(x) = 0 for segment of boundary of shape, obtain
points (x,y) along boundary of shape to give to then calculate area.
*/
vector<point> dfs(functionStruct fs1, vector<point> orderedPoints1)
{
    if((fs1.start.x > fs1.end.x && DELTA > 0) || (fs1.start.x < fs1.end.x && DELTA < 0))
      DELTA *= -1;
    point curPoint, end, next;
    int constantx = 0;
    curPoint = fs1.start;
    end = fs1.end;
    if(fs1.function.compare("constantx") == 0)
        constantx = 1;
    while(!(curPoint == end))//Doing traversal
    {
      if(constantx == 0)
      {
        next.x = curPoint.x + DELTA;//Obtain next search point
        next.y = eval(fs1.function, next.x);
      }
      else
      {
        if((fs1.start.y > fs1.end.y && DELTA > 0) || (fs1.start.y < fs1.end.y && DELTA < 0))
          DELTA *= -1;
        next.x = curPoint.x;
        next.y = curPoint.y + DELTA;
      }
      printPoint(next);
      orderedPoints1.push_back(next);//Add to storage
      curPoint = next;
    }
    return orderedPoints1;
}

//Evaluate a function f(x) at x = 4, return value of function at point (x,y)
double eval(string function, double a)
{
  symbol_table_t symbol_table;
  symbol_table.add_constants();
  symbol_table.add_variable("x", a);
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
    functionStruct fs1, fs2, fs3, fs4;
    point start1, start2, start3, start4;
    vector<point> orderedPoints;
    fs1.function = "constantx";
    fs2.function = "4/x";
    fs3.function = "constantx";
    fs4.function = "1/x";
    start1.x = 1.0;
    start1.y = 1.0;
    start2.x = 1.0;
    start2.y = 4.0;
    start3.x = 2.0;
    start3.y = 2.0;
    start4.x = 2.0;
    start4.y = 0.5;
    fs1.start = start1;
    fs1.end = start2;
    fs2.start = start2;
    fs2.end = start3;
    fs3.start = start3;
    fs3.end = start4;
    fs4.start = start4;
    fs4.end = start1;
    vector<functionStruct> functionVector;
    functionVector.push_back(fs1);
    functionVector.push_back(fs2);
    functionVector.push_back(fs3);
    functionVector.push_back(fs4);
    clock_t clk;
    clk = clock();
    int i;
    for(i = 0; i < functionVector.size(); i++)//Do search
    {
      functionStruct fs1 = functionVector[i];
      orderedPoints = dfs(fs1, orderedPoints);
    }
    double work = calcArea(orderedPoints);//Calculate area
    clk = clock() - clk;
    printf("The net work done by the engine is: %lf\n", work);
    printf("Size of vector: %lu\n", orderedPoints.size());
    printf("Runtime: %lf\n", ((double)clk) / CLOCKS_PER_SEC);
    /*double Qh, tCycle;
    printf("Enter the temperature of the energy absorbed by the engine to calculate efficiency:\n");
    scanf("%lf", &Qh);
    double efficiency = 1.0 - (work / Qh);
    printf("The efficiency of the engine is: %lf\n", efficiency);
    printf("Enter the time it takes to complete 1 cycle to calculate power\n");
    scanf("%lf", &power);
    double power = work / tCycle;
    printf("The power output of this engine is %lf watts", power);*/
    return 0;
}
