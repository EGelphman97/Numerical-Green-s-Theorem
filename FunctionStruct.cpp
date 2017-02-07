//Eric Gelphman
//FunctionStruct class declaration for Cygnus
//February 7, 2017

#include <util>
#include <iostream>
#include <string>
#include "exprk.cpp"
using namespace std;

typedef std::pair<double, double> point;

/* Overload == operator for use with pairs of doubles(points)
 * Equality is based on distance between points, if distance <= epsilon,
 * points are considered equal. If distance >= epsilon, points are not
 * considered equal
 */
inline bool operator == (point const& p, point const& q)
{
    return sqrt(pow(p.first - q.first, 2) + pow(p.second -q.second, 2)) <= EPSILON;
}

inline point operator + (point const& p, point const& q)
{
    return point (p.first + q.first, p.second + q.second);
}

inline point operator - (point const& p, point const& q)
{
  return point (p.first - q.first, p.second - q.second);
}

class functionStruct
{
    string function;
    point start;
    point end;
}

void functionStruct::init()
{
    cout << "Enter a function f(x,y) = 0 or f(x) = 0\n";
    cin >> this->function;
    cout << "Enter the starting point (x,y) for this traversal\n";
    /* TODO */
}
