/*Eric Gelphman
University of California San Diego Jacobs School of Engineering
Department of Mechanical and Aerospace Engineering(MAE)*/
/*Matthew Uffenheimer
University of California Santa Barbara College of Creative Studies(CCS)
Mathematics Division*/
//September 3, 2016
//Pleiades v1.0
//Program to numerically do Green's Theorem - approximates the area of any convex shape

#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

/************************Structure Declarations***********************************************************/
//Point structure for a point in the plane(R^2) (x,y)
struct point{
    double x;//x-coordinate of point
    double y;//y-coordinate of point
};

point lowerLeft;//Lower left point of shape
string pointToString(point);//Function to create a string representation of a point
bool compare(point, point);//Compares two points
double distsqLL(point);//Calculates the square distance between a point and the lower left point
vector<point> modifiedGraham(vector<point>);//Function to create a polygon given a set of points
int findLL(vector<point>);//Function to find the position of the leftmost lowest point in the set
double chenLai(vector<point>);//Function to calculate area given an ordered list of points

//Function to create a string representation of a point
string pointToString(point point1)
{
  return "(" + to_string(point1.x) + "," + to_string(point1.y)+ "); ";
}

/*Function to compare to points. Returns 1 if p is greater than q, returns -1 if p is less than q. A point is greater than
 * another point if it has a greater theta(angular distance from lowerRight). If two points have the same theta, the point that is
 * closer to lowerRight(distance formula) is greater. Two points are considered equal if and only if their theta values are equal
 * and they are the same distance from lowerRight. In this case the function returns 0*/
bool compare(point p, point q)
{
    double thetaLLP = atan2(p.y - lowerLeft.y, p.x - lowerLeft.x);//Angle between p and lowerLeft
    double thetaLLQ = atan2(q.y - lowerLeft.y, q.x - lowerLeft.x);//Angle between q and lowerLeft
    if(thetaLLP > thetaLLQ)//If angle between lowerLeft and p is greater
        return false;
    else if(thetaLLP < thetaLLQ)//if angle between lowerLeft and q is greater
        return true;
    else//Angle between lowerLeft and p = angle between lowerLeft and q
    {
        double distLLP = distsqLL(p);
        double distLLQ = distsqLL(q);
        if(distLLP < distLLQ)//Distance between lowerLeft and q is greater
            return false;
        else if(distLLP > distLLQ)//Distance between lowerLeft and p is greater
            return true;
        else//Two identical points in same set
        {
            cout << "Error, two identical points in the same set!";
            abort();
        }
    }
}

//Function to calculate the square distance between a point and lowerLeft
double distsqLL(point p)
{
    return (p.x - lowerLeft.x) * (p.x - lowerLeft.x) + (p.y - lowerLeft.y) * (p.y - lowerLeft.y);
}

/*Function to order the points on the boundary of the shape in the correct order using a modified version of
Graham's Method*/
vector<point> modifiedGraham(vector<point> points)
{
    int lP = findLL(points);//Finding lowerLeftPos
    lowerLeft = points[lP];
    points.erase(points.begin() + lP);//Remove lowerLeft from set
    sort(points.begin(), points.end(), compare);//Sort points
    vector<point> boundaryPoints;//New vector containing all the points in order
    boundaryPoints.push_back(lowerLeft);//LowerLeft goes first
    int i;
    for(i = 0; i < points.size(); i++)//Adding other points in order
        boundaryPoints.push_back(points[i]);
    return boundaryPoints;
}

//Function to find the position of the lower leftmost point of the point set
int findLL(vector<point> points)
{
    int lowerLeftPos = 0;
    int size = points.size();
    int i;
    for(i = 1; i < size; i++) {//LinearSearch
        if (points[i].y < points[lowerLeftPos].y || (points[i].y == points[lowerLeftPos].y && points[i].x < points[lowerLeftPos].x))
            lowerLeftPos = i;
    }
    return lowerLeftPos;
}

/*Function that actually calculates area, given points along boundary of shape
using variation of Green's Theorem*/
double chenLai(vector<point> orderedPoints)
{
    double area = 0.0;
    int idx = 0;//Index in vector
    int numItr = 0;//Number of iterations(additions). Actual formula is defined using series
    int size = orderedPoints.size();
    while(numItr < size)
    {
        if(numItr < size - 1)
        {
            area += ((orderedPoints[idx].x * orderedPoints[idx + 1].y) - (orderedPoints[idx].y * orderedPoints[idx + 1].x));
        }
        else
        {
            area += ((orderedPoints[idx].x * orderedPoints[0].y) - (orderedPoints[0].x * orderedPoints[idx].y));
        }
        idx++;
        numItr++;
    }
    if(area < 0)//Area cannot be negative
    {
        area *= -1;
    }
    return (area / 2);
}

int main()
{
    vector<point> boundaryPoints;
    while(1)
    {
        string input;
        cout << "Enter point on boundary of shape in this format: (x,y), and enter end when done.\n";
        cin >> input;
        if(input.compare("end") != 0)
        {
            point point1;
            int commaPos = input.find(',');
            point1.x = stod(input.substr(1, commaPos));
            point1.y = stod(input.substr(commaPos + 1, input.length() - 1));
            boundaryPoints.push_back(point1);
        }
        else
        {
            break;
        }
    }
    vector<point> points = modifiedGraham(boundaryPoints);
    int i;
    for(i = 0; i < points.size(); i++)
    {
        cout << pointToString(points[i]);
        cout << "\n";
    }
    double area = chenLai(points);
    cout << "The area of the shape formed by the points is: ";
    cout << area;
    cout << "\n";
    return 0;
}
