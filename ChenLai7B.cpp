/* Eric Gelphman
   University of California San Diego Department of Mathematics
*/
/* Matthew Uffenheimer
   University if California Santa Barbara College of Creative Studies (CCS)
*/
class pair {
	public:
		pair(double x_=0, double y_=0):x(x_),y(y_) {};
		double x;
		double y;
};

pair::operator +(pair p1, pair p2) {
	pair p(p2.x+p1.x, p2.y+p1.y);
	return p;
}

pair::operator -(pair p1, pair p2) {
	pair p(p2.x-p1.x, p2.y-p1.y);
	return p;
}

pair::operator /(pair p1, double a) {
	pair p(p1.x/a, p1.y/a);
	return p;
}


pair origin(0,0);

pair newton(f,pair p, pair d) {
Newton's method traditionally takes as input a function R->R and a point in R, and then finds a root of this function near that point. here we have a function f:R^2->R, a point p, and a unit direction vector d. We restrict f to the ray {p+dr : r in R+} which is essentially a function R+->R, because a line in R^2 is basically just a copy of R+ (R+ is the positive reals). Now we do newton's method as usual on this function.
}

pair anchor(pair p1, pair p2) {
returns the pair p such that (p,p1,p2) is a right isosceles triangle (ordered clockwise)
}

class LinkedList {
	public:
		void insert(pair p,node* n); //inserts a node whose data is "p" after the node "n"
}

LinkedList data;

void first_traversal() { //gets the intercepts of the function with the positive and negative x and y axes
	data.insert(newton(f,origin,pair(1,0)),head);
	data.insert(newton(f,origin,pair(0,1)),tail);
	data.insert(newton(f,origin,pair(-1,0)),tail);
	data.insert(newton(f,origin,pair(0,-1)),tail);
}

void refine(Node* n) {
	point p1 = *n;
	point p2 = *(n->next);
	p = anchor(p1,p2);
	point d = ((p1+p2)/2)-p; //p+d should be the midpoint of (p1,p2), i.e. d is the direction vector from p to the midpoint
	data.insert(newton(f,p,d),n); //calculate the root along the vector which bisects the angle (p1,p,p2). then insert it into the linked list after p1. this way, we can find an essential midpoint of p1 and p2 along the function's zero set.
}
