#include <cmath>

///////////////////////////////////////////////////////////////////////////////
//  A vector in a plane:              
///////////////////////////////////////////////////////////////////////////////
struct point3d
{
	double x, y, z;
};

inline point3d operator-(point3d& b, point3d& a)
{
	point3d ab = {b.x - a.x,	b.y - a.y,	b.z - a.z};
	return ab;
}
inline point3d operator+(point3d& b, point3d& a)
{
	point3d ab = {b.x + a.x,	b.y + a.y,	b.z + a.z};
	return ab;
}
inline point3d operator/(point3d& a, double r)
{
	point3d n = {a.x / r,	a.y / r, a.z / r};
	return n;
}
inline point3d operator*(point3d& a, double r)
{
	point3d n = {a.x * r,	a.y * r, a.z * r};
	return n;
}
///////////////////////////////////////////////////////////////////////////////
//  3d dot product              
///////////////////////////////////////////////////////////////////////////////
inline double dot3d(point3d& v1, point3d& v2)
{
  return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}
///////////////////////////////////////////////////////////////////////////////
// distance between two points
///////////////////////////////////////////////////////////////////////////////
inline double r(point3d& a, point3d& b)
{
	point3d ab = b - a;
	double r2 = dot3d(ab,ab);
	return sqrt(r2);
}
///////////////////////////////////////////////////////////////////////////////
// distance from (0,0)
///////////////////////////////////////////////////////////////////////////////
inline double r(point3d& a)
{
	return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}
