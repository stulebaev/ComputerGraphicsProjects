// triangulation.cpp

#pragma warning(disable:4530) // disable: "C++ exception handler used, but unwind semantics are not enabled. Specify /EHsc"

#include "delaunay.h"
#include <iostream>

int main(void) {
  double data[][2] = { {0.0,-0.5}, {0.5,-0.5}, {0.5,0.0}, {0.46194,-0.191342}, {0.353553,-0.353553}, {0.191342,-0.46194} };
  std::vector<Point2D> points;
  int k;
  for (k = 0; k < sizeof(data)/(sizeof(double)*2); k++)
    points.push_back(Point2D(data[k][0],data[k][1]));
  Delaunay triangulation(points);
  std::cout << "points:"  << std::endl << "[ ";
  for (k = 0; k < triangulation.npts; k++)
    std::cout << "[" << triangulation.pts[k].x << "," << triangulation.pts[k].y << "]; ";
  std::cout << "];" << std::endl;
  std::cout << "mesh:" << std::endl << "[ ";
  for (k = 0; k < triangulation.ntree; k++)
    if (triangulation.thelist[k].stat > 0)
      std::cout << "[" << triangulation.thelist[k].p[0] << " " << triangulation.thelist[k].p[1] << " " << triangulation.thelist[k].p[2] << "]; ";
  std::cout << "];" << std::endl;
  return 0;
}
