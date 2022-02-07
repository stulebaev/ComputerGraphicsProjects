#include <math.h>

#define TRUE 1
#define FALSE 0

#define CALC_ACCURACY 1e-7 // точность проверки совпадения с нулем

double Determinant4(double* m) {
/* Computes the determinant of the matrix of this type:
   |m0  m1  m2  1|
   |m3  m4  m5  1| =
   |m6  m7  m8  1|
   |m9  m10 m11 1|
 m0*(m4*(m8-m11)+m7*(m11-m5)+m10*(m5-m8))+
 m3*(m1*(m11-m8)+m7*(m2-m11)+m10*(m8-m2))+
 m6*(m1*(m5-m11)+m4*(m11-m2)+m10*(m2-m5))+
 m9*(m1*(m8-m5) +m4*(m2-m8) +m7 *(m5-m2))
Input parameter *m points to memory buffer that contains m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11
*/
  double c0 = m[4]*(m[8]-m[11]) + m[7]*(m[11]-m[5]) + m[10]*(m[5]-m[8]);
  double c3 = m[1]*(m[11]-m[8]) + m[7]*(m[2]-m[11]) + m[10]*(m[8]-m[2]);
  double c6 = m[1]*(m[5]-m[11]) + m[4]*(m[11]-m[2]) + m[10]*(m[2]-m[5]);
  double c9 = m[1]*(m[8]-m[5])  + m[4]*(m[2]-m[8])  + m[7] *(m[5]-m[2]);
  return (m[0]*c0 + m[3]*c3 + m[6]*c6 + m[9]*c9);
}

int PointInTetrahedron(double* point, double* tetrahedron) {
/* Determines if a point is within a tetrahedron */

/*    |x1 y1 z1 1|
 D0 = |x2 y2 z2 1|
      |x3 y3 z3 1|
      |x4 y4 z4 1| */
  double D0, D1, D2, D3, D4;
  double tmp[3*4];
  D0 = Determinant4(tetrahedron);
  if (fabs(D0) < CALC_ACCURACY) { /* all points are coplanar (tetrahedron is degenerate) */
    return 0;
  }
/*     |x1 y1 z1 1|
  D4 = |x2 y2 z2 1|
       |x3 y3 z3 1|
       |x  y  z  1| */
  memcpy(tmp, tetrahedron, 9*sizeof(double));
  memcpy(&tmp[9], point, 3*sizeof(double));
  D4 = Determinant4(tmp);
/*     |x  y  z  1|
  D1 = |x2 y2 z2 1|
       |x3 y3 z3 1|
       |x4 y4 z4 1| */
  memcpy(&tmp[9], &tetrahedron[9], 3*sizeof(double));
  memcpy(&tmp[0], point, 3*sizeof(double));
  D1 = Determinant4(tmp);
/*     |x1 y1 z1 1|
  D2 = |x  y  z  1|
       |x3 y3 z3 1|
       |x4 y4 z4 1| */
  memcpy(&tmp[0], &tetrahedron[0], 3*sizeof(double));
  memcpy(&tmp[3], point, 3*sizeof(double));
  D2 = Determinant4(tmp);
/*     |x1 y1 z1 1|
  D3 = |x2 y2 z2 1|
       |x  y  z  1|
       |x4 y4 z4 1| */
  memcpy(&tmp[3], &tetrahedron[3], 3*sizeof(double));
  memcpy(&tmp[6], point, 3*sizeof(double));
  D3 = Determinant4(tmp);
  if ((D0*D1 > 0) && (D0*D2 > 0) && (D0*D3 > 0) && (D0*D4 > 0)) return TRUE;
  else return FALSE;
}

#include <stdio.h>
int main(void) {
  double point[] = { 0.2, -0.5, 0.1 };
  double tetrahedron[] = {
    0.2502119197232506, -0.5232796237309835, 0.2,
    0.236578018614319,  -0.5080378874515262, 0.0,
    0.09283243855059381,-0.4996528869630563, 0.0,
    0.1715723661539595, -0.386120805682164,  0.0
  };
  printf("%d\n", PointInTetrahedron(point, tetrahedron));
  return 0;
}
