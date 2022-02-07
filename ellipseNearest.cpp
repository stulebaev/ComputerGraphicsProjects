/* ellipseNearest.cpp */

#include <cmath>
#include <cstdio>

int findEllipseNearestPoint(double xP, double yP, double a, double b, double* x, double* y)
{
  const int MAX_ITERATIONS = 50; /* maximum number of interations */
  const double EPSILON = 1.0e-6; /* tolerance */

  double t = atan2(a*yP, b*xP); /* initial approximation */
  double ba = (b*b-a*a);
  int k = 1;
  for (; ;) {
    double cost = cos(t), sint = sin(t);
    double f = ba*cost*sint + a*xP*sint - b*yP*cost; /* f(t) */
    double df = ba*(cost*cost-sint*sint) + a*xP*cost + b*yP*sint; /* f'(t) */
    double delta = f/df;
    t -= delta;
    if (fabs(delta) < EPSILON) break;
    k++;
    if (k > MAX_ITERATIONS) return -1;
  }
  *x = a*cos(t);
  *y = b*sin(t);
  return k;
}

int main()
{
  double a = 10.0, b = 5.0; /* ellipse parameters */
  double xP = -10.0, yP = 6.0; /* given point */
  double x, y;

  printf("Number of interations: %d\n", findEllipseNearestPoint(xP,yP, a,b, &x,&y));
  printf("result: (%f,%f)\n", x,y);

  return 0;
}
/*
\documentclass{article}
\usepackage[utf8]{inputenc} 
\usepackage[russian]{babel}

\begin{document}
Эллипс $E(t) = (a\,cos t, b\,sin t)$

Найти точку эллипcа, ближайшую к заданной точке $P(x_p, y_p)$
\bigskip
\\
Задача сводится к минимизации нормы
$$
||E(t)-P||^2 \to min
$$
которая приводит к уравнению:
$$
\frac{d}{dt}\left[ (a\,cos t-x_p)^2 + (b\,sin t-y_p)^2 \right] = 0
$$
\\
Нелинейное уравнение
$$
f(t) = (b^2 - a^2) cos t\,sin t + a x_p sin t - b y_p cos t = 0
$$
можно решить итерационным методом Ньютона:
$$
t_{k+1} = t_k - \frac{f(t_k)}{f'(t_k)}
$$
где
$$
f'(t) = (b^2 - a^2) (cos^2 t - sin^2 t) + a x_p cos t + b y_p sin t
$$
\\
В качестве начального приближения можно выбрать точку пересечения отрезка $PO$ с эллипсом:
$$
t_0 = \arctg\big(\frac{ay_p}{bx_p}\big) = \tt{atan2(ay_p,bx_p)}
$$
\end{document}
*/
