// delaunay.h
#pragma once

#include <cmath> //for sqrt()
#include <vector>

typedef unsigned __int64 Ullong;

// macro-like inline functions
template<class T>
inline T SQR(const T a) { return a*a; }
template<class T>
inline void SWAP(T &a, T &b) { T dum=a; a=b; b=dum; }

// exception handling
#define NRthrow(message) \
{ printf("ERROR: %s\n     in file %s at line %d\n", message,__FILE__,__LINE__); throw(1); }

struct Point2D {
//Simple structure to represent a 2D point
  double x,y; //The coordinates
  Point2D(const Point2D &p) { //Copy constructor
    x = p.x; y = p.y;
  }
  Point2D operator= (const Point2D &p) { //Assigment operator
    x = p.x; y = p.y;
    return *this;
  }
  Point2D operator== (const Point2D &p) const {
    if ((x != p.x) || (y != p.y)) return false;
    return true;
  }
  Point2D(double _x = 0.0, double _y = 0.0) { //Construct by coordinate value
    x = _x; y = _y;
  }
};

struct Circle {
  Point2D center;
  double radius;
  Circle(const Point2D &cen, double rad) : center(cen), radius(rad) {}
};

Circle circumcircle(Point2D a, Point2D b, Point2D c) {
  double a0,a1,c0,c1,det,asq,csq,ctr0,ctr1,rad2;
  a0 = a.x - b.x; a1 = a.y - b.y;
  c0 = c.x - b.x; c1 = c.y - b.y;
  det = a0*c1 - c0*a1;
  if (det == 0.0) NRthrow("no circle thru colinear points");
  det = 0.5/det;
  asq = a0*a0 + a1*a1;
  csq = c0*c0 + c1*c1;
  ctr0 = det*(asq*c1 - csq*a1);
  ctr1 = det*(csq*a0 - asq*c0);
  rad2 = ctr0*ctr0 + ctr1*ctr1;
  return Circle(Point2D(ctr0 + b.x, ctr1 + b.y), sqrt(rad2));
}

struct Ranhash {
  inline Ullong int64(Ullong u) {
    Ullong v = u * 3935559000370003845LL + 2691343689449507681LL;
    v ^= v >> 21; v ^= v << 37; v ^= v >> 4;
    v *= 4768777513237032717LL;
    v ^= v << 20; v ^= v >> 41; v ^= v << 5;
    return  v;
  }
  inline unsigned int int32(Ullong u) {
    return (unsigned int)(int64(u) & 0xffffffff);
  }
  inline double real(Ullong u) {
    return 5.42101086242752217E-20 * int64(u);
  }
};

template<class keyT, class hfnT> struct Hashtable {
  int nhash, nmax, nn, ng;
  std::vector<int> htable, next, garbg;
  std::vector<Ullong> thehash;
  hfnT hash;
  Hashtable(int nh, int nv);
  int iget(const keyT &key);
  int iset(const keyT &key);
  int ierase(const keyT &key);
};

template<class keyT, class hfnT>
Hashtable<keyT,hfnT>::Hashtable(int nh, int nv):
  hash(sizeof(keyT)), nhash(nh), nmax(nv), nn(0), ng(0),
  htable(nh), next(nv), garbg(nv), thehash(nv) {
  for (int j=0; j<nh; j++) { htable[j] = -1; }
}

template<class keyT, class hfnT>
int Hashtable<keyT,hfnT>::iget(const keyT &key) {
  int j,k;
  Ullong pp = hash.fn(&key);
  j = (int)(pp % nhash);
  for (k = htable[j]; k != -1; k = next[k]) {
    if (thehash[k] == pp) {
      return k;
    }
  }
  return -1;
}

template<class keyT, class hfnT>
int Hashtable<keyT,hfnT>::iset(const keyT &key) {
  int j,k,kprev;
  Ullong pp = hash.fn(&key);
  j = (int)(pp % nhash);
  if (htable[j] == -1) {
    k = ng ? garbg[--ng] : nn++ ;
    htable[j] = k;
  } else {
    for (k = htable[j]; k != -1; k = next[k]) {
      if (thehash[k] == pp) {
        return k;
      }
      kprev = k;
    }
    k = ng ? garbg[--ng] : nn++ ;
    next[kprev] = k;
  }
  if (k >= nmax) NRthrow("storing too many values");
  thehash[k] = pp;
  next[k] = -1;
  return k;
}

template<class keyT, class hfnT>
int Hashtable<keyT,hfnT>::ierase(const keyT &key) {
  int j,k,kprev;
  Ullong pp = hash.fn(&key);
  j = (int)(pp % nhash);
  if (htable[j] == -1) return -1;
  kprev = -1;
  for (k = htable[j]; k != -1; k = next[k]) {
    if (thehash[k] == pp) {
      if (kprev == -1) htable[j] = next[k];
      else next[kprev] = next[k];
      garbg[ng++] = k;
      return k;
    }
    kprev = k;
  }
  return -1;
}

template<class keyT, class elT, class hfnT>
struct Hash : Hashtable<keyT, hfnT> {
  using Hashtable<keyT,hfnT>::iget;
  using Hashtable<keyT,hfnT>::iset;
  using Hashtable<keyT,hfnT>::ierase;
  std::vector<elT> els;

  Hash(int nh, int nm) : Hashtable<keyT, hfnT>(nh, nm), els(nm) {}

  void set(const keyT &key, const elT &el) { els[iset(key)] = el; }

  int get(const keyT &key, elT &el) {
    int ll = iget(key);
    if (ll < 0) return 0;
    el = els[ll];
    return 1;
  }

  elT& operator[] (const keyT &key) {
    int ll = iget(key);
    if (ll < 0) {
      ll = iset(key);
      els[ll] = elT();
    }
    return els[ll];
  }

  int count(const keyT &key) {
    int ll = iget(key);
    return (ll < 0 ? 0 : 1);
  }

  int erase(const keyT &key) {
    return (ierase(key) < 0 ? 0 : 1);
  }
};

struct Triel {
//Structure for an element in a descendancy tree of triangles, each having at most three daughters
  Point2D *pts; //Pointers to the array of the points;
  int p[3]; //The triangle's three vertices, always in CCW order
  int d[3]; //Pointers for up to three daughters
  int stat; //Nonzero if this element is "live"
  void setme(int a, int b, int c, Point2D *ptss) {
  //Set the data in a Triel
    pts = ptss;
    p[0] = a; p[1] = b; p[2] = c;
    d[0] = d[1] = d[2] = -1; //The value -1 mean no daughters
    stat = 1;
  }
  int contains(Point2D point) {
    //Return 1 if point is in the triangle, 0 if on boundary, -1 if outside
    //(CCW triangle is assumed)
    double d;
    int i,j,ztest=0;
    for (i=0; i<3; i++) {
      j = (i+1) % 3;
      d = (pts[p[j]].x-pts[p[i]].x)*(point.y-pts[p[i]].y) -
          (pts[p[j]].y-pts[p[i]].y)*(point.x-pts[p[i]].x);
      if (d < 0.0) return -1;
      if (d == 0.0) ztest = 1;
    }
    return (ztest ? 0 : 1);
  }
};

double incircle(Point2D d, Point2D a, Point2D b, Point2D c) {
//Return positive, zero, or negative value if point is respectively inside,
//on, or outside the circle through points a,b, and c
  Circle cc = circumcircle(a,b,c);
  double radd = SQR(d.x-cc.center.x) + SQR(d.y-cc.center.y);
  return (SQR(cc.radius) - radd);
}

struct Nullhash {
//Null hash function. Use a key (assumed to be already hashed) at its own hash
  Nullhash(int nn) {}
  inline Ullong fn(const void *key) const { return *((Ullong*)key); }
};

struct Delaunay {
//Structure for constructing a Delaunay triangulation from a given set of points
  int npts,ntri,ntree,ntreemax,opt; //Number of points, triangles, elements in the Triel list, and maximum of same
  double delx,dely; //Size of the bounding box
  std::vector< Point2D > pts;
  std::vector<Triel> thelist; //The list of Triel elements
  Hash<Ullong,int,Nullhash> *linehash; //Create the hash memories with null hash function
  Hash<Ullong,int,Nullhash> *trihash;
  int *perm;
  Delaunay(std::vector<Point2D> &pvec, int options = 0);
  //Construct the Delaunay triangulation from a vector of points.
  //The variable options is used by some applications
  Ranhash hashfn; //The raw hash function
  //The next four functions are explained in detail below
  void insertapoint(int r);
  int whichcontainspt(const Point2D &p, int strict = 0);
  int storetriangle(int a, int b, int c);
  void erasetriangle(int a, int b, int c, int d0, int d1, int d2);
  static unsigned int jran; //Random number counter
  static const double fuzz, bigscale;
};
const double Delaunay::fuzz = 1.0e-6;
const double Delaunay::bigscale = 1000.0;
unsigned int Delaunay::jran = 14921620;

Delaunay::Delaunay(std::vector<Point2D> &pvec, int options) :
  npts(pvec.size()), ntri(0), ntree(0), ntreemax(10*npts+1000),
  opt(options), pts(npts+3), thelist(ntreemax) {
//Construct Delaunay triangulation from a vector of points pvec.
//If bit 0 in options is nonzero, hash memories used in construction are deleted.
//(Some applications may want to use them and will set options to 1)
  int j;
  double xl,xh,yl,yh;
  linehash = new Hash<Ullong,int,Nullhash>(6*npts+12,6*npts+12);
  trihash = new Hash<Ullong,int,Nullhash>(2*npts+6,2*npts+6);
  perm = new int[npts]; //Permutation for randomizing point order
  xl = xh = pvec[0].x; //Copy points to local store and calculate
  yl = yh = pvec[0].y; //their bounding box
  for (j=0; j<npts; j++) {
    pts[j] = pvec[j];
    perm[j] = j;
    if (pvec[j].x < xl) xl = pvec[j].x;
    if (pvec[j].x > xh) xh = pvec[j].x;
    if (pvec[j].y < yl) yl = pvec[j].y;
    if (pvec[j].y > yh) yh = pvec[j].y;
  }
  delx = xh - xl; //Store bounding box dimensions, then construct
  dely = yh - yl; //the three fictious points and store them
  pts[npts] = Point2D(0.5*(xl+xh), yh+bigscale*dely);
  pts[npts+1] = Point2D(xl-0.5*bigscale*delx, yl-0.5*bigscale*dely);
  pts[npts+2] = Point2D(xl+0.5*bigscale*delx, yl-0.5*bigscale*dely);
  storetriangle(npts,npts+1,npts+2);
  //Create random permutation
  for (j=npts; j>0; j--) SWAP(perm[j-1],perm[hashfn.int64(jran++) % j]);
  for (j=0; j<npts; j++) insertapoint(perm[j]); //All the action is here!
  for (j=0; j<ntree; j++) { //Delete the huge root triangle and
    if (thelist[j].stat > 0) { //all of its connecting edges
      if (thelist[j].p[0] >= npts || thelist[j].p[1] >= npts || thelist[j].p[2] >= npts) {
        thelist[j].stat = -1;
        ntri--;
      }
    }
  }
  if (!(opt & 1)) { //Clean up, unless option bit says to do
    delete[] perm;
    delete trihash;
    delete linehash;
  }
}

void Delaunay::insertapoint(int r) {
//Add the point with index r incrementally to the Delaunay triangulation
  int i,j,k,l,s,tno,ntask,d0,d1,d2;
  Ullong key;
  int tasks[50], taski[50], taskj[50]; //Stacks (3 vertices) for legalizing edges
  for (j=0; j<3; j++) { //Find triangle containing point
    tno = whichcontainspt(pts[r],1); //Fuzz if it lies on an edge
    if (tno >= 0) break; //The desired result: Point is OK
    pts[r].x += fuzz*delx*(hashfn.real(jran++)-0.5);
    pts[r].y += fuzz*dely*(hashfn.real(jran++)-0.5);
  }
  if (j == 3) NRthrow("points degenerate even after fuzzing");
  ntask = 0;
  i = thelist[tno].p[0]; j = thelist[tno].p[1]; k = thelist[tno].p[2];
  //The following line is relevant only when the indicated bit in opt is set.
  //This feature is used by the convex hull application and causes any points
  // already known to be interior to the convex hull to be omitted from the
  // triangulation, saving time (but giving in an incomplete triangulation).
  if (opt & 2 && i < npts && j < npts && k < npts) return;
  d0 = storetriangle(r,i,j); //Create three triangles and queue them for legal edge tests
  tasks[++ntask] = r; taski[ntask] = i; taskj[ntask] = j;
  d1 = storetriangle(r,j,k);
  tasks[++ntask] = r; taski[ntask] = j; taskj[ntask] = k;
  d2 = storetriangle(r,k,i);
  tasks[++ntask] = r; taski[ntask] = k; taskj[ntask] = i;
  erasetriangle(i,j,k,d0,d1,d2); //Erase the old triangle
  while (ntask) { //Legalize edges recursively
    s = tasks[ntask]; i = taski[ntask]; j = taskj[ntask--];
    key = hashfn.int64(j) - hashfn.int64(i); //Look up fourth point
    if (!linehash->get(key,l)) continue; //Case of no triangle on other side
    if (incircle(pts[l],pts[j],pts[s],pts[i]) > 0.0) { //Needs legalizing?
      d0 = storetriangle(s,l,j); //Create two new triangles
      d1 = storetriangle(s,i,l);
      erasetriangle(s,i,j,d0,d1,-1); //and erase old ones
      erasetriangle(l,j,i,d0,d1,-1);
      key = hashfn.int64(i) - hashfn.int64(j); //Erase line in both directions
      linehash->erase(key);
      key = 0 - key; //Unsigned, hence binary minus
      linehash->erase(key);
      //Two new edges now need checking:
      tasks[++ntask] = s; taski[ntask] = l; taskj[ntask] = j;
      tasks[++ntask] = s; taski[ntask] = i; taskj[ntask] = l;
    }
  }
}

int Delaunay::whichcontainspt(const Point2D &p, int strict) {
//Given point p, return index in thelist of the triangle in the
//triangulation that contains it, or return -1 for failure.
//If strict is nonzero, require strict containment, otherwise
//allow the point to lie on an edge
  int i,j,k=0;
  while (thelist[k].stat <= 0) { //Descend in tree until reach a "live" triangle
    for (i=0; i<3; i++) { //Check up to three daughters
      if ((j = thelist[k].d[i]) < 0) continue; //Daughter doesn't exist
      if (strict) {
        if (thelist[j].contains(p) > 0) break;
      } else { //Yes, descend on this branch
        if (thelist[j].contains(p) >= 0) break;
      }
    }
    if (i == 3) return -1; //No daughters contain the point
    k = j; //Set new mother
  }
  return k; //Normal return
}

void Delaunay::erasetriangle(int a, int b, int c, int d0, int d1, int d2) {
//Erase triangle abc in trihash and inactivate it in thelist after setting its daughters
  Ullong key;
  int j;
  key = hashfn.int64(a) ^ hashfn.int64(b) ^ hashfn.int64(c);
  if (trihash->get(key,j) == 0) NRthrow("nonexistent triangle");
  trihash->erase(key);
  thelist[j].d[0] = d0; thelist[j].d[1] = d1; thelist[j].d[2] = d2;
  thelist[j].stat = 0;
  ntri--;
}

int Delaunay::storetriangle(int a, int b, int c) {
//Store a triangle with vertices a, b, c in trihash.
//Store its points in linehash under keys to opposite sides.
//Add it to thelist, returning its index there.
  Ullong key;
  thelist[ntree].setme(a,b,c,&pts[0]);
  key = hashfn.int64(a) ^ hashfn.int64(b) ^ hashfn.int64(c);
  trihash->set(key,ntree);
  key = hashfn.int64(b) - hashfn.int64(c);
  linehash->set(key,a);
  key = hashfn.int64(c) - hashfn.int64(a);
  linehash->set(key,b);
  key = hashfn.int64(a) - hashfn.int64(b);
  linehash->set(key,c);
  if (++ntree == ntreemax) NRthrow("thelist is sized too small");
  ntri++;
  return (ntree-1);
}
