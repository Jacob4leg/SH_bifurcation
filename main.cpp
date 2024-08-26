#include<iostream>
#include<chrono>
#include "capd/capdlib.h"
#include "curve_x.h"
#include "curve_xi.h"

using namespace std;
using namespace chrono;
using namespace capd;

/* Function, which seeks for tightiest enclosure for the zero  */
Interval interval_newton(Interval X, Interval xi, IPoincareMap &pm, int iterations=10) {
    static const interval sqrt2=interval(sqrt(2.0));
    pm.getVectorField().setParameter(0,xi);
    Interval x0 = X.mid();

    C0HORect2Set s({x0,0,(sqr(x0)-1)/sqrt2,0});
    interval y = pm(s,2)[3];

    IMatrix m(4,4),D(4,4);
    C1HORect2Set set({X,0,(sqr(X)-1)/sqrt2,0});
    IVector v = pm(set,m,2);
    D = pm.computeDP(v,m);
    interval derivative = D[3][0] + D[3][2] * sqrt2 * X;
    if(derivative.contains(0))
        return 0;
    interval N = x0-y/derivative;

    interval r;
    if(!intersection(X,N,r))
        return 0;
    if(iterations > 0)
        return interval_newton(r,xi,pm,iterations-1);
    return r;
}

/* The same but for threshold parameter and specified Poincare map */
interval find_tight_enclosure(interval X) {
    interval xi = interval(266291)/131072;
    IMap ivf("par:xi;var:x,y,z,w;fun:y,z,w,x*(1-x^2)-xi*z;");
    IOdeSolver isolver(ivf,20);
    ICoordinateSection isection(4,1);
    IPoincareMap ipm(isolver,isection);
    return interval_newton(X,xi,ipm);
}

int main() {
    cout.precision(17);

    auto start = high_resolution_clock::now();

    // initial values for xi = 0
    double x1 = -1.0849406631521703;
    double x2 = -1.0845890871772264;

    // threshold boxes for xi = 266291./131072
    interval X1(-1.5825372476037, -1.5825328805433);
    interval X2(-1.5824461925798, -1.5824418541956);


    // threshold xi
    double xi_threshold = 266291./131072;

    // tight enclosures for threshold x
    interval X1_end = find_tight_enclosure(X1);
    interval X2_end = find_tight_enclosure(X2);

    cout << "X1_end=" << X1_end << ", diam=" << diam(X1_end) << endl;
    cout << "X2_end=" << X2_end << ", diam=" << diam(X2_end) << endl;
    
    // first curve
    int subdivisions_first_curve = proveCurve_x(x1,X1_end);
    auto after_first_curve = high_resolution_clock::now();
    // // second curve
    int subdivisions_second_curve = proveCurve_x(x2,X2_end);
    auto after_second_curve = high_resolution_clock::now();
    // xi curve
    int subdivisions_xi_curve = proveCurve_xi(xi_threshold, X1_end, X2_end);
    auto stop = high_resolution_clock::now();

    auto duration_first_curve = duration_cast<milliseconds>(after_first_curve - start);
    auto duration_second_curve = duration_cast<milliseconds>(after_second_curve - after_first_curve);
    auto duration_xi_curve = duration_cast<milliseconds>(stop - after_second_curve);
    auto duration_all = duration_cast<milliseconds>(stop - start);
    cout << "Time for first curve: " << double(duration_first_curve.count()) / 1000  << "s" << endl;
    cout << "Time for second curve: " << double(duration_second_curve.count()) / 1000  << "s" << endl;
    cout << "Time for xi curve: " << double(duration_xi_curve.count()) / 1000  << "s" << endl;
    cout << "Total time elapsed: " << double(duration_all.count()) / 1000  << "s" << endl;
    cout << endl;
    cout << "Subdivisions first curve = " << subdivisions_first_curve << endl;
    cout << "Subdivisions second curve = " << subdivisions_second_curve << endl;
    cout << "Subdivisions xi curve = " << subdivisions_xi_curve << endl;

    return 0;
}
