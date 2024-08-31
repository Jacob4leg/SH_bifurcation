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

    IMatrix D{{1,0,0,0},{0,0,0,0},{sqrt2*X,0,0,0},{0,0,0,0}};
    C1HORect2Set set(C1Rect2Set::C0BaseSet({X,0,(sqr(X)-1)/sqrt2,0}),C1Rect2Set::C1BaseSet{D});
    IVector v = pm(set,D,2);
    interval derivative = pm.computeDP(v,D)[3][0];
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
interval find_tight_enclosure(interval X, interval xi) {
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

    // threshold xi
    double xi_threshold = 266291./131072;
    
    // first curve
    auto [subdivisions_first_curve, X1_threshold] = proveCurve_x(x1,xi_threshold);
    X1_threshold = find_tight_enclosure(X1_threshold,xi_threshold);
    cout << "X1_threshold=" << X1_threshold << endl;
    auto after_first_curve = high_resolution_clock::now();

    // second curve
    auto [subdivisions_second_curve, X2_threshold] = proveCurve_x(x2,xi_threshold);
    X2_threshold = find_tight_enclosure(X2_threshold,xi_threshold);
    cout << "X2_threshold=" << X2_threshold << endl;
    auto after_second_curve = high_resolution_clock::now();

    // xi curve
    auto [subdivisions_xi_curve, second_der] = proveCurve_xi(xi_threshold, X1_threshold, X2_threshold);
    auto stop = high_resolution_clock::now();

    auto duration_first_curve = duration_cast<milliseconds>(after_first_curve - start);
    auto duration_second_curve = duration_cast<milliseconds>(after_second_curve - after_first_curve);
    auto duration_xi_curve = duration_cast<milliseconds>(stop - after_second_curve);
    auto duration_all = duration_cast<milliseconds>(stop - start);

    cout.precision(2);
    cout << "Time for first curve: " << double(duration_first_curve.count()) / 1000  << "s" << endl;
    cout << "Time for second curve: " << double(duration_second_curve.count()) / 1000  << "s" << endl;
    cout << "Time for xi curve: " << double(duration_xi_curve.count()) / 1000  << "s" << endl;
    cout << "Total time elapsed: " << double(duration_all.count()) / 1000  << "s" << endl;
    cout << endl;
    cout.precision(17);
    cout << "Subdivisions first curve = " << subdivisions_first_curve << endl;
    cout << "Subdivisions second curve = " << subdivisions_second_curve << endl;
    cout << "Subdivisions xi curve = " << subdivisions_xi_curve << endl;
    cout << "X1_threshold = " << X1_threshold << endl;
    cout << "X2_threshold = " << X2_threshold << endl; 
    cout << "second derivative xi''(x) = " << second_der << endl;

    return 0;
}
