#include<iostream>
#include<chrono>
#include "capd/capdlib.h"
#include "curve.h"

using namespace std;
using namespace chrono;
using namespace capd;

double findOrbit_x(double xi, double x, DPoincareMap& pm){
  pm.getVectorField().setParameter(0,xi);
  DMatrix D{{1,0,0,0},{0,0,0,0},{sqrt(2.)*x,0,0,0},{0,0,0,0}}, T(4,4);
  DVector u{x,0,(x*x-1)/sqrt(2.),0};
  for(int i=0;i<2;++i){
    u = pm(u,T);
    D = T*D;
  }
  D = pm.computeDP(u,D);
  return x - u[3]/D[3][0];
}

double findOrbit_xi(double xi, double x, DPoincareMap& pm) {
    DMatrix D({{1,0,0,0,0},{0,0,0,0,0},{sqrt(2.0)*x,0,0,0,0},{0,0,0,0,0},{0,0,0,0,1}}), T(5,5);
    DVector u({x,0,(x*x-1)/sqrt(2.0),0,xi});
    for(int i=0; i<2;++i) {
        u = pm(u,T);
        D = T * D;
    }
    D = pm.computeDP(u,D);
    return xi - u[3] / D[3][4];
}

int counter = 0;
tuple<bool,interval> proveOrbit_x(interval xi, double x, IPoincareMap& pm){
    static const interval sqrt2=interval(sqrt(2.0));
    static interval S = interval(-1,1)*8e-6;
    static interval maxS = interval(-1,1)*5e-4;

    pm.getVectorField().setParameter(0,xi);
    interval x0 = x;
    C0HORect2Set s({x0,0,(sqr(x0)-1)/sqrt2,0});
    interval y;

    try{
        y = pm(s,2)[3];
    }
    catch(poincare::PoincareException<C0HORect2Set>) {
        return {false,0};
    }

    interval X = x0 + S;
    IMatrix m(4,4),D(4,4);
    C1HORect2Set set({X,0,(sqr(X)-1)/sqrt2,0});
    IVector v1,v2;

    try{
        v1 = pm(set);
        v2 = pm(set,m);
    }
    catch(poincare::PoincareException<C1HORect2Set>) {
        S *= 0.95;
        return {false,0};
    }
    D = pm.computeDP(v2,m);
    interval derivative = D[3][0] + D[3][2] * sqrt2 * X;
    if(derivative.contains(0)) {
        return {false,0};
    }
        
    interval N = -y/derivative;
    interval end = interval(266291)/131072;
    
    bool is_subset = subset(N,S);
    bool geometric_ok = X<-1 and v1[0]>1 and v2[0]<1 and v2[0]>-1;
    
    if(is_subset and geometric_ok)
        if(counter++% 1000==0 || xi.contains(end)) cout << "xi=" << xi << ", X=" << X << ", N=" << N << ", S=" << S << endl;

    return {is_subset and geometric_ok and intersection(maxS,N*interval(-1,1)*1.03,S),X};
}

bool proveMaximum(interval X, interval XI, IC2PoincareMap &pm, interval S=interval(-1,1)*1e-10) {
    static const Interval sqrt2 = Interval(sqrt(2.0));

    Interval x0 = X.mid();
    Interval xi0 = XI.mid();

    IMatrix DPhi(5,5), D(5,5);
    IHessian HPhi(5,5), H(5,5);


    C1Rect2Set s({x0,0,(sqr(x0) - 1)/sqrt2, 0, xi0});
    IVector u = pm(s,DPhi,2);
    D = pm.computeDP(u,DPhi);

    IVector Fx({u[3],D[3][0] + D[3][2] * sqrt2 * x0});

    DPhi.clear();
    D.clear();

    C2Rect2Set set({X,0,(sqr(X) - 1)/sqrt2,0,XI});
    IVector U = pm(set,DPhi,HPhi,2);
    pm.computeDP(U,DPhi,HPhi,D,H);

    Interval G_xi = D[3][4];

    Interval G_x = D[3][0] + D[3][2] * sqrt2 * X;

    Interval xi_x = - G_x / G_xi;
    Interval G_x_x = H(3,0,0) * 2 + sqrt2 * (H(3,0,2) * 2 + D[3][2]);
    Interval G_xi_xi = H(3,4,4) * 2;
    Interval G_x_xi = H(3,0,4) * 2;
    Interval xi_x_x = -(G_x_x + (2 * G_x_xi + G_xi_xi * xi_x) * xi_x) / G_xi;

    IMatrix DF({{G_x,G_xi},{G_x_x,G_x_xi}});

    IVector N = - matrixAlgorithms::gauss(DF,Fx);

    cout << "N=" << N << ", S=" << IVector{S,S} << endl;
    cout << "xi''(x)=" << xi_x_x << endl;
 
    return subset(N,IVector{S,S}) and (xi_x_x < 0);    
}

int counter_xi = 0;
tuple<bool,interval,interval> proveOrbit_xi(double xi, interval x, IC2PoincareMap &pm, interval X_end) {
    static const interval sqrt2=interval(sqrt(2.0));
    static interval S = interval(-1,1)*8e-6;
    static interval maxS = interval(-1,1)*5e-4;

    interval xi0 = xi;
    C0HORect2Set s({x,0,(sqr(x)-1)/sqrt2,0,xi0});
    interval y;

    try{
        y = pm(s,2)[3];
    }
    catch(poincare::PoincareException<C0HORect2Set>) {
        return {false,0,0};
    }

    interval XI = xi0 + S;
    IMatrix DPhi(5,5), D(5,5);
    IHessian HPhi(5,5), H(5,5);
    C2Rect2Set set({x,0,(sqr(x) - 1)/sqrt2,0,XI});
    IVector v1,v2;

    try{
        v1 = pm(set,DPhi,HPhi);
        v2 = pm(set,DPhi,HPhi);
    }
    catch(poincare::PoincareException<C2Rect2Set>) {
        S *= 0.95;
        return {false,0,0};
    }
    pm.computeDP(v2,DPhi,HPhi,D,H);

    Interval fxi = v2[3];
    Interval G_xi = D[3][4];

    Interval N = y / G_xi;

    Interval G_x = D[3][0] + D[3][2] * sqrt2 * x;

    Interval xi_x = - G_x / G_xi;

    Interval G_x_x = H(3,0,0) * 2 + sqrt2 * (H(3,0,2) * 2 + D[3][2]);
    Interval G_xi_xi = H(3,4,4) * 2;
    Interval G_x_xi = H(3,0,4) * 2;
    Interval xi_x_x = -(G_x_x + (2 * G_x_xi + G_xi_xi * xi_x) * xi_x) / G_xi;
    interval X_max = -1.5824941113082425 + interval(-1,1)*1e-10;
    interval XI_max = 2.0316516135713902 + interval(-1,1)*1e-10;

    if(x.contains(X_max) && XI.contains(XI_max)) {
        cout << "checking maximum" << endl;
        pm.getSolver().setAbsoluteTolerance(2e-13);
        pm.getSolver().setRelativeTolerance(2e-13);
        bool is_maximum = proveMaximum(X_max,XI_max,pm);
        if(is_maximum)
            cout << "maximum in XI=" << XI << ", x=" << x << endl;
        else
            cout << "Problems in determining if is max" << endl;
        pm.getSolver().setAbsoluteTolerance(2e-7);
        pm.getSolver().setRelativeTolerance(2e-7);
    }

    bool is_subset = subset(N,S);
    bool geometric_ok = x<-1 and v1[0]>1 and v2[0]<1 and v2[0]>-1;

    if(is_subset and geometric_ok)
        if(counter_xi++%100 == 0 || x.contains(X_end)) cout << "XI=" << XI << ", x=" << x << ", N=" << N << ", S=" << S << ", xi''(x)=" << xi_x_x << endl;
    

    return {is_subset and geometric_ok and (xi_x_x < 0) and intersection(maxS,N*interval(-1,1)*1.001,S),XI,xi_x_x};
}


void prove_curve_x(double x, Interval end_X_to_check) {
    cout.precision(14);

    DMap vf("par:xi;var:x,y,z,w;fun:y,z,w,x*(1-x^2)-xi*z;");
    DOdeSolver solver(vf,20);
    DCoordinateSection section(4,1);
    DPoincareMap pm(solver,section);

    IMap ivf("par:xi;var:x,y,z,w;fun:y,z,w,x*(1-x^2)-xi*z;");
    IOdeSolver isolver(ivf,8);
    ICoordinateSection isection(4,1);
    IPoincareMap ipm(isolver,isection);
    isolver.setAbsoluteTolerance(2e-7);
    isolver.setRelativeTolerance(2e-7);

    long double L=0, xi=0;
    long double delta = 1e-5;

    interval set0 = x + interval(-1,1);
    interval r;
    interval end = interval(266291)/131072;

    while(L<=end){

        double xi = L+0.5*delta;
        x = findOrbit_x(xi,x,pm);

        if(xi >= 2.031){
            isolver.setAbsoluteTolerance(2e-8);
            isolver.setRelativeTolerance(2e-8);
        }

        if(L + delta > end) {
            cout << "last step" << endl;
            interval end_interval = interval(L,end.rightBound());
            auto res_orbit_check = proveOrbit_x(end_interval,x,ipm);
            bool proven_orbit = get<0>(res_orbit_check);
            interval set = get<1>(res_orbit_check);

            if(!proven_orbit) {
                cout << "Problems to prove x orbit in last step" << endl;
            }

            if(!intersection(set0,set,r)) {
                cout << "x curve not glued correctly" << endl;
            }

            bool contains_threshold_box = end_interval.contains(end) and set.contains(end_X_to_check);

            if(!contains_threshold_box) {
                cout << "threshold box not in parameterization" << endl;
            }
            break;
        }
        auto res_orbit_check = proveOrbit_x(interval(L,L+delta),x,ipm);
        bool proven_orbit = get<0>(res_orbit_check);
        interval set = get<1>(res_orbit_check);

        if(!(proven_orbit))  
            delta*=0.95;
        else {
            if(!intersection(set0,set,r)) {
                cout << "x curve not glued correctly" << endl;
                cout << "set0=" << set0 << endl;
                cout << "set=" << set << endl;
                break;
            }
            L = L+delta;
            delta*=1.001;
            set0 = set;
        }
    }
    cout << "subdivision:" << counter << endl;
}


void prove_curve_xi(double xi, interval X_start, interval X_end) {
    DMap vf("var:x,y,z,w,xi;fun:y,z,w,x*(1-x^2)-xi*z,0;");
    DOdeSolver solver(vf,20);
    DCoordinateSection section(5,1);
    DPoincareMap pm(solver,section);

    IMap ivf("var:x,y,z,w,xi;fun:y,z,w,x*(1-x^2)-xi*z,0;");
    IC2OdeSolver isolver(ivf,8);
    ICoordinateSection isection(5,1);
    IC2PoincareMap ipm(isolver,isection);
    isolver.setAbsoluteTolerance(2e-7);
    isolver.setRelativeTolerance(2e-7);

    interval xi_threshold = interval(266291)/131072;
    interval set0 = xi + interval(-1,1);
    interval xi_x_x(-50000);
    interval X,r;
    double x = X_start.leftBound();

    long double L=x;
    long double delta = 1e-8;

    bool first_iter = true;

    while(x <= X_end) {
        x = L + 0.5 * delta;
        xi = findOrbit_xi(xi,x,pm);

        if(x + delta > X_end) {
            cout << "last step" << endl;
            X = interval(L,X_end.rightBound());
            auto res_orbit_check = proveOrbit_xi(xi,X,ipm,X_end);
            bool proven_orbit = get<0>(res_orbit_check);
            interval set = get<1>(res_orbit_check);
            interval second_der = get<2>(res_orbit_check);

            

            bool contains_second_threshold_box = X.contains(X_end) and set.contains(xi_threshold);
            if(!proven_orbit)
                cout << "problems to prove xi curve in last step" << endl;
            if(!intersection(set0, set, r))
                cout << "xi curve not glued correctly" << endl;
            if(contains_second_threshold_box) 
                cout << "threshold box in parameterization" << endl;
            else
                cout << "threshold box not in parameterization" << endl;

            

            if(second_der.left() < xi_x_x.left())
                xi_x_x.setLeftBound(second_der.leftBound());
            if(second_der.right() > xi_x_x.right())
                xi_x_x.setRightBound(second_der.rightBound());

            cout << "enclosure of second derivative=" << xi_x_x << endl;
            break;
        }

        X = interval(L,L+delta);
        auto res_orbit_check = proveOrbit_xi(xi,interval(L,L+delta),ipm,X_end);
        bool proven_orbit = get<0>(res_orbit_check);
        interval set = get<1>(res_orbit_check);
        interval second_der = get<2>(res_orbit_check);
        

        if(first_iter) {
            bool contains_first_threshold_box = X.contains(X_start) and set.contains(xi_threshold);
            if(contains_first_threshold_box)
                cout << "threshold box in parameterization" << endl;
            else
                cout << "threshold box not in parameterization" << endl;
            first_iter = false;
        }

        if(!proven_orbit) {
            delta *= 0.95;
        }
        else {
            if(!intersection(set0,set,r)) {
                cout << "curve xi not glued correctly" << endl;
                cout << "set0=" << set0 << endl;
                cout << "set=" << set << endl;
            }
            if(second_der.left() < xi_x_x.left())
                xi_x_x.setLeftBound(second_der.leftBound());
            if(second_der.right() > xi_x_x.right())
                xi_x_x.setRightBound(second_der.rightBound());

            L = L+delta;
            delta*=1.001;
            set0 = set;
        }

    }
    cout << "subdivision: " << counter_xi << endl;
}

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

    double x1 = -1.0849406631521703;
    double x2 = -1.0845890871772264;

    interval X1(-1.5825372476037, -1.5825328805433);
    interval X2(-1.5824461925798, -1.5824418541956);


    double xi_threshold = 266291./131072;

    interval X1_end = find_tight_enclosure(X1);
    interval X2_end = find_tight_enclosure(X2);

    cout << "X1_end=" << X1_end << ", diam=" << diam(X1_end) << endl;
    cout << "X2_end=" << X2_end << ", diam=" << diam(X2_end) << endl;
    

    int subdivisions_first_curve = 0;

    // first curve
    prove_curve_x(x1,X1_end);
    auto after_first_curve = high_resolution_clock::now();
    // // second curve
    subdivisions_first_curve = counter;
    counter = 0;
    prove_curve_x(x2,X2_end);
    auto after_second_curve = high_resolution_clock::now();
    // xi curve
    prove_curve_xi(xi_threshold, X1_end, X2_end);
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
    cout << "Subdivisions second curve = " << counter << endl;
    cout << "Subdivisions xi curve = " << counter_xi << endl;


    return 0;
}
