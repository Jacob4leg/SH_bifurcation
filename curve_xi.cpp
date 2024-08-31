#include "curve_xi.h"

using namespace std;
using namespace capd;

/* function for finding approximate solution of G(\xi,x) = 0, where x is treated as the parameter */
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


/* Function, which uses second derivative, to prove the existence of unique maximum of \xi(x). */
bool proveMaximum(interval X, interval XI, IC2PoincareMap &pm, interval S=interval(-1,1)*1e-10) {
    static const Interval sqrt2 = Interval(sqrt(2.0));

    Interval x0 = X.mid();
    Interval xi0 = XI.mid();

    IMatrix DPhi(5,5), D(5,5); // monodromy matrix and derivative
    IHessian HPhi(5,5), H(5,5); // Hessian

    C1Rect2Set s({x0,0,(sqr(x0) - 1)/sqrt2, 0, xi0});
    IVector u = pm(s,DPhi,2);
    D = pm.computeDP(u,DPhi);

    IVector Fx({u[3],D[3][0] + D[3][2] * sqrt2 * x0});

    DPhi.clear();
    D.clear();

    C2Rect2Set set({X,0,(sqr(X) - 1)/sqrt2,0,XI});
    IVector U = pm(set,DPhi,HPhi,2);
    pm.computeDP(U,DPhi,HPhi,D,H);

    // below is the collection of all needed derivatives
    Interval G_xi = D[3][4];
    Interval G_x = D[3][0] + D[3][2] * sqrt2 * X;
    Interval xi_x = - G_x / G_xi;
    Interval G_x_x = H(3,0,0) * 2 + sqrt2 * (H(3,0,2) * 2 + D[3][2]);
    Interval G_xi_xi = H(3,4,4) * 2;
    Interval G_x_xi = H(3,0,4) * 2;
    Interval xi_x_x = -(G_x_x + (2 * G_x_xi + G_xi_xi * xi_x) * xi_x) / G_xi;

    IMatrix DF({{G_x,G_xi},{G_x_x,G_x_xi}}); // derivative of the function H(\xi,x) = (G(\xi,x),G_x(\xi,x))

    IVector N = - matrixAlgorithms::gauss(DF,Fx);

    cout << "N=" << N << ", S=" << IVector{S,S} << endl;
    cout << "xi''(x)=" << xi_x_x << endl;
    cout << "(xi^*, x^*) belongs to the cartesian product of " << IVector{xi0,x0} - N << endl;
 
    return subset(N,IVector{S,S}) and (xi_x_x < 0);    
}

/*
Function which uses the theorem for Interval Newton Operator to prove the existence of the smooth curve within small interval \Xi = [\xi-eps,\xi+eps].
Each point of the curve corresponds to the unique periodic orbit of the Swift-Hohenberg equation. Here x is treated as the parameter
If the function returns false, then we shrink x range and try again.
*/
tuple<bool,interval,interval> proveOrbit_xi(double xi, interval x, IC2PoincareMap &pm, interval X_end) {
    static const interval sqrt2=interval(sqrt(2.0));
    static interval S = interval(-1,1)*8e-6;
    static interval maxS = interval(-1,1)*5e-4;

    interval xi0 = xi;
    
    C0HORect2Set s({x,0,(sqr(x)-1)/sqrt2,0,xi0});
    interval XI = xi0 + S;
    interval y, N, xi_x_x;
    IVector v1, v2;

    try{
        y = pm(s,2)[3];
        
        IMatrix DPhi(5,5), D(5,5);
        IHessian HPhi(5,5), H(5,5);
        C2Rect2Set set({x,0,(sqr(x) - 1)/sqrt2,0,XI});
        v1 = pm(set,DPhi,HPhi);
        v2 = pm(set,DPhi,HPhi);
        pm.computeDP(v2,DPhi,HPhi,D,H); // computation of all derivatives of the order <= 2

        // all needed derivatives
        Interval fxi = v2[3];
        Interval G_xi = D[3][4];
        N = y / G_xi;
        Interval G_x = D[3][0] + D[3][2] * sqrt2 * x;
        Interval xi_x = - G_x / G_xi;
        Interval G_x_x = H(3,0,0) * 2 + sqrt2 * (H(3,0,2) * 2 + D[3][2]);
        Interval G_xi_xi = H(3,4,4) * 2;
        Interval G_x_xi = H(3,0,4) * 2;
        xi_x_x = -(G_x_x + (2 * G_x_xi + G_xi_xi * xi_x) * xi_x) / G_xi;

        // box, which contains maximum of the function \xi(x)
        interval X_max = -1.5824941113082425 + interval(-1,1)*1e-10;
        interval XI_max = 2.0316516135713902 + interval(-1,1)*1e-10;

        // we verify if max box is contained in the enclosure of the curve
        if(x.contains(X_max) && XI.contains(XI_max)) {
            cout << "checking maximum" << endl;
            // we need to decrease the tolerance of ODE solver for max checking purpose
            pm.getSolver().setAbsoluteTolerance(2e-13);
            pm.getSolver().setRelativeTolerance(2e-13);
            bool is_maximum = proveMaximum(X_max,XI_max,pm); // checking if this box contains maximum
            if(is_maximum)
                cout << "maximum in XI=" << XI << ", x=" << x << endl;
            else
                cout << "Problems in determining if is max" << endl;
            pm.getSolver().setAbsoluteTolerance(2e-7);
            pm.getSolver().setRelativeTolerance(2e-7);
        }
    }catch(exception e) {
        S *= 0.95; // we shrink \xi range if we failed
        return {false,0,0};
    }

    bool geometric_ok = x<-1 and v1[0]>1 and v2[0]<1 and v2[0]>-1; // geometric condition
    return {subset(N,S) and geometric_ok and (xi_x_x < 0) and intersection(maxS,N*interval(-1,1)*1.001,S),XI,xi_x_x};
}

/*
Function determines the existence of the curve in range x \in X_*.
Returns number of subintervals
*/
tuple<int,interval> proveCurve_xi(double xi, interval X_start, interval X_end) {
    int counter_subintervals = 0;
    // definition of double Poincare map
    DMap vf("var:x,y,z,w,xi;fun:y,z,w,x*(1-x^2)-xi*z,0;");
    DOdeSolver solver(vf,20);
    DCoordinateSection section(5,1);
    DPoincareMap pm(solver,section);

    // definition of interval Poincare map
    IMap ivf("var:x,y,z,w,xi;fun:y,z,w,x*(1-x^2)-xi*z,0;");
    IC2OdeSolver isolver(ivf,8);
    ICoordinateSection isection(5,1);
    IC2PoincareMap ipm(isolver,isection);
    isolver.setAbsoluteTolerance(2e-7);
    isolver.setRelativeTolerance(2e-7);

    interval xi_threshold = xi;
    interval set0 = xi + interval(-1,1);
    interval xi_x_x(-50000); // for storing whole enclosure of second derivative. -50000 is chosen, because first subinterval contains this value
    interval X,r;
    double x = X_start.leftBound();

    double L=x;
    double delta = 1e-8;

    bool first_step = true;

    while(x <= X_end) {
        x = L + 0.5 * delta;
        xi = findOrbit_xi(xi,x,pm); // we find an approximate zero of G

        bool last_step = (L + delta > X_end);
        double R = capd::min(L+delta, X_end.rightBound());
        interval X(L,R);
        auto [proven_orbit,set,second_der] = proveOrbit_xi(xi,X,ipm,X_end); // we use Interval Newton Method to determine the existence of the smooth curve within some box

        if(!proven_orbit) {
            delta *= 0.95; // if we failed, then we shrink x range 
        }
        else {
            if(!intersection(set0,set,r)) { // checking if curves on neighbouring intervals glue correctly into one curve
                cout << "curve xi not glued correctly" << endl;
                cout << "set0=" << set0 << endl;
                cout << "set=" << set << endl;
            }
            if(counter_subintervals++%100 == 0 || X.contains(X_end)) 
                cout << "XI=" << set << ", X=" << X << ", xi''(x) in " << second_der << endl;

            xi_x_x = intervalHull(second_der,xi_x_x); // updating the enclosure of the second derivative
            L = L+delta;
            delta*=1.001;
            set0 = set;
        }

        if(first_step) { // here we check if first threshold box is contained in the first subinterval of the enclosure of the curve
            bool contains_first_threshold_box = X.contains(X_start) and set.contains(xi_threshold);
            if(contains_first_threshold_box)
                cout << "threshold box in parameterization" << endl;
            else
                cout << "threshold box not in parameterization" << endl;
            first_step = false;
        }
        if(last_step) { // here we check if second threshold box is contained in the last subinterval of the enclosure of the curve
            bool contains_second_threshold_box = X.contains(X_end) and set.contains(xi_threshold);
            if(contains_second_threshold_box)
                cout << "threshold box in parameterization" << endl;
            else
                cout << "threshold box not in parameterization" << endl;
            break;
        }
    }
    cout << "Whole curve computed" << endl;
    return {counter_subintervals,xi_x_x};
}