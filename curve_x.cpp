#include "curve_x.h"

using namespace std;
using namespace capd;

/* function for finding an approximate solution of G(\xi,x) = 0, where \xi is treated as the parameter */
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

int counter = 0; // for counting the number of needed subintervals for enclosing the curves
/*
Function which uses the theorem for Interval Newton Operator to prove the existence of the smooth curve within small interval X = [x-eps,x+eps].
Each point of the curve corresponds to the unique periodic orbit of the Swift-Hohenberg equation. Here \xi is treated as the parameter.
If the function returns false, then we shrink \xi range and try again.
*/
tuple<bool,interval> proveOrbit_x(interval xi, double x, IPoincareMap& pm){
    static const interval sqrt2=interval(sqrt(2.0));
    static interval S = interval(-1,1)*8e-6; // interval [-eps,eps]
    static interval maxS = interval(-1,1)*5e-4;

    pm.getVectorField().setParameter(0,xi); 
    interval x0 = x;
    C0HORect2Set s({x0,0,(sqr(x0)-1)/sqrt2,0}); // initial vector as an input to the vector field
    interval y;

    try{ // we try to integrate vector field with given \xi range
        y = pm(s,2)[3];
    }
    catch(poincare::PoincareException<C0HORect2Set>) { // if finding the value of Poincare map failed, we shrink the range of \xi parameter
        return {false,0};
    }

    interval X = x0 + S; // interval [x-eps,x+eps]
    IMatrix D{{1,0,0,0},{0,0,0,0},{sqrt2*X,0,0,0},{0,0,0,0}};
    C1HORect2Set set(C1Rect2Set::C0BaseSet({X,0,(sqr(X)-1)/sqrt2,0}),C1Rect2Set::C1BaseSet{D});
    IVector v1,v2;

    try{ // again we try to integrate vector field, but here x is an interval
        v1 = pm(set);
        v2 = pm(set,D);
    }
    catch(poincare::PoincareException<C1HORect2Set>) {
        S *= 0.95; // we shrink [-eps,eps] interval if integration failed
        return {false,0};
    }
    interval derivative = pm.computeDP(v2,D)[3][0];
    if(derivative.contains(0)) { // if derivative contains zero, then we try with smaller \xi range
        return {false,0};
    }
        
    interval N = -y/derivative;
    bool is_subset = subset(N,S); // Interval Newton Operator condition
    bool geometric_ok = X<-1 and v1[0]>1 and v2[0]<1 and v2[0]>-1; // geometric conditions for periodic orbits
    
    interval end = interval(266291)/131072; // \xi threshold value
    if(is_subset and geometric_ok)
        if(counter++% 1000==0 || xi.contains(end)) cout << "xi=" << xi << ", X=" << X << ", N=" << N << ", S=" << S << endl;

    // returns the result and X range needed to verify if the curve is glued correctly
    return {is_subset and geometric_ok and intersection(maxS,N*interval(-1,1)*1.03,S),X};
}

/*
Function determines the existence of the curve in range \xi \in [0, \xi_*].
Returns number of subintervals and the last X^N set
*/
tuple<int,interval> proveCurve_x(double x, interval end) {  
    counter = 0;
    // definition of double Poincare map
    DMap vf("par:xi;var:x,y,z,w;fun:y,z,w,x*(1-x^2)-xi*z;");
    DOdeSolver solver(vf,20);
    DCoordinateSection section(4,1);
    DPoincareMap pm(solver,section);

    // definition of interval Poincare map
    IMap ivf("par:xi;var:x,y,z,w;fun:y,z,w,x*(1-x^2)-xi*z;");
    IOdeSolver isolver(ivf,8);
    ICoordinateSection isection(4,1);
    IPoincareMap ipm(isolver,isection);
    isolver.setAbsoluteTolerance(2e-7);
    isolver.setRelativeTolerance(2e-7);

    double L=0, xi=0;
    double delta = 1e-5;

    interval set0 = x + interval(-1,1);
    interval r;

    while(L<=end){

        double xi = L+0.5*delta;
        x = findOrbit_x(xi,x,pm); // we find an approximate zero of G

        if(xi >= 2.031){ // if we are close to bifurcation, we decrease the tolerance of ODE solver
            isolver.setAbsoluteTolerance(2e-8);
            isolver.setRelativeTolerance(2e-8);
        }

        double R = capd::min(L+delta, end.rightBound());
        interval XI(L,R);
        auto [proven_orbit,set] = proveOrbit_x(XI,x,ipm); // we use Interval Newton Method to determine the existence of the smooth curve within some box
        
        if(!(proven_orbit))  
            delta*=0.95; // if we failed, then we shrink \xi range
        else {
            if(!intersection(set0,set,r)) { // checking if curves of neighbouring intervales glue correctly into one curve
                cout << "x curve not glued correctly" << endl;
                cout << "set0=" << set0 << endl;
                cout << "set=" << set << endl;
                break;
            }
            L = L+delta;
            delta*=1.001; // if we succeded, then we slightly increase \xi range
            set0 = set;
        }
    }
    cout << "Whole curve computed" << endl;
    return {counter,set0};
}