#include "curve_x.h"

using namespace std;
using namespace capd;

/**
 * @brief Finds an approximate solution of G(\xi,x) = 0, where \xi is treated as the parameter
 * 
 * @param double parameter of vector field
 * @param double initial condition for ODE solver
 * @param DPoincareMap reference to Poincare map object
 * 
 * @returns approximate zero x
 */
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

/**
 * @brief proves the existence of a smooth curve in interval product xi \times X
 * 
 * @param interval xi range for parameter
 * @param double x middle point of interval X
 * @param IPoincareMap reference to interval Poincare map object
 * 
 * @returns tuple {res, X} - res is a bool value, which tells if proof was successful, X is a range of initial values, which encloses the curve
 */
tuple<bool,interval> proveOrbit_x(interval xi, double x, IPoincareMap& pm){
    static const interval sqrt2=interval(sqrt(2.0));
    static interval S = interval(-1,1)*8e-6; // interval [-eps,eps]
    static interval maxS = interval(-1,1)*5e-4;

    pm.getVectorField().setParameter(0,xi); 
    interval x0 = x;
    C0HORect2Set s({x0,0,(sqr(x0)-1)/sqrt2,0}); // initial vector as an input to the vector field
    interval X = x0 + S; // interval [x-eps,x+eps]
    interval y, N;
    IVector v1, v2;

    try{
        y = pm(s,2)[3];
        IMatrix D{{1,0,0,0},{0,0,0,0},{sqrt2*X,0,0,0},{0,0,0,0}};
        C1HORect2Set set(C1Rect2Set::C0BaseSet({X,0,(sqr(X)-1)/sqrt2,0}),C1Rect2Set::C1BaseSet{D});
        v1 = pm(set);
        v2 = pm(set,D);
        interval derivative = pm.computeDP(v2,D)[3][0];
        N = -y/derivative;
    }catch(exception e) {
        S *= 0.95; // we shrink [-eps,eps] interval if integration failed
        return {false,0};
    }
    bool geometric_ok = X<-1 and v1[0]>1 and v2[0]<1 and v2[0]>-1; // geometric conditions for periodic orbits
    // returns the result and X range needed to verify if the curve is glued correctly
    return {subset(N,S) and geometric_ok and intersection(maxS,N*interval(-1,1)*1.03,S),X};
}

/**
 * @brief determines the existence of the curve in range \xi \in [0,\xi_*]
 * 
 * @param double initial point for determining the curve, approximate zero of G
 * @param interval interval, which contains value \xi_*
 * 
 * @returns tuple {num_of_subintervals, X_*} - the amount of needed subintervals and enclosure of the point x(\xi_*)
 */
tuple<int,interval> proveCurve_x(double x, interval end) {  
    int counter_subintervals = 0;
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
            if(counter_subintervals++% 1000==0 || XI.contains(end)) 
                cout << "XI=" << XI << ", X=" << set << endl;

            L = L+delta;
            delta*=1.001; // if we succeded, then we slightly increase \xi range
            set0 = set;
        }
    }
    cout << "Whole curve computed" << endl;
    return {counter_subintervals,set0};
}