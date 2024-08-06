#include<iostream>
#include "capd/capdlib.h"
#include "curve.h"

using namespace std;
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
    // IMatrix D{{1,0,0,0},{0,0,0,0},{sqrt2*X,0,0,0},{0,0,0,0}};
    IMatrix m(4,4),D(4,4);
    C1HORect2Set set({X,0,(sqr(X)-1)/sqrt2,0});
    // C1Rect2Set set(C1Rect2Set::C0BaseSet({X,0,(sqr(X)-1)/sqrt2,0}),C1Rect2Set::C1BaseSet{D});
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
    if(counter++% 100==0 || xi.contains(2031635./1000000)) cout << "xi=" << xi << ", X=" << X << ", N=" << N << ", S=" << S << endl;

    bool is_subset = subset(N,S);
    // bool geometric_ok = X<-1 and v1[0]>1 and v2[0]<1 and v2[0]>-1;
    bool geometric_ok = true;


    return {is_subset && geometric_ok && intersection(maxS,N*interval(-1,1)*1.03,S),X};
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

bool proveOrbit_xi(double xi, interval x, IC2PoincareMap &pm) {
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
        return false;
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
        return false;
    }
    pm.computeDP(v2,DPhi,HPhi,D,H);

    Interval fxi = v2[3];
    Interval G_xi = D[3][4];

    Interval N = fxi / G_xi;

    Interval G_x = D[3][0] + D[3][2] * sqrt2 * x;

    // cout << "G(x,xi)=" << fxi << endl;
    // cout << "dG/dxi=" << G_xi << endl;
    // cout << "dG/dx=" << G_x << endl;

    Interval xi_x = - G_x / G_xi;
    // cout << "dxi/dx=" << xi_x << endl;
    // cout << "dxi/dx=" << xi_x << endl;
    Interval G_x_x = H(3,0,0) + sqrt2 * (H(3,0,2) + D[3][2]);
    Interval G_xi_xi = H(3,4,4);
    Interval G_x_xi = H(3,0,4);

    Interval xi_x_x = -(G_x_x + (2 * G_x_xi + G_xi_xi * xi_x) * xi_x) / G_xi;
    // cout << "d^2xi/dx^2=" << xi_x_x << endl;
    return subset(N,S) and (xi_x_x < 0) and intersection(maxS,N*interval(-1,1)*1.03,S);
}


void fun(double x) {
    cout.precision(14);

    DMap vf("par:xi;var:x,y,z,w;fun:y,z,w,x*(1-x^2)-xi*z;");
    DOdeSolver solver(vf,20);
    DCoordinateSection section(4,1);
    DPoincareMap pm(solver,section);

    IMap ivf("par:xi;var:x,y,z,w;fun:y,z,w,x*(1-x^2)-xi*z;");
    IOdeSolver isolver(ivf,7);
    ICoordinateSection isection(4,1);
    IPoincareMap ipm(isolver,isection);
    isolver.setAbsoluteTolerance(2e-7);
    isolver.setRelativeTolerance(2e-7);

    long double L=0, xi=0;
    long double delta = 1e-5;
    // double x = -1.0845890871772264;

    interval set0 = x + interval(-1,1);
    interval r;
    interval end = interval(2031635)/1000000;

    while(L<=end){

        double xi = L+0.5*delta;
        x = findOrbit_x(xi,x,pm);

        if(L + delta > end) {
            auto res_orbit_check = proveOrbit_x(interval(L,end.rightBound()),x,ipm);
            bool proven_orbit = get<0>(res_orbit_check);
            interval set = get<1>(res_orbit_check);

            if(!intersection(set0,set,r)) {
                cout << "x curve not glued correctly" << endl;
                break;
            }

            if(!proven_orbit) {
                cout << "NOT COOL" << endl;
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

int main() {
    // xi=[2.0316354390905, 2.0316354456206], N=[-2.2880109633725e-06, 2.2849418823175e-06], S=[-2.34688121145e-06, 2.34688121145e-06], step=[6.5300120866141e-09, 6.5300120866141e-09]

    // fun(-1.0849406631521703);

    fun(-1.0845890871772264);
    return 0;
    // good end_xi = 2.0316303819895247, x=-1.5825456493039376

    cout.precision(17);
    double x1 = -1.0849406631521703;
    double x2 = -1.0845890871772264;
    double xi = 0.0;
    // cout << "first_curve" << endl;
    // determine_curve_x(x1,xi);
    // xi=[2.0316266806418004, 2.0316266806917382]      x=-1.5825492625615398

    // double x0 = -1.5825492625615398;
    // double xi0 = interval(2.0316266806418004, 2.0316266806917382).mid().leftBound();
    // determine_curve_xi(x0,xi0);

    // probably max max_xi=[2.0316506135713217, 2.031652613571322]

    cout << "second_curve" << endl;
    determine_curve_x(x2,xi);
}
