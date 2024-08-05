#include "curve.h"

using namespace std;
using namespace capd;

long double find_orbit(long double x, long double xi, LDPoincareMap &pm, bool mode /*1 for x, 0 for xi*/) {
    if(mode) {
        pm.getVectorField().setParameter("xi",xi);
    }
    
    int dim = pm.getVectorField().dimension();
    static long double sqrt2 = sqrt(2.);
    double eps = 1e-13;

    LDMatrix m1(dim, dim);
    LDMatrix m2(dim, dim);
    LDMatrix D(dim,dim);
    
    long double d;
    do {
        LDVector u0 = mode ? LDVector({x,0, (x*x - 1) / sqrt2, 0}) : LDVector({x,0, (x*x - 1) / sqrt2, 0, xi});
        LDVector u1 = pm(u0, m1);
        LDVector u2 = pm(u1, m2);
        
        D = pm.computeDP(u2,m2) * pm.computeDP(u1,m1);

        long double fx = u2[3];
        long double derivative = mode ? D[3][0] + D[3][2] * sqrt2 * x : D[3][4];

        d = fx / derivative;

        if(mode) { // we modify x
            x = x - d;
        }
        else{ // we modify xi
            xi = xi - d;
        }
    } while(abs(d) > eps);

    return mode ? x : xi;
}


// tuple<bool, Interval> interval_newton(Interval x, Interval xi, IPoincareMap &pm, bool mode, int iterations) {
//     int dim = pm.getVectorField().dimension();
//     static Interval two(2);
//     static Interval sqrt2 = sqrt(two);

//     IMatrix m(dim,dim);

//     IVector u0 = mode ? IVector({x,0, (power(x,2) - 1) / sqrt2, 0}) : IVector({x,0, (power(x,2) - 1) / sqrt2, 0, xi});


//     C1HORect2Set s(u0); 
//     IVector u2 = pm(s,m,2);
//     IMatrix D = pm.computeDP(u2,m);

//     Interval fx = u2[3];

//     Interval derivative = mode ? D[3][0] + D[3][2] * sqrt2 * x : D[3][4];

//     Interval N = - fx / derivative;

//     Interval X = mode ? x : xi;

//     Interval r;

//     if(intersection(N,X,r)) {
//         if(mode)
//             return interval_newton(X,xi,pm,mode,iterations - 1);
//         else
//             return interval_newton(x,X,pm,mode,iterations - 1);
//     }
// }


tuple<bool,Interval> check_if_periodic(Interval x, Interval xi, IPoincareMap &pm, bool mode /*1 for x, 0 for xi*/) {
    if(mode)
        pm.getVectorField().setParameter("xi",xi);
    
    static interval S = interval(-1,1)*1e-7;
    static interval maxS = interval(-1,1)*5e-4;


    // cout << "x=" << x << endl;
    // cout << "xi=" << xi << endl;
    int dim = pm.getVectorField().dimension();
    static Interval two(2);
    static Interval sqrt2 = sqrt(two);
    
    IMatrix m(dim,dim);


    IVector u0 = mode ? IVector({x,0, (power(x,2) - 1) / sqrt2, 0}) : IVector({x,0, (power(x,2) - 1) / sqrt2, 0, xi});

    C0HORect2Set s0(u0);
    
    IVector u1 = pm(s0);
    
    

    if(mode) {
        x = x + S;
    }
    else {
        xi = xi + S;
    }
    u0 = mode ? IVector({x,0, (power(x,2) - 1) / sqrt2, 0}) : IVector({x,0, (power(x,2) - 1) / sqrt2, 0, xi});

    C1HORect2Set s(u0); 
    IVector u2 = pm(s,m,2);
    IMatrix D = pm.computeDP(u2,m);

    bool geometry_ok = u0[0] < -1 && u0[2] > 0 &&
                        u1[0] > 1 && u1[2] < 0 &&
                        u2[0] > -1 && u2[0] < 1 && u2[2] > 0;
    
    if(!geometry_ok) {
        cout << "Wrong geometry for xi=" << xi << "\t x=" << x << endl;
        cout << "u0=" << u0[0] << "\t u0''=" << u0[2] << endl;
        cout << "u1=" << u1[0] << "\t u1''=" << u1[2] << endl;
        cout << "u2=" << u2[0] << "\t u2''=" << u2[2] << endl;
        
        return {false, 0};
    }

    Interval fx = u2[3];
    Interval derivative = mode ? D[3][0] + D[3][2] * sqrt2 * x : D[3][4];
    
    if(derivative.contains(0)) {
        cout << "derivative contains zero" << endl;
        return {false, 0};
    }


    Interval N = - fx / derivative;

    // Interval N = mode ? x.mid() - d : xi.mid() - d;

    // cout << "S=" << S << endl;
    // cout << "N=" << N << endl;
    // cout << subset(N,S) << endl;

    Interval res_set = mode ? x : xi;

    return {subset(N,S) && intersection(maxS,N*interval(-1,1)*1.03,S), res_set};
}

void determine_curve_x(long double x, long double xi) {
    int dim = 4;
    string vf_str = "par:xi;var:x,y,z,w;fun:y,z,w,x*(1-x^2)-xi*z;";

    LDMap vf(vf_str);
    LDOdeSolver solver(vf,20);
    LDCoordinateSection section(dim,1);
    LDPoincareMap pm(solver,section);

    IMap ivf(vf_str);
    IOdeSolver isolver(ivf,7);
    isolver.setAbsoluteTolerance(2e-7);
    isolver.setRelativeTolerance(2e-7);
    // isolver.setAbsoluteTolerance(1e-10);
    // isolver.setRelativeTolerance(1e-10);
    ICoordinateSection isection(dim,1);
    IPoincareMap ipm(isolver, isection);

    double L = xi;
    double step = 1e-7;
    
    int counter = 0;

    interval x_set = x;
    interval r;

    while(xi <= interval(204)/100) {

        xi =  L + 0.5 * step;

        x = find_orbit(x,xi,pm,1);

        bool is_periodic;
        auto period_checking_res = check_if_periodic(x,interval(L,L+step),ipm,1);

        is_periodic = get<0>(period_checking_res);
    
        

        // try {
        //     is_periodic = check_if_periodic(x,interval(L,L+step),ipm,1);
        // }
        // catch(poincare::PoincareException<C0HORect2Set>) {
        //     is_periodic = false;
        // }

        if(is_periodic) {
            if(!intersection(x_set,get<1>(period_checking_res),r)) {
                cout << "NOPE" << endl;
                return;
            }
            x_set = get<1>(period_checking_res);

            if(counter % 1000 == 0)
                cout << "xi=" << interval(L,L+step) << "\t x=" << x << endl;
            L += step;
            step *= 1.001;
            // if(xi >= 2.0315)
            //     step = 1e-7;
            counter++;
        }
        else{
            step *= 0.95;
            // if(xi >= 2.0315)
            //     step = 1e-7;
        }
    }
    cout << "xi=" << interval(L,L+step) << "\t x=" << x << endl;
    cout << counter << endl;
}

void determine_curve_xi(long double x, long double xi) {
    int dim = 5;
    string vf_str = "var:x,y,z,w,xi;fun:y,z,w,x*(1-x^2)-xi*z,0;";

    LDMap vf(vf_str);
    LDOdeSolver solver(vf,20);
    LDCoordinateSection section(dim,1);
    LDPoincareMap pm(solver,section);

    IMap ivf(vf_str);
    IOdeSolver isolver(ivf,7);
    isolver.setAbsoluteTolerance(2e-7);
    isolver.setRelativeTolerance(2e-7);
    // isolver.setAbsoluteTolerance(1e-10);
    // isolver.setRelativeTolerance(1e-10);
    ICoordinateSection isection(dim,1);
    IPoincareMap ipm(isolver, isection);

    double L = x;
    double step = 1e-7;
    

    int counter = 0;

    while(true) {
        // double xi_temp = xi;

        x = L + 0.5 * step;
        xi = find_orbit(x,xi,pm,0);
        bool is_periodic;
        auto period_checking_res = check_if_periodic(interval(L,L+step),xi,ipm,0);
        is_periodic = get<0>(period_checking_res);

        // bool is_max = xi < xi_temp;
        if(is_periodic) {
            if(counter % 1000 == 0)
                cout << "x=" << interval(L,L+step) << "\t xi=" << xi + interval(-1,1) * 1e-6 << endl;
            // if(is_max && xi_temp > 2.0316500)
            //     cout << "max=" << xi_temp + interval(-1,1) * 1e-6 << endl;
            L += step;
            step *= 1.001;
            // if(xi >= 2.0315)
            //     step = 1e-7;
            counter++;
        }
        else{
            step *= 0.95;
            // if(xi >= 2.0315)
            //     step = 1e-7;
        }
    }
}