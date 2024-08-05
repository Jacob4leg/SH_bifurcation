#include<iostream>
#include "capd/capdlib.h"
#include "curve.h"

using namespace std;
using namespace capd;
// double findOrbit(double xi, double x, DPoincareMap& pm){
//   pm.getVectorField().setParameter(0,xi);
//   DMatrix D{{1,0,0,0},{0,0,0,0},{sqrt(2.)*x,0,0,0},{0,0,0,0}}, T(4,4);
//   DVector u{x,0,(x*x-1)/sqrt(2.),0};
//   for(int i=0;i<2;++i){
//     u = pm(u,T);
//     D = T*D;
//   }
//   D = pm.computeDP(u,D);
//   return x - u[3]/D[3][0];
// }

// int counter = 0;
// bool proveOrbit(interval xi, double x, IPoincareMap& pm){
//   static const interval sqrt2=interval(sqrt(2.0));
//   static interval S = interval(-1,1)*8e-6;
//   static interval maxS = interval(-1,1)*5e-4;
//   pm.getVectorField().setParameter(0,xi);
//   interval x0 = x;
//   C0Rect2Set s({x0,0,(sqr(x0)-1)/sqrt2,0});
//   interval y = pm(s,2)[3];
//   interval X = x0 + S;
//   IMatrix D{{1,0,0,0},{0,0,0,0},{sqrt2*X,0,0,0},{0,0,0,0}};
//   C1Rect2Set set(C1Rect2Set::C0BaseSet({X,0,(sqr(X)-1)/sqrt2,0}),C1Rect2Set::C1BaseSet{D});
//   IVector v1 = pm(set);
//   IVector v = pm(set,D);
//   interval N = -y/pm.computeDP(v,D)[3][0];
//   if(++counter% 1024==0) cout << "xi=" << xi << ", N=" << N << ", S=" << S << endl;
//   return subset(N,S) and X<-1 and v1[0]>1 and v[0]<1 and v[0]>-1 and intersection(maxS,N*interval(-1,1)*1.03,S);
// }


int main() {
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
