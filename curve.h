#include<iostream>
#include "capd/capdlib.h"

#ifndef _CURVE_H_
#define _CURVE_H_

long double find_orbit(long double, long double, capd::LDPoincareMap&, bool);

bool check_if_periodic(capd::Interval, capd::Interval, capd::IPoincareMap&, bool, bool);

void determine_curve_x(long double, long double);

void determine_curve_xi(long double, long double);

#endif