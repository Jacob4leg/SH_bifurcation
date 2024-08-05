#include<iostream>
#include "capd/capdlib.h"
#include<tuple>

#ifndef _CURVE_H_
#define _CURVE_H_

long double find_orbit(long double, long double, capd::LDPoincareMap&, bool);

std::tuple<bool,capd::Interval> check_if_periodic(capd::Interval, capd::Interval, capd::IPoincareMap&, bool);

void determine_curve_x(long double, long double);

void determine_curve_xi(long double, long double);

#endif