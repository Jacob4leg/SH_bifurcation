#include<iostream>
#include "capd/capdlib.h"
#include<tuple>

#ifndef _NEWTON_OPERATORS_H_
#define _NEWTON_OPERATORS_H_

long double find_orbit(long double, long double, capd::LDPoincareMap&, bool);

std::tuple<bool, capd::Interval> interval_newton(capd::Interval, capd::Interval, capd::IPoincareMap&, bool);


std::tuple<bool,capd::Interval> check_if_periodic(capd::Interval, capd::Interval, capd::IPoincareMap&, bool);


#endif