/*
 *  AtmScaling.h
 *  Created by Jeremy Riousset on 10/25/07.
 */

#ifndef ATMSCALING_H
#define ATMSCALING_H

#include <cmath>
#include <iostream>
#include <list>
using namespace std;
class Scaling
{
public:
    static double StdAtmosphere(double h);
    static double StdAtmosphere(double hh, list<double> _alt, list<double> _ng);
};

#endif
