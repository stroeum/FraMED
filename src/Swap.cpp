/*
 *  Swap.cpp
 *  Created by Jeremy Riousset on 10/25/07.
 */

#include "Swap.h"

void Swap::DBL(double& x1, double& x2)
{
	double tmp(x1);
	x1 = x2;
	x2 = tmp;
}

void Swap::INT(int& n1, int& n2)
{
	int tmp(n1);
	n1 = n2;
	n2 = tmp;
}
