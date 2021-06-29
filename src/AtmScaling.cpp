/*
 *  AtmScaling.cpp
 *  Created by Jeremy Riousset on 10/25/07.
 */

#include "AtmScaling.h"

double Scaling::StdAtmosphere(double _h)
{
    /**********************************************************************************/
    /*                            !!! Altitude in km !!!                                  */
    /**********************************************************************************/
    _h*=1e-3;

    double _dh,_result;
    int _UpperBound, _LowerBound, _Median;

    /**********************************************************************************/
    /* US Standard Atmosphere:                                                          */
    /*    We replaced original 2.5e19cm-3 at 0 km altitude by 2.688e19 cm-3 which is our*/
    /*    reference number density at temperature 273 K.                                  */
    /* This will be our altitude scaling. We prefer this to the exp(-z/8.4km) found in*/
    /* litertature, because the last one is temperature dependent.                      */
    /**********************************************************************************/

    double _Altitude[]={
        0.e+0,        5.e+0,        1.e+1,        1.5e+1,        2.e+1,
        2.5e+1,        3.e+1,        3.5e+1,        4.e+1,        4.5e+1,
        5e+1,        5.5e+1,        6e+1,        6.4e+1,        6.8e+1,
        7.2e+1,        7.6e+1,        8e+1,        8.4e+1,        8.8e+1,
        9.2e+1,        9.6e+1,        1e+2,        1.08e+2,    1.14e+2,
        1.2e+2,        1.26e+2,    1.32e+2,    1.4e+2,        1.5e+2};

    double _NeutralDensity[]={
        2.688e+19,    1.53e+19,    8.59e+18,    4.05e+18,    1.85e+18,
        8.33e+17,    3.83e+17,    1.76e+17,    8.31e+16,    4.088e+16,
        2.13e+16,    1.181e+16,    6.33e+15,    3.93e+15,    2.39e+15,
        1.39e+15,    7.72e+14,    4.03e+14,    1.99e+14,    9.48e+13,
        4.37e+13,    2.07e+13,    1.04e+13,    3.18e+12,    1.43e+12,
        6.61e+11,    3.4e+11,    1.91e+11,    9.7e+10,    4.92e+10};

    if (_h<0 || _h>150e3)
    {
        printf("h is outside of 0-150 km interval.\n");
        exit(2);
    }

    /**********************************************************************************/
    /* We have 30 samples partitioned in 30 intervals. First, we derive the upper/lower */
    /* bounds (_UpperBound/_LowerBound) of the interval                                  */
    /**********************************************************************************/

    _LowerBound=0;
    _UpperBound=29;
    while(_UpperBound!=_LowerBound+1)
    {
        _Median=(_UpperBound+_LowerBound)/2;
        _dh=_h-_Altitude[_Median];
        if(_dh<0)  _UpperBound=_Median;
        if(_dh>=0) _LowerBound=_Median;
    }

    /**********************************************************************************/
    /* Now, we do an interpolation to estimate the neutral density at any height      */
    /* between 0 and 150 km                                                              */
    /**********************************************************************************/

    _result = (_NeutralDensity[_LowerBound]/_NeutralDensity[0])*
        pow(_NeutralDensity[_UpperBound]/_NeutralDensity[_LowerBound],
            (_h-_Altitude[_LowerBound])/(_Altitude[_UpperBound]-_Altitude[_LowerBound]));
    //    double result=pow(10.,log10(_NeutralDensity[_LowerBound]/_NeutralDensity[0]) +
    //              (log10(_NeutralDensity[_UpperBound])-log10(_NeutralDensity[_LowerBound]))
    //              *(_h-_Altitude[_LowerBound])/(_Altitude[_UpperBound]-_Altitude[_LowerBound]));
    //    double c = 100*abs(_result-result)/_result;
    return _result;

    /**********************************************************************************/
    /* We note that result == _result analytically. However it does not seem to be the*/
    /* case here. The tests we tried with Matlab, show the same difference. We guess  */
    /* either an error in our reasoning or an approximation due to the use of log10 in*/
    /* the first case. We choose the first solution.                                  */
    /**********************************************************************************/
}

double Scaling::StdAtmosphere(double _h, list<double> _alt, list<double> _ng)
{
    /**********************************************************************************/
    /*                            !!! Altitude in km !!!                                  */
    /**********************************************************************************/
    _h*=1e-3;
    double _dh,_result;
    int _UpperBound, _LowerBound, _Median;
    int _N;
    if(_alt.size()!=_ng.size()) {
        printf("ee:: Dismension mismatch when reading atmospheric profiles.");
        exit(2);
    } else {
        _N = (int)_alt.size();
    }
    double _Altitude[_N];
    double _NeutralDensity[_N];
    copy( _alt.begin(), _alt.end(), _Altitude );
    copy( _ng.begin() , _ng.end(), _NeutralDensity );

    
    if (_h<_Altitude[0] || _h>_Altitude[_N-1])
    {
        printf("h=%6.2f km is outside of %6.2f-%6.2f km interval.\n",_h,_Altitude[0],_Altitude[_N-1] );
        exit(2);
    }

    /**********************************************************************************/
    /* We have _N samples partitioned in 30 intervals. First, we derive the upper/lower */
    /* bounds (_UpperBound/_LowerBound) of the interval                                  */
    /**********************************************************************************/

    _LowerBound=0;
    _UpperBound=_N-1;
    while(_UpperBound!=_LowerBound+1)
    {
        _Median=(_UpperBound+_LowerBound)/2;
        _dh=_h-_Altitude[_Median];
        if(_dh<0)  _UpperBound=_Median;
        if(_dh>=0) _LowerBound=_Median;
    }

    /**********************************************************************************/
    /* Now, we do an interpolation to estimate the neutral density at any height      */
    /**********************************************************************************/

    _result = (_NeutralDensity[_LowerBound]/_NeutralDensity[0])*
        pow(_NeutralDensity[_UpperBound]/_NeutralDensity[_LowerBound],
            (_h-_Altitude[_LowerBound])/(_Altitude[_UpperBound]-_Altitude[_LowerBound]));
    //    double result=pow(10.,log10(_NeutralDensity[_LowerBound]/_NeutralDensity[0]) +
    //              (log10(_NeutralDensity[_UpperBound])-log10(_NeutralDensity[_LowerBound]))
    //              *(_h-_Altitude[_LowerBound])/(_Altitude[_UpperBound]-_Altitude[_LowerBound]));
    //    double c = 100*abs(_result-result)/_result;
    return _result;

    /**********************************************************************************/
    /* We note that result == _result analytically. However it does not seem to be the*/
    /* case here. The tests we tried with Matlab, show the same difference. We guess  */
    /* either an error in our reasoning or an approximation due to the use of log10 in*/
    /* the first case. We choose the first solution.                                  */
    /**********************************************************************************/
}
