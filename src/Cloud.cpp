/* Cloud.cpp */
#include "Cloud.h"

/**************************************************************************************/
/* Load tripole model                                                                 */
/**************************************************************************************/
void Cloud::LoadTripole(FILE * file, const double II1, const double II2, const double IIscreen)
{
    clock_t	_StartTime	= clock();
    clock_t	_EndTime	= clock();
    clock_t	_runTime	= clock();
    ListCharge::iterator  it;														// Table with all parameters of the charge configuation
    CriticalFields _OverShotEc((1+Var::ThresholdOvershoot)*Var::Ec.getParams()[0],
                               Var::Ec.getParams()[1],Var::Ec.getParams()[2],Var::Ec.getParams()[3],
                               Var::d,Var::N,1);
    int		_choice		= 2;	// Define choice of method to derive bd time
    double	TT			= 0;	// Bd Time
    double	rr			= 0;	// E/[(1+ThresholdOvershoot)*Einit]
    int		kk			= 0;	// Nb of iterations
    
    if(_choice == 0) // Bisection Method
    {
        double	TTl		= 0;				// _s
        double	TTr		= 180;				// _s
        double	TTav	= (TTl+TTr)/2;		// _s
        double	rr_tmp	= 0;				// return value of E/[(1+ThresholdOvershoot)*Einit]
        int		rr_l	= 0;				// =1 if rr_tmp>1 (ie E>Einit) ; =-1 else
        int		rr_r	= 0;				// =1 if rr_tmp>1 (ie E>Einit) ; =-1 else
        int		rr_av	= 0;				// =1 if rr_tmp>1 (ie E>Einit) ; =-1 else
        
        rr_tmp	= Cloud::EoverEk(TTl, _OverShotEc);
        rr_l	= (rr_tmp>=1)*1+(rr_tmp<1)*-1;
        rr_tmp	= Cloud::EoverEk(TTr, _OverShotEc);
        rr_r	= (rr_tmp>=1)*1+(rr_tmp<1)*-1;
        
        if(TTl > TTr) {Swap::DBL(TTl,TTr); Swap::INT(rr_l,rr_l);}
        while(fabs(2*(TTr-TTl)/(TTr+TTl))>Var::epsilon)
            //	while(rr_r>1e-5)
        {
            kk++;
            TTav	= (TTl+TTr)/2;
            rr_tmp	= Cloud::EoverEk(TTav, _OverShotEc);
            rr_av	= (rr_tmp>=1)*1+(rr_tmp<1)*-1;
            if(rr_l*rr_av > 0) { TTl = TTav; rr_l = rr_av; }
            else {TTr = TTav; rr_r = rr_av;}
        }
        TT = (TTl+TTr)/2;
    }
    else if(_choice == 1)	//Adaptive timestep
    {
        double	dt		= 60;
        double	rr_tmp	= 0;
        
        while(dt>Var::epsilon)
        {
            kk++;
            rr_tmp	= Cloud::EoverEk(TT+dt, _OverShotEc);
            rr		= (rr_tmp>=1)*1+(rr_tmp<1)*-1;
            while(rr>=0 && dt>Var::epsilon)
            {
                dt /=2;
                rr_tmp	= Cloud::EoverEk(TT+dt, _OverShotEc);
                rr		= (rr_tmp>=1)*1+(rr_tmp<1)*-1;
            }
            TT +=dt;
            //printf("[TT dt] = [%lf %f]\n",TT,dt);
        }
    }
    else if(_choice == 2)	// Linearity
    {
        double rr = 1;
        TT	= 1; // _s
        rr = Cloud::EoverEk(TT, _OverShotEc);
        TT *= 1/rr;
    }
    
    
    
    /**********************************************************************************/
    /* Estimate how much the initiation threshold is exceeded						  */
    /**********************************************************************************/
    rr = Cloud::EoverEk(TT, _OverShotEc);
    
    /**********************************************************************************/
    /* display results																  */
    /**********************************************************************************/
    IO::print(file,"T = " + to_string(TT) + " ; max(E/((1+alpha)*Ek)) = "  + to_string(rr) + "\n");
    IO::print(file,"   Nb of iterations    = " + to_string(kk+1) + "\n");
    IO::print(file,"ii:\t Initiation Threshold is exceeded by " + to_string((rr*(1+Var::ThresholdOvershoot)-1)*100) + "%\n");
    IO::print(file,"ii:\t Final values of the charge domains:\n");

    printf("ii:\t\t [Q1 Q2 Q3 Q4]       = [ ");
    for (it=Var::ChargeCfg.begin() ; it!=Var::ChargeCfg.end() ; it++)
        cout<<it->getParams()[0]<<" ";
    printf("] C\n");
    
    _EndTime = clock();
    _runTime = _EndTime - _StartTime;
    
    IO::print(file,"ii:\t Run time for layers load: " + to_string((double)_runTime/CLOCKS_PER_SEC) + "s\n");
}


/**************************************************************************************/
double Cloud::EoverEk(const double tt, CriticalFields _OverShotEc)
{
    double					_MaxRatio(0);
    double					_ratio(0);
    int						nn;														// Counter of iteration
    double					QQ, XXq,YYq,ZZq, RRq1,RRq2,RRq3;
    Charge					CC(Var::d,Var::N);										// Temporary charge structure for SOR derivation
    ListCharge::iterator	it;														// Table with all parameters of the charge configuation
    
    nn=0; QQ=0;
    for (it=Var::ChargeCfg.begin() ; it!=Var::ChargeCfg.end() ; it++)
    {
        /* retrieve tripole geometrical parameters */
        XXq = it->getParams()[1];
        YYq = it->getParams()[2];
        ZZq = it->getParams()[3];
        RRq1= it->getParams()[4];
        RRq2= it->getParams()[5];
        RRq3= it->getParams()[6];
        
        /* evaluate layers charge */
        if(nn==0) QQ = -Var::I2*tt;
        if(nn==1) QQ = (-Var::I1+Var::I2)*tt;
        if(nn==2) QQ = (Var::I1*tt);
        if(nn==3) QQ = -Var::Iscreen*tt - ( -Var::I2*tt + (-Var::I1+Var::I2)*tt + (Var::I1*tt));
        if(nn>=4) QQ = 0;
        
        /* update layers */
        it->disk(QQ, XXq,YYq,ZZq, RRq1,RRq3, Var::d,Var::N);
        CC+=*it;
        nn++;
    }
    
    /* solve fields */
    Var::SOR.init(Var::phi,Var::epsilon, Var::MaxStep, Var::d, Var::N, CC, Var::Un);
    BC::Apply(Var::BCtype,Var::phi,CC.rho,Var::d,Var::N);
    Var::SOR.Solve(Var::d,Var::N,Var::Un,Var::phi);
    
    /* Find ddeltaE */
    for(int kk=0; kk<Var::N.z; kk++) for(int jj=1; jj<Var::N.y-1; jj++) for(int ii=0; ii<Var::N.x; ii++)
    {
        _ratio = foo::Eijk(ii,jj,kk,Var::phi,Var::d,Var::N)[0]/_OverShotEc.initiation[kk];
        if(_ratio>=_MaxRatio)	_MaxRatio	= _ratio;
    };
    return _MaxRatio;
}
/**************************************************************************************/
