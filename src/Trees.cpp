/* Trees.cpp */
#include "Trees.h"

/**************************************************************************************/
/* Initiate Tree (Find initiation point)											  */
/**************************************************************************************/
bool Tree::Initiate(FILE * file, const int InitiationType, Point& InitiationPoint)
{
    bool	isInitiated = false;
    if(InitiationType == RANDOM)
    {
        list<Point>				_ListOfCandidates;
        list<Point> :: iterator it;
        Point					_Candidate;
        int						_CntPoints(0);
        int						_RandPoint;
        int						kk=0;
        
        for(int ii=1; ii<Var::N.x-1; ii++) for(int jj=1; jj<Var::N.y-1; jj++) for(int kk=1; kk<Var::N.z-1; kk++)
        {
            /* Choose the points which are likely to launch a flash	*/
            if(foo::Eijk(ii,jj,kk, Var::phi,Var::d,Var::N)[0] >= (1+0*Var::ThresholdOvershoot)*Var::Ec.initiation[kk])
            {
                _Candidate.i=ii;
                _Candidate.j=jj;
                _Candidate.k=kk;
                _ListOfCandidates.push_back(_Candidate);
                _CntPoints++;
            }
        };
        IO::print(file,"ii:\t Number of possible candidate endpoints: " + to_string((int)_CntPoints) + "\n");
        //IO::print(file,"ii:\t Maximum threshold overreach: " + to_string(maxPerExc) + "%\n");
        
        if(_CntPoints>0)
        {
            /* Pick a point randomly in this list */
            _RandPoint = rand()%_CntPoints+1;
            for(it = _ListOfCandidates.begin() ; it!=_ListOfCandidates.end() ; it++)
            {
                if(kk==_RandPoint)
                {
                    InitiationPoint = *it;
                    Var::Un[InitiationPoint.i][InitiationPoint.j][InitiationPoint.k] = 1;
                    Var::phi0 = Var::phi[InitiationPoint.i][InitiationPoint.j][InitiationPoint.k];
                    break;
                }
                kk++;
            }
            isInitiated = true;
        }
        else
        {
            isInitiated = false;
            cout<<"ii:\t Cannot initiate flash (insufficient field)\n";
        }
    }
    else if(InitiationType == AT_EMAX)
    {
        /******************/
        /* Init at Emax   */
        /******************/
        double Emax = 0;
        for(int ii=1; ii<Var::N.x-1; ii++) for(int jj=1; jj<Var::N.y-1; jj++) for(int kk=1; kk<Var::N.z-1; kk++)
            if(foo::Eijk(ii,jj,kk, Var::phi,Var::d,Var::N)[0] >= Var::Ec.initiation[kk])
                if(foo::Eijk(ii,jj,kk, Var::phi,Var::d,Var::N)[0] >= Emax)
                {
                    Emax=foo::Eijk(ii,jj,kk, Var::phi,Var::d,Var::N)[0];
                    InitiationPoint.init(ii,jj,kk);
                };
        Var::Un[InitiationPoint.i][InitiationPoint.j][InitiationPoint.k] = 1;
        Var::phi0 = Var::phi[InitiationPoint.i][InitiationPoint.j][InitiationPoint.k];
        if(foo::Eijk(InitiationPoint.i,InitiationPoint.j,InitiationPoint.k, Var::phi,Var::d,Var::N)[0] < Var::Ec.initiation[InitiationPoint.k])
        {
            isInitiated = false;
            IO::print(file,"ii:\t Cannot initiate flash (insufficient field)\n");
        }
        else {isInitiated = true;}
    }
    else if(InitiationType == AT_PREDEF_POS)
    {
        /******************/
        /* Fix Init Point */
        /******************/
        if(foo::Eijk(InitiationPoint.i, InitiationPoint.j, InitiationPoint.k, Var::phi,Var::d,Var::N)[0] >= Var::Ec.initiation[InitiationPoint.k])
        {
            Var::Un[InitiationPoint.i][InitiationPoint.j][InitiationPoint.k] = 1;
            Var::phi0 = Var::phi[InitiationPoint.i][InitiationPoint.j][InitiationPoint.k];
            isInitiated = true;
        }
        else
        {
            IO::print(file,"ii:\t Cannot initiate flash (insufficient field)\n");
            isInitiated = false;
        }
    }
    else if(InitiationType == AT_REL_EMAX)
    {
        /*************************/
        /* Init at Emax relative */
        /*************************/
        double Emax = 0;
        for(int ii=1; ii<Var::N.x-1; ii++) for(int jj=1; jj<Var::N.y-1; jj++) for(int kk=1; kk<Var::N.z-1; kk++)
            if(foo::Eijk(ii,jj,kk, Var::phi,Var::d,Var::N)[0] >= Var::Ec.initiation[kk])
                if(foo::Eijk(ii,jj,kk, Var::phi,Var::d,Var::N)[0]/Var::Ec.initiation[kk] >= Emax)
                {
                    Emax=foo::Eijk(ii,jj,kk, Var::phi,Var::d,Var::N)[0]/Var::Ec.initiation[kk];
                    InitiationPoint.init(ii,jj,kk);
                };
        Var::Un[InitiationPoint.i][InitiationPoint.j][InitiationPoint.k] = 1;
        Var::phi0 = Var::phi[InitiationPoint.i][InitiationPoint.j][InitiationPoint.k];
        if(foo::Eijk(InitiationPoint.i,InitiationPoint.j,InitiationPoint.k, Var::phi,Var::d,Var::N)[0] < Var::Ec.initiation[InitiationPoint.k])
        {
            isInitiated = false;
            IO::print(file,"ii:\t Cannot initiate flash (insufficient field)\n");
        }
        else {isInitiated = true;}
    }
    else
    {
        IO::print(file,"ee:\t Wrong Choice for initiation procedure\n");
        exit(3);
    }
    /* Display coordinates of initiation point */
    IO::print(file, "++:\t Discharge initiated at [" + to_string(InitiationPoint.i*Var::d.x/1e3) + " km, " + to_string(InitiationPoint.j*Var::d.y/1e3) + " km, " + to_string((InitiationPoint.k*Var::d.z+Var::z_gnd)/1e3) + " km]\n");
    //	SOR.Solve(d,N,Un,phi);
    
    /* Derive field and potential before discharge along the central vertical axis */
    for(int kk=0 ; kk<Var::N.z ; kk++)
    {
        Var::EzNum[kk]	= foo::Eijk((Var::N.x-1)/2,(Var::N.y-1)/2,kk,Var::phi,Var::d,Var::N)[3];
        Var::phiNum[kk]	= Var::phi((Var::N.x-1)/2,(Var::N.y-1)/2,kk);
    }
    
    /* Store Initiation/Ambient data */
    IO::write(Var::phiNum,	(char*)"phiNumBF.dat");
    IO::write(Var::EzNum,	(char*)"EnumBF.dat");
	IO::write(Var::phi_amb, (char*)"phiAmb.dat");
    foo::GlobalE(Var::phi_amb, Var::d, Var::N, -2); // Stores ambient electric field values
    
	CMatrix3D				rrho_amb(Var::N.x, Var::N.y, Var::N.z);
    CMatrix2D				rrho_amb_yz(Var::N.y, Var::N.z);
    CMatrix2D				rrho_amb_xz(Var::N.x, Var::N.z);
    
	rrho_amb    = Var::rho;
    for(int kk=0 ; kk<Var::N.z ; kk++)
    {
        for(int ii=0 ; ii<Var::N.x ; ii++)
            rrho_amb_xz[ii][kk]= foo::rhoijk(ii,(Var::N.y-1)/2,kk,Var::phi_amb,Var::d,Var::N)*1e9;
        
        for(int jj=0 ; jj<Var::N.y ; jj++)
        {
            rrho_amb_yz[jj][kk]= foo::rhoijk( (Var::N.x-1)/2,jj,kk, Var::phi_amb,	Var::d,Var::N)*1e9;
            for(int ii=0 ; ii<Var::N.x ; ii++)
                rrho_amb(ii,jj,kk)	*= 1e9;
        }
    }
	IO::write(rrho_amb,     (char*)"rhoAmb.dat");
    IO::write(rrho_amb_yz,	(char*)"rhoAmbYZ.dat");
    IO::write(rrho_amb_xz,	(char*)"rhoAmbXZ.dat");
    
    return isInitiated;
}
/**********************************************************************************/


/**********************************************************************************/
/* Grow the Tree																  */
/**********************************************************************************/
void Tree::Grow(FILE * file, bool AddNew)
{
    bool					_isGndConnected;									    // Check connection to the ground
    char					_strRho3D[50];											// File name for storage of charge density
    char					_strPhi3D[50];											// File name for storage of electrical potential
    char					_strUn3D[50];											// File name for storage of fixed potential cells
    int						_CntLinks(0);											// Current link iteration
    int						_endi;													// Tracks the most recently added link's x-position
    int						_endj;													// Tracks the most recently added link's y-position
    int						_endk;													// Tracks the most recently added link's z-position
	double					_BndErrorTmp = 0;										// Error at the boundaries at the current step
    double					_EsEnergyTmp = 0;										// Electrostatic energy at the current step
    double					_rhoDiffEnd = 0;										// Charge density difference before and after growth of current link at the end node
    double					_rhoDiffNeg = 0;										// Cumulative negative charge density difference before and after growth of current link
    double					_rhoDiffPos = 0;										// Cumulative positive charge density difference before and after growth of current link
    CMatrix1D				_TotalPotentialTmp(Var::N.z);							// Total potential on the central vertical axis at the current step
    CMatrix1D				_TotalEfieldTmp(Var::N.z);								// Total eField on the central vertical axis at the current step
    CMatrix2D				_phi2D_cha(Var::N.y,Var::N.z);							// Potential induced by the channel in the y-z plane
    CMatrix3D				_phi_cha_tmp;											// Potential due to the channel in the domain before update of the boundary conditions
    CMatrix3D				_rho_cha(Var::N.x,Var::N.y,Var::N.z);					// 3-D Matrix for channel induced charge density
    CMatrix3D				_rho(Var::N.x,Var::N.y,Var::N.z);					    // 3-D Matrix for channel induced charge density
    CMatrix3D				_rhoBefore(Var::N.x,Var::N.y,Var::N.z);					// 3-D Matrix for global charge density before growth of current link
    CMatrix3D				_rhoAfter(Var::N.x,Var::N.y,Var::N.z);					// 3-D Matrix for global charge density after growth of current link
    
	//CMatrix1D				rrho3D(Var::N.x*Var::N.y*Var::N.z);						// Total charge density everywhere at the current step
    int						nn = 0;													// Blind variable
    ListLink::iterator		it;													    // index on the list of established links
    Charge					CC;														// Charge transfer at the current step
    Potential				PP;														// Channel potential at the current step
    Vector					pp;														// x-, y-, z-components and norm of the dipole moment at the current step
    Vector					_BndError;												// Potential at the boundary before update, after update and relative difference between the two values
    
    nn = 0;
    
	for(int kk=0 ; kk<Var::N.z ; kk++)
    {
        _TotalEfieldTmp[kk]	= foo::Eijk((Var::N.x-1)/2,(Var::N.y-1)/2,kk,Var::phi,Var::d,Var::N)[3];
        _TotalPotentialTmp[kk]	= Var::phi((Var::N.x-1)/2,(Var::N.y-1)/2,kk);
        if(Var::isEsEnergyCalculated == true && Var::step3d == 0)
            for(int jj=0 ; jj<Var::N.y ; jj++)
                for(int ii=0 ; ii<Var::N.x ; ii++)
				{
                    _EsEnergyTmp += eps0*pow(foo::Eijk(ii,jj,kk,Var::phi,Var::d,Var::N)[0],2)/2*Var::d.x*Var::d.y*Var::d.z;
					_rhoBefore[ii][jj][kk] = foo::rhoijk(ii,jj,kk,Var::phi,Var::d,Var::N)*1e+9; //_nC
				}
        if(Var::isEsEnergyCalculated == false && Var::step3d != 0)
            for(int jj=0 ; jj<Var::N.y ; jj++) for(int ii=0 ; ii<Var::N.x ; ii++)
            {
                _rho[ii][jj][kk]		= foo::rhoijk(ii,jj,kk,Var::phi,Var::d,Var::N)*1e+9; //_nC
				_rhoBefore[ii][jj][kk]  = _rho[ii][jj][kk];
                nn++;
            };
        if(Var::isEsEnergyCalculated == true && Var::step3d != 0)
            for(int jj=0 ; jj<Var::N.y ; jj++) for(int ii=0 ; ii<Var::N.x ; ii++)
            {
                _EsEnergyTmp += eps0*pow(foo::Eijk(ii,jj,kk,Var::phi,Var::d,Var::N)[0],2)/2*Var::d.x*Var::d.y*Var::d.z;
                _rho[ii][jj][kk] = foo::rhoijk(ii,jj,kk,Var::phi,Var::d,Var::N)*1e+9; //_nC
                _rhoBefore[ii][jj][kk]  = _rho[ii][jj][kk];
				nn++;
            };
    }
    if(Var::step3d != 0)
    {
        IO::write(Var::step3d,(char*)"step3d.dat");
        //IO::write(rrho3D,					(char*)"rho3d0.dat");
    }
  
    pp = foo::DipoleMoment(Var::QchannelPlus,Var::phi_cha,Var::Un,Var::L,Var::d,Var::N);

    Var::ChannelPotential.push_back(Var::phi0);
    Var::DischargeDipoleMoment.push_back(pp);
    Var::CarriedCharge.push_back(Var::QchannelPlus);
    Var::EsEnergy.push_back(_EsEnergyTmp);
    Var::TotalEfield.push_back(_TotalEfieldTmp);
    Var::TotalPotential.push_back(_TotalPotentialTmp);
    
    
    /* Clear lists for multiple simulation runs */
    //Var::EstablishedLinks.clear();
    //Var::NumberOfCandidates.clear();
    //Var::MaximumCandidateOverreach.clear();
    
    
    _CntLinks = Var::NumLinks;
    while(AddNew==true && _CntLinks>=0 )//&& _CntLinks<55)
    {
        AddNew	= Tree::AddNewLink(file,Var::d,Var::N,Var::Un,Var::phi, Var::Ec,Var::Vd, Var::InitiationPoint,Var::EstablishedLinks,
                                   Var::isBndXingPossible,  Var::isRsDeveloped,
                                   Var::isLinkXingPossible, Var::isQMinimized);
        
        /* Check Connection to the ground */
        _isGndConnected = false;
        for(it = Var::EstablishedLinks.begin(); it != Var::EstablishedLinks.end() ; it++)
            if(it->end.k == 0)
            {
                _isGndConnected = true;
                break;
            }
        if(Var::isQMinimized==false)
        {
			
            BC::Update(Var::isFlashAccoutedInBC,Var::BCtype,Var::phi,Var::rhoAmbMin,Var::rhoAmbMax,Var::d,Var::N);
            Var::SOR.Solve(Var::d,Var::N,Var::Un,Var::phi);
			Var::phi_cha	= Var::phi-Var::phi_amb;
        };
        if(Var::isQMinimized==true)
        {
            if(_isGndConnected == false)
            {
                Var::phi0 = Tree::fMinSearch(file, Var::phi0, Var::QchannelPlus, Var::Vmin,Var::Vmax , Var::epsilon,Var::MaxStep, Var::phi_cha,Var::phi_amb ,Var::Un, Var::InitiationPoint,Var::EstablishedLinks, Var::d,Var::N);
            }
            if (_isGndConnected == true)
            {
                Var::phi0 = Tree::EqualizeAtGroundPotential(Var::epsilon,Var::MaxStep, Var::phi_cha,Var::phi_amb ,Var::Un, Var::InitiationPoint,Var::EstablishedLinks, Var::d,Var::N);
            }
            Var::ChannelPotential.push_back(Var::phi0);
            
            if(Var::isBCerrorCalculated == false)
            {
                BC::Update(Var::isFlashAccoutedInBC,Var::BCtype,Var::phi_cha,Var::rhoAmbMin,Var::rhoAmbMax,Var::d,Var::N);
                Var::phi	= Var::phi_amb+Var::phi_cha;
            }
            
            if(Var::isBCerrorCalculated ==true)
            {
                _phi_cha_tmp	= Var::phi_cha;
                BC::Update(Var::isFlashAccoutedInBC,Var::BCtype,Var::phi_cha,Var::rhoAmbMin,Var::rhoAmbMax,Var::d,Var::N);
                Var::phi	= Var::phi_amb+Var::phi_cha;
                
                _BndError.x = 0;
                _BndError.y = 0;
                _BndError.z = 0;
                for(int kk=0; kk<Var::N.z; kk++) for(int jj=0; jj<Var::N.y; jj++) for(int ii=0; ii<Var::N.x; ii++)
                    if( ii == 0 || ii == Var::N.x-1 ||jj == 0 || jj == Var::N.y-1 || kk == 0 || kk == Var::N.z-1)
                    {
                        _BndErrorTmp = fabs(Var::phi_cha(ii,jj,kk) - _phi_cha_tmp(ii,jj,kk));
                        if(_BndErrorTmp >= fabs(_BndError.x-_BndError.y))
                        {
                            _BndError.x = _phi_cha_tmp(ii,jj,kk);
                            _BndError.y = Var::phi_cha(ii,jj,kk);
                            _BndError.z = Var::phi_amb(ii,jj,kk);;
                        }
                    }
                Var::BndUpdateErrors.push_back(_BndError);
            }

            /**************************************************************************/
            /* The update of the BC introduce sharp gradient and artificial charge    */
            /* close to the boundary, which slightly violate the charge conservation. */
            /* Because potential of the channel at stage N is influenced by the BC due*/
            /* to the channel at step N-1. It is not possible to completely compensate*/
            /* this effect. It is not possible to derive the potential of the channel */
            /* and the potential at the boundary at the same time.					  */
            /* We can try an approximation, using the potential of the channel derived*/
            /* at step N, using BC at N-1. We then derive new BC, and re-run an SOR   */
            /* algorithm to minimize the error.										  */
            /**************************************************************************/
            // Solution //
            /*
             PP.init(Var::phi_cha,Var::Un);
             Var::SOR.init(Var::phi_cha, Var::epsilon,Var::MaxStep, Var::d, Var::N, PP, Var::Un);
             Var::SOR.Solve(Var::d,Var::N,Var::Un,Var::phi_cha);
             Var::phi = Var::phi_amb+Var::phi_cha;
             */
            
            _EsEnergyTmp = 0;
            if(Var::step3d != 0 && _CntLinks%Var::step3d==0)
            {
                sprintf(_strRho3D,"rho3d%d.dat", _CntLinks);
                sprintf(_strPhi3D,"phi3d%d.dat", _CntLinks);
                sprintf(_strUn3D ,"Un3d%d.dat" , _CntLinks);
                nn          = 0;
            }
			for(it=Var::EstablishedLinks.begin(); it!=Var::EstablishedLinks.end();it++)
			{   // Tracks the coordinates of the most recently added end-node
				_endi = it->end.i;
				_endj = it->end.j;
				_endk = it->end.k;
			}
            for(int kk=0 ; kk<Var::N.z ; kk++)
            {
                _TotalEfieldTmp[kk]	= foo::Eijk((Var::N.x-1)/2,(Var::N.y-1)/2,kk,Var::phi,Var::d,Var::N)[3];
                _TotalPotentialTmp[kk]	= Var::phi((Var::N.x-1)/2,(Var::N.y-1)/2,kk);
                if(Var::isEsEnergyCalculated == true && Var::step3d == 0)
                    for(int jj=0 ; jj<Var::N.y ; jj++) for(int ii=0 ; ii<Var::N.x ; ii++)
					{
                        _EsEnergyTmp += eps0*pow(foo::Eijk(ii,jj,kk,Var::phi,Var::d,Var::N)[0],2)/2*Var::d.x*Var::d.y*Var::d.z;
						_rhoAfter[ii][jj][kk] = foo::rhoijk(ii,jj,kk,Var::phi,Var::d,Var::N)*1e+9; //_nC
						if((_rhoAfter[ii][jj][kk] - _rhoBefore[ii][jj][kk])>1e-9)
							_rhoDiffPos+= _rhoAfter[ii][jj][kk] - _rhoBefore[ii][jj][kk];
						else if ((_rhoAfter[ii][jj][kk] - _rhoBefore[ii][jj][kk])<-1e-9)
							_rhoDiffNeg+= _rhoAfter[ii][jj][kk] - _rhoBefore[ii][jj][kk];	
						if((ii==_endi) && (jj==_endj) && (kk==_endk))
							_rhoDiffEnd = _rhoAfter[ii][jj][kk] - _rhoBefore[ii][jj][kk];	
						_rhoBefore[ii][jj][kk] = _rhoAfter[ii][jj][kk];
					}
                if(Var::isEsEnergyCalculated == false && Var::step3d != 0 && _CntLinks%Var::step3d==0)
                    for(int jj=0 ; jj<Var::N.y ; jj++) for(int ii=0 ; ii<Var::N.x ; ii++)
                    {
                        _rho[ii][jj][kk] = foo::rhoijk(ii,jj,kk,Var::phi,Var::d,Var::N)*1e+9; //_nC
                       	if((_rhoAfter[ii][jj][kk] - _rhoBefore[ii][jj][kk])>1e-9)
							_rhoDiffPos+= _rhoAfter[ii][jj][kk] - _rhoBefore[ii][jj][kk];
						else if ((_rhoAfter[ii][jj][kk] - _rhoBefore[ii][jj][kk])<-1e-9)
							_rhoDiffNeg+= _rhoAfter[ii][jj][kk] - _rhoBefore[ii][jj][kk];	
						if((ii==_endi) && (jj==_endj) && (kk==_endk))
							_rhoDiffEnd = _rhoAfter[ii][jj][kk] - _rhoBefore[ii][jj][kk];	
						_rhoBefore[ii][jj][kk] = _rho[ii][jj][kk];
						nn++;
						
                    };
                if(Var::isEsEnergyCalculated == true && Var::step3d != 0)
                    for(int jj=0 ; jj<Var::N.y ; jj++) for(int ii=0 ; ii<Var::N.x ; ii++)
                    {
                        _EsEnergyTmp += eps0*pow(foo::Eijk(ii,jj,kk,Var::phi,Var::d,Var::N)[0],2)/2*Var::d.x*Var::d.y*Var::d.z;
                        if(_CntLinks%Var::step3d==0)
                        {
                            _rho[ii][jj][kk] = foo::rhoijk(ii,jj,kk,Var::phi,Var::d,Var::N)*1e+9; //_nC
                            nn++;
                        }
						_rhoAfter[ii][jj][kk] = foo::rhoijk(ii,jj,kk,Var::phi,Var::d,Var::N)*1e+9; //_nC
						if((_rhoAfter[ii][jj][kk] - _rhoBefore[ii][jj][kk]) > 1e-9)
							_rhoDiffPos+= _rhoAfter[ii][jj][kk] - _rhoBefore[ii][jj][kk];
						if(((_rhoAfter[ii][jj][kk] - _rhoBefore[ii][jj][kk]) < -1e-9))
							_rhoDiffNeg+= _rhoAfter[ii][jj][kk] - _rhoBefore[ii][jj][kk];	
						if((ii==_endi) && (jj==_endj) && (kk==_endk))
							_rhoDiffEnd = _rhoAfter[ii][jj][kk] - _rhoBefore[ii][jj][kk];	
						_rhoBefore[ii][jj][kk] = _rhoAfter[ii][jj][kk];
                    };
            }
            if(Var::step3d != 0 && _CntLinks%Var::step3d==0) {
                Var::rho = foo::Globalrho(Var::phi,Var::d,Var::N);
                IO::write(_rho    , _strRho3D);
                //IO::write(Var::rho, _strRho3D);
                IO::write(Var::phi, _strPhi3D);
                IO::write(Var::Un,  _strUn3D );
                foo::GlobalE(Var::phi, Var::d, Var::N, _CntLinks); //NB: This also stores all components of Ex,Ey,Ez.
            }
            
        };
        
		/* Update and export results */
		pp = foo::DipoleMoment(Var::QchannelPlus,Var::phi_cha,Var::Un,Var::L,Var::d,Var::N);
		Var::DischargeDipoleMoment.push_back(pp);
		Var::CarriedCharge.push_back(Var::QchannelPlus);

		
		Var::TransportedRhoEnd.push_back(_rhoDiffEnd);
		Var::TransportedRhoNeg.push_back(_rhoDiffNeg);
		Var::TransportedRhoPos.push_back(_rhoDiffPos);
		_rhoDiffEnd = 0;
		_rhoDiffNeg = 0;
		_rhoDiffPos = 0;
		if(Var::isQMinimized==true){
			Var::EsEnergy.push_back(_EsEnergyTmp);
			Var::TotalEfield.push_back(_TotalEfieldTmp);
			Var::TotalPotential.push_back(_TotalPotentialTmp);
		}

        /* Write interim results after every tenth link is added. */
        if((_CntLinks % Var::step3d) == 0)
            Tree::StoreData(file);
		if((Var::curType == PROPAGATING) || (Var::curType != INTRA_CLOUD))
       	 	_CntLinks++;
        cout<<"ii:\t Nb of established links: "<<_CntLinks<<endl;
    }
    cout<<"Total Number of links: "<<_CntLinks<<endl;
    _rho_cha			= foo::Globalrho(Var::phi_cha,Var::d,Var::N);
    Var::QchannelPlus	= foo::ChannelChargePositive(_rho_cha,Var::Un,Var::d,Var::N);
    Var::QchannelMinus	= foo::ChannelChargeNegative(_rho_cha,Var::Un,Var::d,Var::N);
    Var::Qtot_af		= foo::TotalCharge(foo::Globalrho(Var::phi,Var::d,Var::N),Var::d,Var::N);
    
    bool FFlagInit=false;
    for(int ii=1; ii<Var::N.x-1; ii++) for(int jj=1; jj<Var::N.y-1; jj++) for(int kk=1; kk<Var::N.z-1; kk++)
        if(foo::Eijk(ii,jj,kk, Var::phi,Var::d,Var::N)[0] >= (1+Var::ThresholdOvershoot)*Var::Ec.initiation[kk])
        {
            //			cout<<"Initiation Threshold exceeded at: ["<<ii*d.x<<" "<<jj*d.y<<" "<<kk*d.z<<"]\n";
            FFlagInit=true;
            break;
        };
    if(FFlagInit==false) cout<<"*** No further initiation possible. ***\n";
    else{cout<<"*** Further initiation possible. ***\n";}
    if (Var::isBCerrorCalculated == true)
        cout<<"Error on the boundaries\n"<<setw(12)<<"phiCha(N)"<<setw(12)<<"phiCha(N+1)"<<setw(16)<<"phiAmb"<<endl;
};
/**************************************************************************************/


/**************************************************************************************/
/* Store main data																	  */
/**************************************************************************************/

/* SAM:	Remember, the IO class takes care of putting the results in the correct folder.
 *		Therefore, specifying file names such as "results/ChannelPotential.dat" is unnecessary.
 *
 * btw, SAM = "Slight Additions and Modifications"
 */
void Tree::StoreData(FILE * file)
{
    //clock_t		_StartTime = clock();
    //clock_t		_EndTime;
    //clock_t		_RunTime;
    
    CMatrix2D				EE2D(	Var::N.y, Var::N.z);
    CMatrix2D				EEx2D(	Var::N.y, Var::N.z);
    CMatrix2D				EEy2D(	Var::N.y, Var::N.z);
    CMatrix2D				EEz2D(	Var::N.y, Var::N.z);
    CMatrix2D				pphi2D(	Var::N.y, Var::N.z);
    ListCMatrix1D			CChargeLayers;
    ListCharge::iterator	it;
	double					tempEfield = 0;
    
    IO::print(file, "..: Writing results\n");
    
    // Derive required values before storage //
    Var::E		= foo::GlobalE(Var::phi,		Var::d,Var::N,-1);
    Var::rho	= foo::Globalrho(Var::phi,		Var::d,Var::N);
    
    for(int kk=0 ; kk<Var::N.z ; kk++)
    {
        Var::EzNum[kk]	= foo::Eijk((Var::N.x-1)/2,(Var::N.y-1)/2,kk,Var::phi,Var::d,Var::N)[3];
        Var::phiNum[kk]	= Var::phi((Var::N.x-1)/2,(Var::N.y-1)/2,kk);
        for(int jj=0 ; jj<Var::N.y ; jj++)
        {
            EE2D[jj][kk]	   =   foo::Eijk( (Var::N.x-1)/2,jj,kk, Var::phi,		Var::d,Var::N)[0];
            EEx2D[jj][kk]	   =   foo::Eijk( (Var::N.x-1)/2,jj,kk, Var::phi,		Var::d,Var::N)[1];
            EEy2D[jj][kk]	   =   foo::Eijk( (Var::N.x-1)/2,jj,kk, Var::phi,		Var::d,Var::N)[2];
            EEz2D[jj][kk]	   =   foo::Eijk( (Var::N.x-1)/2,jj,kk, Var::phi,		Var::d,Var::N)[3];
            pphi2D[jj][kk]	   =    Var::phi( (Var::N.x-1)/2,jj,kk);
            for(int ii=0 ; ii<Var::N.x ; ii++)
            {
                Var::rho(ii,jj,kk)	*= 1e9;
            	if(Var::E[ii][jj][kk]>tempEfield)
					tempEfield = Var::E[ii][jj][kk];
            }
        }
    }
    for (it=Var::ChargeCfg.begin() ; it!=Var::ChargeCfg.end() ; it++)
        CChargeLayers.push_back(it->getParams());
    
    // Store altitude of the ground plane //
    IO::write(Var::Ec.getParams()[3], (char*)"z_gnd.dat");
    
    // Store Values in main plane //
    IO::write(Var::rho,		(char*)"rho3d.dat");
    IO::write(EE2D,			(char*)"E2d.dat");
    IO::write(EEx2D,		(char*)"Ex2d.dat");
    IO::write(EEy2D,		(char*)"Ey2d.dat");
    IO::write(EEz2D,		(char*)"Ez2d.dat");
    IO::write(pphi2D,		(char*)"phi2d.dat");
    IO::write(Var::Un,      (char*)"Un3d.dat" );
    
    if(Var::curType != PROPAGATING){
		if(Var::curType == JET)
			IO::write((char*)"Jet", 	        (char*)"Type_Result.txt");
		else if(Var::curType == CLOUD_TO_GROUND)
			IO::write((char*)"Cloud-to-Ground", (char*)"Type_Result.txt");
		else if(Var::curType == HORIZONTAL)
			IO::write((char*)"Horizontal",	    (char*)"Type_Result.txt");
		else if(Var::curType == INTRA_CLOUD)
			IO::write((char*)"Intracloud",      (char*)"Type_Result.txt");
	}
    /**************************************************************************
     double	QchannelPlus;
     double	QchannelMinus;
     double	tmp;
     ListDouble LL1(IO::read("rho.dat"));
     ListDouble LL2(IO::read("rhoAmb.dat"));
     ListDouble::iterator itD1;
     ListDouble::iterator itD2;
     
     QchannelPlus  = 0;
     QchannelMinus = 0;
     tmp			  = 0;
     
     itD2= LL2.begin();
     for (itD1=LL1.begin() ; itD1!=LL1.end() ; itD1++)
     {
     tmp = (*itD1 - *itD2)*1e-9;
     if(tmp>=0)
     QchannelPlus  += tmp*(float)Var::d.x*(float)Var::d.y*(float)Var::d.z;
     if(tmp<=0)
     QchannelMinus += tmp*(float)Var::d.x*(float)Var::d.y*(float)Var::d.z;
     itD2++;
     }
     cout.precision(10);
     cout<<"Positive Charge in the channel (rounded)      : "<<QchannelPlus<<" C\n";
     cout<<"Negative Charge in the channel (rounded)      : "<<QchannelMinus<<" C\n"<<endl;
     
     QchannelPlus  = 0;
     QchannelMinus = 0;
     tmp			  = 0;
     for(int ii=1 ; ii<Var::N.x-1 ; ii++) for(int jj=1 ; jj<Var::N.y-1 ; jj++) for(int kk=1 ; kk<Var::N.z-1 ; kk++)
     {
     tmp = Var::rho(ii,jj,kk)-rrho_amb(ii,jj,kk);
     if (Var::Un(ii,jj,kk)!=0)
     {
     if(tmp>=0)
     QchannelPlus  += tmp*Var::d.x*Var::d.y*Var::d.z*1e-9;
     if(tmp<=0)
     QchannelMinus += tmp*Var::d.x*Var::d.y*Var::d.z*1e-9;
     }
     };
     cout<<"Positive Charge in the channel                : "<<QchannelPlus<<" C\n";
     cout<<"Negative Charge in the channel                : "<<QchannelMinus<<" C\n"<<endl;
     **************************************************************************/
    
    // Store channel potential, charge and dipole moment //
    IO::write(Var::ChannelPotential,		(char*)"ChannelPotentials.dat");
    IO::write(Var::CarriedCharge,			(char*)"CarriedCharge.dat");
	IO::write(Var::TransportedRhoEnd,		(char*)"TransportedRhoEnd.dat");
	IO::write(Var::TransportedRhoNeg,		(char*)"TransportedRhoNeg.dat");
	IO::write(Var::TransportedRhoPos,		(char*)"TransportedRhoPos.dat");
    IO::write(Var::DischargeDipoleMoment,	(char*)"DischargeDipoleMoment.dat");
    
    // Store values of the electrostatic energy in the domain (if calculated) //
    if(Var::isEsEnergyCalculated == true)
        IO::write(Var::EsEnergy,			(char*)"EsEnergy.dat");
    
    // Store values of the boundary errors in the domain (if calculated) //
    if(Var::isBCerrorCalculated == true)
        IO::write(Var::BndUpdateErrors,		(char*)"BndError.dat");
    
    // Store d, N and Initiation point //
    IO::write(Var::d.x, Var::d.y, Var::d.z,	(char*)"dxyz.dat");
    IO::write(Var::N.x, Var::N.y, Var::N.z,	(char*)"Nxyz.dat");
    IO::write(Var::InitiationPoint.i*Var::d.x, Var::InitiationPoint.j*Var::d.y, Var::InitiationPoint.k*Var::d.z, Var::InitR, (char*)"InitPoint.dat");
    
    // Store charge center parameters //
    IO::write(CChargeLayers,	(char*)"ChargeLayers.dat");
    
    // Store Values of the fields and potential on the principal axis //
    IO::write(Var::phiNum,			(char*)"phiNumAF.dat");
    IO::write(Var::EzNum,			(char*)"EnumAF.dat");
    IO::write(Var::Ec.initiation,	(char*)"Einitiation.dat");
    IO::write(Var::Ec.negative,		(char*)"EthNegative.dat");
    IO::write(Var::Ec.positive,		(char*)"EthPositive.dat");
    IO::write(Var::Vd.negative,		(char*)"VdNegative.dat");
    IO::write(Var::Vd.positive,		(char*)"VdPositive.dat");
    IO::write(Var::TotalPotential,	(char*)"TotalPotential.dat");
    IO::write(Var::TotalEfield,		(char*)"TotalEfield.dat");
    
    // Store Values of the field and potential everywhere //
    IO::write(Var::phi,				(char*)"phi.dat");
    //IO::write(Var::E,				(char*)"E.dat");
	if(tempEfield>Var::MaximumEfield){
		Var::MaximumEfield = tempEfield;
		IO::write(tempEfield,(char*)"MaximumEfield.dat");
	};
	// Store the list of Established Links //
    IO::write(Var::EstablishedLinks,(char*)"EstablishedLinks.dat");
}
/**************************************************************************************/


bool Tree::init(char * ffolder, char * ffile, ListLink& EEstablishedLinks)
{
    char fileName[200];
    
    char *  line        = NULL;
    size_t  line_length = 0;
    ssize_t nread;
    Link    LL;
    FILE *	file;
    IO::setPathName(ffolder);
    IO::getPathName(fileName);
    strcat(fileName, ffile);
    file = fopen(fileName, "r");
    
    if (file)
        while ((nread = getline(&line, &line_length, file)) != -1) {
            sscanf(line,"%5d %5d %5d %5d %5d %5d %lf %lf %lf %lf\n", &(LL.start.i), &(LL.start.j), &(LL.start.k), &(LL.end.i), &(LL.end.j), &(LL.end.k), &(LL.l), &(LL.efield), &(LL.deltaV), &(LL.proba) );
            EEstablishedLinks.push_back(LL);
        }
    
    free(line);
    fclose(file);
    return 0;
}

/**************************************************************************************/
bool Tree::AddNewLink(FILE * file, ResGrid _d, SizeGrid _N,
				CMatrix3D& UUn, CMatrix3D& _phi,
				CriticalFields& EEc, VoltageDrops& VVd,
				const Point& InitiationPoint, ListLink& EEstablishedLinks,
				bool iisBndXingPossible,  bool iisRsDeveloped,
				bool iisLinkXingPossible, bool iisQMinimized)
{
	clock_t		_StartTime = clock();
	clock_t		_EndTime;
	clock_t		_RunTime;

	Point				_start;
	Point				_end;
	Link				_Candidate;
	Link				_ChosenLink;
	ListLink			_ListOfCandidates;
	ListLink::iterator	it1;
	ListLink::iterator	it2;
	double				_CounterOfCandidates;
	double				_SumProba;
	double				rr;							// random number in [0,1]

	double				linkAlt;					// Altitude (in meters) of endpoint of current link.
	bool				flag = false;				// is given the value "true" if the
													// candidate is crossing an establ-
													// ished link. Assume no crossing
													// link at start.
	bool				isFlash3D = true;			// equal true for  a 3-D flash
													// equal false for a 1-D flash used for testing purposes
	int                 cpt=0;

	/* SAM variables. */
	double overreach;								// The amound by which a candidate exceeds the electric field
													//  propagation threshold.
	double max_overreach(0);						// max_overreach = max(overreach), over all candidates.
	ListLink::iterator it=EEstablishedLinks.begin();

    IO::print(file,"..: Attempting to add new link\n");

	/* Seed the random number generator. */
	srand((int)(time(NULL) * 101));

	_start = InitiationPoint;
    
	if(isFlash3D == true) do{
	/**********************************************************************************/
	/* At this point we are sure the starting point is available for linking		  */
	/**********************************************************************************/
	for(_end.i = _start.i-1 ; _end.i <= _start.i+1 ; _end.i++)
	{
		for(_end.j = _start.j-1 ; _end.j <= _start.j+1 ; _end.j++)
		{
			for(_end.k = _start.k-1 ; _end.k <= _start.k+1 ; _end.k++)
			{
				if( _end.i >= 0 && _end.i <= _N.x-1 &&
					_end.j >= 0 && _end.j <= _N.y-1 &&
					_end.k >= 0 && _end.k <= _N.z-1 )
					if(UUn[_end.i][_end.j][_end.k] == 0)
					{
					/**********************************************************************************/
					/* At this point we are sure the ending point is available for linking			  */
					/**********************************************************************************/
						_Candidate.init(_start,_end, 0, _d,_phi[_start.i][_start.j][_start.k], _phi[_end.i][_end.j][_end.k]);
					/**********************************************************************************/
					/* Now we establish a candidate link and check if crossing any established link	  */
					/**********************************************************************************/
						if(_Candidate.efield >= EEc.positive[_Candidate.end.k] ||
						   _Candidate.efield <= EEc.negative[_Candidate.end.k])
						{
							if(iisLinkXingPossible == false)
							{
								flag = false;
								if(_Candidate.type == x || _Candidate.type == y || _Candidate.type == z)
									flag = false;

								if(_Candidate.type == xy)
									for(it1 = EEstablishedLinks.begin() ; it1 != EEstablishedLinks.end() ; it1++)
										if(_Candidate.start.k == (*it1).start.k && _Candidate.start.k == (*it1).end.k &&
										   (_Candidate.start.i + _Candidate.end.i) == ((*it1).start.i + (*it1).end.i) &&
										   (_Candidate.start.j + _Candidate.end.j) == ((*it1).start.j + (*it1).end.j) )
											flag = true;

								if(_Candidate.type == yz)
									for(it1 = EEstablishedLinks.begin() ; it1 != EEstablishedLinks.end() ; it1++)
										if(_Candidate.start.i == (*it1).start.i && _Candidate.start.i == (*it1).end.i &&
										   (_Candidate.start.j + _Candidate.end.j) == ((*it1).start.j + (*it1).end.j) &&
										   (_Candidate.start.k + _Candidate.end.k) == ((*it1).start.k + (*it1).end.k) )
											flag = true;

								if(_Candidate.type == xz)
									for(it1 = EEstablishedLinks.begin() ; it1 != EEstablishedLinks.end() ; it1++)
										if(_Candidate.start.j == (*it1).start.j && _Candidate.start.j == (*it1).end.j &&
										   (_Candidate.start.i + _Candidate.end.i) == ((*it1).start.i + (*it1).end.i) &&
										   (_Candidate.start.k +_Candidate. end.k) == ((*it1).start.k + (*it1).end.k) )
											flag = true;

								if(_Candidate.type == xyz)
									for(it1 = EEstablishedLinks.begin() ; it1 != EEstablishedLinks.end() ; it1++)
										if( (_Candidate.start.i + _Candidate.end.i) == ((*it1).start.i + (*it1).end.i) &&
											(_Candidate.start.j + _Candidate.end.j) == ((*it1).start.j + (*it1).end.j) &&
											(_Candidate.start.k + _Candidate.end.k) == ((*it1).start.k + (*it1).end.k) )
											flag = true;

								if(flag == false)
								{
								/**********************************************************************************/
								/* At this stage _Candidate contains a viable candidate for further propagation	  */
								/**********************************************************************************/

								/* SAM procedure.  Used to determine the greatest amount by which a
								 * candidate exceeds the propagation threshold.
								 */
									if(_Candidate.efield >= EEc.positive[_Candidate.end.k])
										overreach = 100 * (_Candidate.efield - EEc.positive[_Candidate.end.k])/EEc.positive[_Candidate.end.k];
									else
										overreach = 100 * (_Candidate.efield - EEc.negative[_Candidate.end.k])/EEc.negative[_Candidate.end.k];

									if(overreach > max_overreach)
										max_overreach = overreach;


									// This sets the voltage drop wrt. the initiation point //
									if(_Candidate.start == InitiationPoint)
										_Candidate.deltaV =
											(_Candidate.efield>=0)*VVd.positive[_Candidate.end.k]*_Candidate.l+
											(_Candidate.efield<=0)*VVd.negative[_Candidate.end.k]*_Candidate.l;
                                    else {
										_Candidate.deltaV = it->deltaV +
											(_Candidate.efield>=0)*VVd.positive[_Candidate.end.k]*_Candidate.l+
											(_Candidate.efield<=0)*VVd.negative[_Candidate.end.k]*_Candidate.l;
                                    }
							/*		// This section has only been designed for testing purposes //
									if(_Candidate.start == InitiationPoint)
										_Candidate.deltaV =
										(_Candidate.efield>=0)*1+
										(_Candidate.efield<=0)*-1;
									else
										_Candidate.deltaV = it->deltaV +
										(_Candidate.efield>=0)*1+
										(_Candidate.efield<=0)*-1;
							*/		// Add the Link to the list of candidates //
									_ListOfCandidates.push_back(_Candidate);
									// Increment counter of candidates //
									cpt++;
								}
							}
							else
							{
								_ListOfCandidates.push_back(_Candidate);
								cpt++;
							}
						}
					}
				}
			}
		}

		if(_start==InitiationPoint)
		{
			it = EEstablishedLinks.begin();
			_start = it->end;
		}
		else
		{
			it++;
			_start = it->end;
		}
	}while(it!=EEstablishedLinks.end());

	if(isFlash3D == false) do{
	/**********************************************************************************/
	/* At this point we are sure the starting point is available for linking		  */
	/**********************************************************************************/
		_end.i = _start.i;
		_end.j = _start.j;
		for(_end.k = _start.k-1 ; _end.k <= _start.k+1 ; _end.k++)
			if(	_end.k >= 0 && _end.k <= _N.z-1 )
				if(UUn[_end.i][_end.j][_end.k] == 0)
				{
					/**********************************************************************************/
					/* At this point we are sure the ending point is available for linking			  */
					/**********************************************************************************/
					_Candidate.init(_start,_end, 0, _d,_phi[_start.i][_start.j][_start.k], _phi[_end.i][_end.j][_end.k]);
					/**********************************************************************************/
					/* Now we establish a candidate link and check if crossing any established link	  */
					/**********************************************************************************/
						if(_Candidate.efield >= EEc.positive[_Candidate.end.k] ||
						   _Candidate.efield <= EEc.negative[_Candidate.end.k])
						{
							if(iisLinkXingPossible == false)
							{
								flag = false;
								if(flag == false)
								{
								/**********************************************************************************/
								/* At this stage _Candidate contains a viable candidate for further propagation	  */
								/**********************************************************************************/
									// This sets the voltage drop wrt. the initiation point //
									if(_Candidate.start == InitiationPoint)
										_Candidate.deltaV =
											(_Candidate.efield>=0)*VVd.positive[_Candidate.end.k]*_Candidate.l+
											(_Candidate.efield<=0)*VVd.negative[_Candidate.end.k]*_Candidate.l;
									else
										_Candidate.deltaV = it->deltaV +
											(_Candidate.efield>=0)*VVd.positive[_Candidate.end.k]*_Candidate.l+
											(_Candidate.efield<=0)*VVd.negative[_Candidate.end.k]*_Candidate.l;
									// Add the Link to the list of candidates //
									_ListOfCandidates.push_back(_Candidate);
									// Increment counter of candidates //
									cpt++;
								}
							}
							else
							{
								_ListOfCandidates.push_back(_Candidate);
								cpt++;
							}
						}
					};

		if(_start==InitiationPoint)
		{
			it = EEstablishedLinks.begin();
			_start = it->end;
		}
		else
		{
			it++;
			_start = it->end;
		}
	}while(it!=EEstablishedLinks.end());

	if(max_overreach > ANOMALOUS_OVERREACH)
	    IO::print(file, "AA: Anomalous maximum candidate overreach of " + to_string(max_overreach) + "% encountered\n");

	/**********************************************************************************/
	/* Derive probability of propagation for each candidate.						  */
	/**********************************************************************************/
	_SumProba = 0;
	_CounterOfCandidates = 0;
	for (it1 = _ListOfCandidates.begin() ; it1 != _ListOfCandidates.end() ; it1++)
	{
		if( (*it1).efield >= EEc.positive[(*it1).end.k] )
			(*it1).proba = pow(fabs((*it1).efield -EEc.positive[(*it1).end.k]),eta);
		else if ( (*it1).efield <= EEc.negative[(*it1).end.k] )
			(*it1).proba = pow(fabs((*it1).efield -EEc.negative[(*it1).end.k]),eta);
		else{
			IO::print(file, "ee:\t This link should not exist. There is an error in the code!!!\n");
            IO::print(file, "ee:\t Bad link encountered\n");
            IO::print(file, "ii:\t\t Number of links added so far: " + to_string((int)Var::NumLinks) + "\n");
            IO::print(file, "xx:\t Program exiting after fatal error\n");
			exit(1);
		}
		_SumProba	+= (*it1).proba;
		_CounterOfCandidates++;
	}

	if(_CounterOfCandidates == 0)
	{
		_EndTime = clock();
		_RunTime = _EndTime-_StartTime;
        IO::print(file, "ii:\t No more candidates\n");
        IO::print(file, "ii:\t Run time for Link addition: " + to_string((double)_RunTime/CLOCKS_PER_SEC) + " s\n");

        IO::print(file, "ii:\t Simulation has run out of candidates\n");
        IO::print(file, "ii:\t\t Number of links added so far: " + to_string((int)Var::NumLinks) + "\n");

		if(Var::maxAlt > ANOMALOUS_HEIGHT)
		{
			Var::curType = JET;
            IO::print(file, "ii:\t Discharge exceeded anomalous height of " + to_string(ANOMALOUS_HEIGHT/1e3) + " km\n");
            IO::print(file, "..:\t\t Classifying discharge as a jet\n");
		}
		else
		{
			Var::curType = INTRA_CLOUD;
            IO::print(file, "ii:\t Classifying discharge as intracloud\n");
		}
		return false;
    } else {
        Var::NumLinks++;
    }

	/* Keep track of the number of candidates and maximum candidate overreach.
	 * Keeping track of the number of candidates was the motivation to define the
	 * 'ListInt' data type.
	 */
    IO::print(file, "ii:\t Number of candidates: " + to_string((int)_CounterOfCandidates) + "\n");
    Var::NumberOfCandidates.push_back(_CounterOfCandidates);
    
    IO::print(file, "ii:\t Maximum candidate overreach: " + to_string(max_overreach) + "%\n");
	Var::MaximumCandidateOverreach.push_back(max_overreach);

	/**********************************************************************************/
	/* Convert probability to a	segment of length between 0 and 1, using the proba as */
	/* defined till now would increase the difficulty to choose a weightered by the   */
	/* probas. Hence we convert the proba to ease the choice of a link. Link i is	  */
	/* given the proba sum(p(j), j=1..i), consequently each proba is still between 0  */
	/* and 1.																		  */
	/* We randomly choose a point in this interval, we choose the link which proba is */
	/* the smallest greater than the random value previously defined				  */
	/**********************************************************************************/
	it2 = _ListOfCandidates.begin();
	(*it2).proba /= _SumProba;
	for(it1=_ListOfCandidates.begin() ; it2 != _ListOfCandidates.end() ; it1++)
	{
		it2 = it1;
		it2++;
		(*it2).proba /= _SumProba;
		(*it2).proba += (*it1).proba;
	}

	rr = random()/(double)RAND_MAX;
	it1 = _ListOfCandidates.begin();
	while(it1 != _ListOfCandidates.end() && (*it1).proba <= rr) {it1++;};
	_ChosenLink = *it1;
	_ListOfCandidates.clear();

	/**********************************************************************************/
	/* The chosen link is always the one propagating in the highest electric field,   */
	/* (i.e., the one with the highest propability to propagate). This way, the		  */
	/* stochasticity is eleminated.													  */
	/**********************************************************************************/
/*
		rr = 0;
		it1 = _ListOfCandidates.begin();
		for(it1 = _ListOfCandidates.begin() ; it1 != _ListOfCandidates.end() ; it1++)
		{
			if( it1->proba >= rr)
			{
				rr			= it1->proba;
				_ChosenLink = *it1;
			};
		}
		_ListOfCandidates.clear();
*/
	/**********************************************************************************/
	/* We finally fix the link in Un and add it to the list of established links and  */
	/* update phi.																	  */
	/**********************************************************************************/
	UUn[_ChosenLink.end.i][_ChosenLink.end.j][_ChosenLink.end.k] =
		UUn[_ChosenLink.start.i][_ChosenLink.start.j][_ChosenLink.start.k];
	if (iisQMinimized == true)
	{
		// The two following lines are not important since potential is rederived
		// afterwards to account for potential drop and overall neutrality in command
		// pphi0 = fMinSearch(pphi0, QQchannelPlus, VVmin,VVmax , eepsilon,MMaxStep, pphi_cha,pphi_amb ,UUn, InitiationPoint,EEstablishedLinks, _d,_N);
		// in main.cpp
//		_phi[_ChosenLink.end.i][_ChosenLink.end.j][_ChosenLink.end.k] =
//		_phi[_ChosenLink.start.i][_ChosenLink.start.j][_ChosenLink.start.k];
	}
	else if (iisQMinimized == false)
	{
		_phi[_ChosenLink.end.i][_ChosenLink.end.j][_ChosenLink.end.k] =
		_phi[_ChosenLink.start.i][_ChosenLink.start.j][_ChosenLink.start.k] -
		( (_ChosenLink.efield<0) * (EEc.negative[_ChosenLink.end.k]) +
		  (_ChosenLink.efield>0) * (EEc.positive[_ChosenLink.end.k]) ) * _ChosenLink.l;
	}
	EEstablishedLinks.push_back(_ChosenLink);

	/**********************************************************************************/
	/* Send authorization of propagation if a link has been established and that it   */
	/* does not reach the upper/lower boundary.										  */
	/* If we have _N.z points, the matrix is indexed from 0 to _N.z-1, as we here     */
	/* consider that boundaries cannot be linked, the accessible range is 1.._N.z-2   */
	/**********************************************************************************/
	_EndTime = clock();
	_RunTime = _EndTime - _StartTime;

	/**********************************************************************************/
	/* Check boundary crossing														  */
	/**********************************************************************************/
	if(_N.IsOnBoundary(_ChosenLink.end))
		if(iisBndXingPossible == false)
		{
            IO::print(file,"ii:\t Channel reached a boundary at [" + to_string(_ChosenLink.end.i) + " " + to_string(_ChosenLink.end.j) + " " + to_string(_ChosenLink.end.k) + "]\n");
            IO::print(file,"ii:\t Run time for Link addition: " + to_string((double)_RunTime/CLOCKS_PER_SEC) + " s\n");
            IO::print(file, "ii:\t\t Number of links added so far: " + to_string((int)Var::NumLinks) + "\n");

			if(_ChosenLink.end.k == (Var::N.z - 1))
			{
				Var::curType = JET;
				IO::print(file, "ii:\t Discharge reached top of domain\n");
                IO::print(file, "..:\t\t Classifying discharge as a jet\n");
			}
			else if(_ChosenLink.end.k == 0)
			{
				Var::curType = CLOUD_TO_GROUND;
                IO::print(file, "ii:\t Discharge reached bottom of domain\n");
                IO::print(file, "..:\t\t Classifying discharge as a cloud-to-ground strike\n");
			}
			else
			{
				Var::curType = HORIZONTAL;
                IO::print(file, "ii:\t Discharge reached a side of the domain\n");
                IO::print(file, "..:\t\t Classifying discharge as \"horizontal\"\n");
			}

			return false;
		};

	/**********************************************************************************/
	/* Allow/Prevent Return Stroke DeVeloPmenT										  */
	/**********************************************************************************/
	linkAlt = _ChosenLink.end.GetZ();
    IO::print(file, "ii:\t Link termination position: [" + to_string((_ChosenLink.end.i*Var::d.x)/1e3) + " " + to_string((_ChosenLink.end.j*Var::d.y)/1e3) + " " + to_string((_ChosenLink.end.k*Var::d.z+Var::z_gnd)/1e3)+ "] km\n");
	if(linkAlt > Var::maxAlt)
		Var::maxAlt = linkAlt;

	if(linkAlt > ANOMALOUS_HEIGHT)
	    IO::print(file, "AA:\t\t Altitude of link indicates possibility of jet development\n");
	
	if(_ChosenLink.end.k == 0 || _ChosenLink.end.k == _N.z-1)
		if(iisRsDeveloped == true)
		{
			for(int ii=0 ; ii<_N.x ; ii++) for(int jj=0 ; jj<_N.y ; jj++)
			{
				// Turn channel potential to 0 //
				for(int kk=0 ; kk<_N.z ; kk++)
					if(UUn[ii][jj][kk] == UUn[_ChosenLink.end.i][_ChosenLink.end.j][_ChosenLink.end.k])
						_phi[ii][jj][kk] = 0;
				// Prevent changes in ground/ionospheric potential //
				UUn[ii][jj][_ChosenLink.end.k] = UUn[_ChosenLink.end.i][_ChosenLink.end.j][_ChosenLink.end.k];
				_phi[ii][jj][_ChosenLink.end.k] = 0;
			};
            IO::print(file, "ii:\t Return stroke developed!\n");
			IO::print(file, "ii:\t\t Number of links added so far: " + to_string((int)Var::NumLinks) + "\n");
		};
	
    IO::print(file, "ii:\t Run time for Link addition: " + to_string((double)_RunTime/CLOCKS_PER_SEC) + " s\n");
	return true;
}

/**************************************************************************************/
/* Adjust potential in the channel													  */
/**************************************************************************************/
double Tree::Qchannel(const double VV,
					  const double eepsilon, const int MMaxStep,
					  CMatrix3D& pphi_cha, CMatrix3D& pphi_amb, CMatrix3D& UUn,
					  const Point& InitiationPoint, ListLink& EEstablishedLinks,
					  ResGrid _d, const SizeGrid& _N)
{
	ListLink::iterator it;

	Potential	P1;//(pphi_cha,UUn);
	SorSolution	SSOR;//(pphi_cha, eepsilon,MMaxStep, _d, _N, P1, UUn);

	pphi_cha[InitiationPoint.i][InitiationPoint.j][InitiationPoint.k]	= VV - pphi_amb[InitiationPoint.i][InitiationPoint.j][InitiationPoint.k];
	for(it=EEstablishedLinks.begin() ; it!= EEstablishedLinks.end() ; it++)
		pphi_cha[it->end.i][it->end.j][it->end.k]	= VV - it->deltaV - pphi_amb[it->end.i][it->end.j][it->end.k];

	P1.init(pphi_cha,UUn);
	SSOR.init(pphi_cha, eepsilon,MMaxStep, _d, _N, P1, UUn);
	SSOR.Solve(_d,_N,UUn,pphi_cha);

	double QQ = foo::ChannelCharge(foo::Globalrho(pphi_cha,_d,_N),UUn,_d,_N);
	return QQ;
}
/**************************************************************************************/

/**************************************************************************************/
/* Derive Channel Potential to minimize total Charge - Dichotomy vs. Nelder-Mead	  */
/**************************************************************************************/
double Tree::fMinSearch(FILE * file, const double VV, const double QQchannelPlus,
						const double VVmin, const double VVmax,
						const double eepsilon, const int MMaxStep,
						CMatrix3D& pphi_cha, CMatrix3D& pphi_amb, CMatrix3D& UUn,
						const Point& InitiationPoint, ListLink& EEstablishedLinks,
						ResGrid _d, const SizeGrid& _N)
{
	double choice = 0;
    if (choice == 0) // Bisection Method //
	{
		clock_t _StartTime = clock();
		double	CC;
		double	CCl	= VVmin;
		double	CCr = VVmax;
		double  CCav((CCr+CCl)/2);
		double	QQr = Qchannel( CCr , eepsilon,MMaxStep, pphi_cha,pphi_amb,UUn, InitiationPoint,EEstablishedLinks, _d,_N);
		double	QQl = Qchannel( CCl , eepsilon,MMaxStep, pphi_cha,pphi_amb,UUn, InitiationPoint,EEstablishedLinks, _d,_N);
		double	QQav= Qchannel( CCav, eepsilon,MMaxStep, pphi_cha,pphi_amb,UUn, InitiationPoint,EEstablishedLinks, _d,_N);
		int		kk	= 0;

		if(CCl > CCr)
		{
			Swap::DBL(CCl,CCr);
			Swap::DBL(QQl,QQr);
		}

//		while(fabs(2*(CCr-CCl)/(CCr+CCl))>eepsilon)
		while(fabs(QQav/QQchannelPlus)>eepsilon && fabs(2*(CCr-CCl)/(CCr+CCl))>eepsilon)
		{
			kk++;
			CCav = (CCr+CCl)/2;
			QQav = Qchannel( CCav, eepsilon,MMaxStep, pphi_cha,pphi_amb,UUn, InitiationPoint,EEstablishedLinks, _d,_N);
			if( QQl*QQav > 0 )
			{
				CCl = CCav;
				QQl = QQav;
			}
			else
			{
				CCr = CCav;
				QQr = QQav;
			}
//			cout<<"[Vl Vr] = ["<<setw(12)<<CCl<<" "<<setw(12)<<CCr<<"] ; [Ql Qr] = ["<<setw(12)<<QQl<<setw(13)<<QQr<<"] ; Nb of iterations = "<<kk<<endl;
//			cout<<"|QQav/QQchannelPlus| = "<<setw(12)<<fabs(QQav/QQchannelPlus)<<endl;
		}
//		if(fabs(QQav/QQchannelPlus)<=eepsilon)		cout<<"==============================>> break on CHARGE   "<<endl;
//		if(fabs(2*(CCr-CCl)/(CCr+CCl))<=eepsilon)	cout<<"==============================>> break on POTENTIAL"<<endl;

		CC		= (fabs(QQr)>fabs(QQl))*CCl+(fabs(QQl)>fabs(QQr))*CCr; //choose value which minimize QQ
//		QQav	= (fabs(QQr)>fabs(QQl))*QQl+(fabs(QQl)>fabs(QQr))*QQr;
		QQav = Qchannel( CCav, eepsilon,MMaxStep, pphi_cha,pphi_amb,UUn, InitiationPoint,EEstablishedLinks, _d,_N);
		clock_t _EndTime = clock();
		clock_t _RunTime = _EndTime - _StartTime;

//		double QQcha_cha = ChannelCharge(Globalrho(pphi_cha,_d,_N),UUn,_d,_N);
//		double QQcha_tot = TotalCharge(Globalrho(pphi_cha,_d,_N),_d,_N);
//		cout<<"QQcha_cha = "<<QQcha_cha<<endl;
//		cout<<"QQcha_tot = "<<QQcha_tot<<endl;
//		cout<<"Qdiff     = " << QQcha_tot-QQcha_cha<<endl;

        IO::print(file,"ii:\t Run time for Qminimization: " + to_string((double)_RunTime/CLOCKS_PER_SEC) + "s\n");
//		cout<<"Vc = "<<setw(12)<<CC<<" ; Qc/Qc+ = "<<setw(12)<<QQav<<" ; Nb of iterations = "<<kk<<endl;
		return CC;
	}
	else // Nelder-Mead //
	{
		clock_t _StartTime = clock();
		double	CC;
		double	x1(VVmin)	, f1(0);
		double	x2(VVmax)	, f2(0);
//		double	x1((1-1e-5)*VV)	, f1;
//		double	x2((1+1e-5)*VV)	, f2;
		double	xr(0)		, fr(0);
		double	xe(0)		, fe(0);
		double	xc(0)		, fc(0);
		double	xcc(0)		, fcc(0);
		int		kk = 0;


		while(fabs(2*(x1-x2)/(x1+x2))>eepsilon /*&& kk<25*/)
//		while(fabs(Qchannel((x1+x2)/2, eepsilon,MMaxStep, pphi_cha,pphi_amb,UUn, InitiationPoint,EEstablishedLinks, _d,_N))>=1e-6)
		{
			kk++;
			if		(kk != 1 && x1 == xr )	{f1 = fr ;}
			else if	(kk != 1 && x1 == xe )	{f1 = fe ;}
			else if	(kk != 1 && x1 == xc )	{f1 = fc ;}
			else if	(kk != 1 && x1 == xcc)	{f1 = fcc;}
			else							{f1 = fabs(Qchannel(x1, eepsilon,MMaxStep, pphi_cha,pphi_amb,UUn, InitiationPoint,EEstablishedLinks, _d,_N));}

			if		(kk != 1 && x2 == xr )	{f2 = fr ;}
			else if	(kk != 1 && x2 == xe )	{f2 = fe ;}
			else if	(kk != 1 && x2 == xc )	{f2 = fc ;}
			else if	(kk != 1 && x2 == xcc)	{f2 = fcc;}
			else							{f2 = fabs(Qchannel(x2, eepsilon,MMaxStep, pphi_cha,pphi_amb,UUn, InitiationPoint,EEstablishedLinks, _d,_N));}
			// Step 1: Order //
			if(f1 > f2)
			{
				Swap::DBL(x1,x2);
				Swap::DBL(f1,f2);
			}
			// Step 2: Reflect //
			xr = 2.*x1-x2;
			if		(kk != 1 && xr == x1 )	{fr = f1 ;}
			else if	(kk != 1 && xr == x2 )	{fr = f2 ;}
			else if	(kk != 1 && xr == xe )	{fr = fe ;}
			else if	(kk != 1 && xr == xc )	{fr = fc ;}
			else if	(kk != 1 && xr == xcc)	{fr = fcc;}
			else							{fr = fabs(Qchannel(xr, eepsilon,MMaxStep, pphi_cha,pphi_amb,UUn, InitiationPoint,EEstablishedLinks, _d,_N));}
			// Step 3: Expand //
			if(fr < f1)
			{
				xe = 3.*x1-2.*x2;
				if		(kk != 1 && xe == x1 )	{fe = f1 ;}
				else if	(kk != 1 && xe == x2 )	{fe = f2 ;}
				else if	(kk != 1 && xe == xr )	{fe = fr ;}
				else if	(kk != 1 && xe == xc )	{fe = fc ;}
				else if	(kk != 1 && xe == xcc)	{fe = fcc;}
				else							{fe = fabs(Qchannel(xe, eepsilon,MMaxStep, pphi_cha,pphi_amb,UUn, InitiationPoint,EEstablishedLinks, _d,_N));}

				if(fe<fr) {x2 = xe;}
				else {x2 = xr;} // (fe>=fr)
			}
			// Step 4: Contract //
			else
			{
				// Step 4a) Outside //
				if(f1<=fr && fr<f2)
				{
					xc = 3/2.*x1-1/2.*x2;
					if		(kk != 1 && xc == x1 )	{fc = f1 ;}
					else if	(kk != 1 && xc == x2 )	{fc = f2 ;}
					else if	(kk != 1 && xc == xr )	{fc = fr ;}
					else if	(kk != 1 && xc == xe )	{fc = fe ;}
					else if	(kk != 1 && xc == xcc)	{fc = fcc;}
					else							{fc = fabs(Qchannel(xc, eepsilon,MMaxStep, pphi_cha,pphi_amb,UUn, InitiationPoint,EEstablishedLinks, _d,_N));}

					if(fc<=fr) {x2=xc;}
					// Step 5: Shrink Step //
					else {x2 = x1 + 1/2.*(x2-x1);}
				}
				// Step 4b) Inside //
				if(fr>=f2)
				{
					xcc = 1/2.*x1+1/2.*x2;
					if		(kk != 1 && xcc == x1)	{fcc = f1 ;}
					else if	(kk != 1 && xcc == x2)	{fcc = f2 ;}
					else if	(kk != 1 && xcc == xr)	{fcc = fr ;}
					else if	(kk != 1 && xcc == xe)	{fcc = fe ;}
					else if	(kk != 1 && xcc == xc)	{fcc = fc;}
					else							{fcc = fabs(Qchannel(xcc, eepsilon,MMaxStep, pphi_cha,pphi_amb,UUn, InitiationPoint,EEstablishedLinks, _d,_N));}

					if(fcc<f2) {x2=xcc;}
					// Step 5: Shrink Step //
					else {x2 = x1 + 1/2*(x2-x1);}
				}
			}
//			cout<<"[x1 x2] = ["<<setw(12)<<x1<<" "<<setw(12)<<x2<<"] ; Q = "<<setw(12)<<fabs(Qchannel((x1+x2)/2, eepsilon,MMaxStep, pphi_cha,pphi_amb,UUn, InitiationPoint,EEstablishedLinks, _d,_N))<<" ; Nb of iterations = "<<kk<<endl;
		}
//		cout<<"*** FinaL ***"<<endl<<"[x1 x2] = ["<<setw(12)<<x1<<" "<<setw(12)<<x2<<"] ;\nQ = "<<setw(12)<<ChannelCharge(Globalrho(pphi_cha,_d,_N), UUn, _d,_N)<<endl;
//		cout<<"Qt = "<<setw(12)<<TotalCharge(Globalrho(pphi_cha+pphi_amb,_d,_N), _d,_N)<<endl;
		clock_t _EndTime = clock();
		clock_t _RunTime = _EndTime - _StartTime;

        IO::print(file,"ii:\t Run time for Qminimization: " + to_string((double)_RunTime/CLOCKS_PER_SEC) + "s\n");
        IO::print(file,"Vc = " + to_string(x1) + " ; Qc = " + to_string(f1) + " ; Nb of iterations = " + to_string(kk) + "\n");
        return CC=x1;
	}
}
/**************************************************************************************/

/**************************************************************************************/
double Tree::EqualizeAtGroundPotential(const double eepsilon, const int MMaxStep,
									   CMatrix3D& pphi_cha, CMatrix3D& pphi_amb, CMatrix3D& UUn,
									   const Point& InitiationPoint, ListLink& EEstablishedLinks,
									   ResGrid _d, const SizeGrid& _N)
{
	ListLink::iterator	it;
	Potential			P1;		//(pphi_cha,UUn);
	SorSolution			SSOR;	//(pphi_cha, eepsilon,MMaxStep, _d, _N, P1, UUn);

	pphi_cha[InitiationPoint.i][InitiationPoint.j][InitiationPoint.k]	= -pphi_amb[InitiationPoint.i][InitiationPoint.j][InitiationPoint.k];
	for(it=EEstablishedLinks.begin() ; it!= EEstablishedLinks.end() ; it++)
		pphi_cha[it->end.i][it->end.j][it->end.k]	= -pphi_amb[it->end.i][it->end.j][it->end.k];
	P1.init(pphi_cha,UUn);
	SSOR.init(pphi_cha, eepsilon,MMaxStep, _d, _N, P1, UUn);
	SSOR.Solve(_d,_N,UUn,pphi_cha);
	return 0; // Total potential of the channel phi = phi_amb + phi_cha = phi_amb + (-phi_amb) = 0;
}
/**************************************************************************************/
