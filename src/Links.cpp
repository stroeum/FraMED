/* Links.cpp */

#include "Links.h"

/**************************************************************************************/
/* Point																			  */
/**************************************************************************************/
bool Point::operator==(const Point& Pt) const
{
	if (i == Pt.i && j == Pt.j && k == Pt.k)
		return true;
	else
		return false;
}  // Operator ==

bool Point::operator!=(const Point& Pt) const
{
	if (i != Pt.i || j != Pt.j || k != Pt.k)
		return true;
	else
		return false;
}  // Operator !=

Point& Point::operator=(const Point& Pt)
{
	i = Pt.i;
	j = Pt.j;
	k = Pt.k;
	return *this;
} // Operator =
/**************************************************************************************/

/**************************************************************************************/
/* Link																				  */
/**************************************************************************************/
Link::Link(const Link& LL)
{Link::init(LL);}

Link::Link(Point sstart, Point eend, double ddeltaV, StepsSizes dd, double pphiStart, double pphiEnd)
{Link::init(sstart,eend, ddeltaV, dd, pphiStart,pphiEnd);}

Link::Link(int sstarti, int sstartj, int sstartk, int eendi, int eendj, int eendk, double ddeltaV, StepsSizes dd, double pphiStart, double pphiEnd)
{Link::init(sstarti,sstartj,sstartk, eendi,eendj,eendk, ddeltaV, dd, pphiStart, pphiEnd);}

Link::Link(Point sstart, Point eend, double ddeltaV, StepsSizes dd, CMatrix3D pphi)
{Link::init(sstart, eend, ddeltaV, dd, pphi);}

Link::Link(int sstarti, int sstartj, int sstartk, int eendi, int eendj, int eendk, double ddeltaV, StepsSizes dd, CMatrix3D pphi)
{Link::init(sstarti,sstartj,sstartk, eendi,eendj,eendk, ddeltaV, dd, pphi);}

void Link::init(const Link& LL)
{
	start	= LL.start;
	end		= LL.end;
	type	= LL.type;
	l		= LL.l;
	deltaV	= LL.deltaV;
	efield	= LL.efield;
	proba   = LL.proba;
}

void Link::init(Point sstart, Point eend, double ddeltaV, StepsSizes dd, double pphiStart, double pphiEnd)
{
	start	= sstart;
	end		= eend;
	
	if (start == end)
		cout<<"Link can be only be established if starting and ending points are differents. ;)\n";
	
	//	l = sqrt(pow( (end.i-start.i)*dd.x,2) + pow((end.j-start.j)*dd.y,2) + pow( (end.k-start.k)*dd.z,2) );
	
	if(start.k == end.k)
	{
		if		(start.j == end.j)	{type = x;		l = dd.x;  }				// Link along x
		else if (start.i == end.i)	{type = y;		l = dd.y;  }				// Link along y
		else						{type = xy;		l = dd.xy; }				// Link in plane xy
	}
	else if (start.j == end.j)
	{
		if		(start.i == end.i)	{type = z;		l = dd.z;  }
		else if (start.k == end.k)	{type = x;		l = dd.x;  }
		else						{type = xz;		l = dd.xz; }
	}
	else if (start.i == end.i)
	{
		if		(start.j == end.j)	{type = z;		l = dd.z;  }
		else if (start.k == end.k)	{type = y;		l = dd.y;  }
		else						{type = yz;		l = dd.yz; }
	}
	else							{type = xyz;	l = dd.xyz;}
	
	deltaV = ddeltaV;
	efield = -(pphiEnd-pphiStart)/l;
} // Endof init

void Link::init(int sstarti, int sstartj, int sstartk, int eendi, int eendj, int eendk, double ddeltaV, StepsSizes dd, double pphiStart, double pphiEnd)
{
	start.i = sstarti;
	start.j = sstartj;
	start.k = sstartk;
	
	end.i = eendi;
	end.j = eendj;
	end.k = eendk;
	
	Link::init(start,end, ddeltaV, dd, pphiStart, pphiEnd);
} // Endof init

void Link::init(Point sstart, Point eend, double ddeltaV, StepsSizes dd, CMatrix3D pphi)
{
	double pphiStart = pphi[start.i][start.j][start.k];
	double pphiEnd	 = pphi[end.i][end.j][end.k];

	Link::init(start,end, ddeltaV, dd, pphiStart, pphiEnd);
} // Endof init

void Link::init(int sstarti, int sstartj, int sstartk, int eendi, int eendj, int eendk, double ddeltaV, StepsSizes dd, CMatrix3D pphi)
{	start.i = sstarti;
	start.j = sstartj;
	start.k = sstartk;
	
	end.i = eendi;
	end.j = eendj;
	end.k = eendk;
	
	Link::init(start,end, ddeltaV, dd, pphi);
} // Endof init


void Link::fixLink(CMatrix3D& UUn)
{UUn[end.i][end.j][end.k] = UUn[start.i][start.j][start.k];}

bool Link::operator==(const Link& LLink) const
{
	if (start == LLink.start && end == LLink.end)
		return true;
	else if (start == LLink.end && end == LLink.start)
		return true;
	else
		return false;
} // Operator ==

bool Link::operator!=(const Link& LLink) const
{
	if (start == LLink.start && end == LLink.end)
		return false;
	else if (start == LLink.end && end == LLink.start)
		return false;
	else
		return true;
} // Operator !=

Link& Link::operator=(const Link& LL)
{
	start	= LL.start;
	end		= LL.end;
	l		= LL.l;
	efield	= LL.efield;
	type	= LL.type;
	deltaV	= LL.deltaV;
	proba	= LL.proba;
	return *this;
} // operator=

bool Link::isCrossing(const Link& LL) const
{
	/**********************************************************************************/
	/* Note: The two next conditions are an introduction compared to V. Pasko's code. */
	/* However, they are pointless since link between 2 points already linked is      */
	/* forbidden.																	  */
	/**********************************************************************************/
	
	if (start == LL.start && end == LL.end)
		return true;								// twice the same link
	if (start == LL.end && end == LL.start)
		return true;								// twice the same link
		
	/**********************************************************************************/
	/* Note: With the correction above, when a 1-D, 2-D, or 3-D link is tested with   */
	/* itself the value return is 1. In Pasko et al. 2000, this was only true 2-D,    */ 
	/* 3-D links																	  */
	/**********************************************************************************/
	
	if(type == x || type == y || type == z) 
	   return false;								// impossible to cross a 1-D link
	
	if(type == xy && LL.type == xy && start.k == LL.start.k)
	{
		if( (start.i + end.i) == (LL.start.i + LL.end.i) &&
			(start.j + end.j) == (LL.start.j + LL.end.j))
			return true;
		else
			return false;
	}												// 2 links in same xy plane

	if(type == xz && LL.type == xz && start.j == LL.start.j)
	{
		if( (start.i + end.i) == (LL.start.i + LL.end.i) &&
			(start.k + end.k) == (LL.start.k + LL.end.k))
			return true;
		else
			return false;
	}
	
	if(type == yz && LL.type == yz && start.i == LL.start.i)
	{
		if( (start.j + end.j) == (LL.start.j + LL.end.j) &&
			(start.k + end.k) == (LL.start.k + LL.end.k))
			return true;
		else
			return false;
	}
	
	if(type == xyz && LL.type == xyz)
	{
		if( (start.i + end.i) == (LL.start.i + LL.end.i) &&
			(start.j + end.j) == (LL.start.j + LL.end.j) &&
			(start.k + end.k) == (LL.start.k + LL.end.k) )
			return true;							// return true if the 2 links are identical
		else
			return false;
	}
	return false;
}

bool Link::isCrossing(ListLink& LL) const
{
	ListLink::iterator it;

	if(type == x || type == y || type == z)
	{return false;}
	
	if(type == xy)
	{
		for(it = LL.begin() ; it != LL.end() ; it++)
			if(start.k == (*it).start.k && start.k == (*it).end.k)
				if( (start.i + end.i) == ((*it).start.i + (*it).end.i) &&
					(start.j + end.j) == ((*it).start.j + (*it).end.j))
					return true;
		return false;
	}												// 2 links in same xy plane
		
	if(type == yz)	
	{
		for(it = LL.begin() ; it != LL.end() ; it++)
			if(start.i == (*it).start.i && start.i == (*it).end.i)
				if( (start.j + end.j) == ((*it).start.j + (*it).end.j) &&
					(start.k + end.k) == ((*it).start.k + (*it).end.k))
					return true;				
		return false;
	}
	
	if(type == xz)
	{
		for(it = LL.begin() ; it != LL.end() ; it++)
			if(start.j == (*it).start.j && start.j == (*it).end.j)
				if( (start.i + end.i) == ((*it).start.i + (*it).end.i) &&
					(start.k + end.k) == ((*it).start.k + (*it).end.k))
					return true;
		return false;
	}
	
	if(type == xyz)
	{
		for(it = LL.begin() ; it != LL.end() ; it++)
			if( (start.i + end.i) == ((*it).start.i + (*it).end.i) &&
				(start.j + end.j) == ((*it).start.j + (*it).end.j) &&
				(start.k + end.k) == ((*it).start.k + (*it).end.k) )
				return true;
		return false;
	}
	return false;
}

/*
bool Link::isCrossing(ListLink& LL) const
{
	ListLink::iterator it;
	for(it = LL.begin() ; it != LL.end() ; it++)
	{
*/
		/******************************************************************************/
		/* Note: The two next conditions are an introduction compared to V. Pasko's	  */
		/* code. However, they are pointless since link between 2 points already	  */
		/* linked is forbidden														  */
		/******************************************************************************/
/*		
		if (start == (*it).start && end == (*it).end)
			return true;								// twice the same link
		if (start == (*it).end && end == (*it).start)
			return true;								// twice the same link
*/		
		/******************************************************************************/
		/* Note: With the correction above, when a 1-D, 2-D, or 3-D link is tested	  */
		/* with itself the value return is 1. In Pasko et al. 2000, this was only true*/
		/* 2-D, 3-D links															  */
		/******************************************************************************/
/*		
		if(type == xy && (*it).type == xy && start.k == (*it).start.k)
		{
			if( (start.i + end.i) == ((*it).start.i + (*it).end.i) &&
				(start.j + end.j) == ((*it).start.j + (*it).end.j))
				return true;
		}												// 2 links in same xy plane
		
		if(type == xz && (*it).type == xz && start.j == (*it).start.j)
		{
			if( (start.i + end.i) == ((*it).start.i + (*it).end.i) &&
				(start.k + end.k) == ((*it).start.k + (*it).end.k))
				return true;
		}
		
		if(type == yz && (*it).type == yz && start.i == (*it).start.i)
		{
			if( (start.j + end.j) == ((*it).start.j + (*it).end.j) &&
				(start.k + end.k) == ((*it).start.k + (*it).end.k))
				return true;
		}
		
		if(type == xyz && (*it).type == xyz)
		{
			if( (start.i + end.i) == ((*it).start.i + (*it).end.i) &&
				(start.j + end.j) == ((*it).start.j + (*it).end.j) &&
				(start.k + end.k) == ((*it).start.k + (*it).end.k) )
				return true;							// return true if the 2 links are identical
		}
	}
	return false;
}
*/
ostream & operator<< (ostream & os, const Link & LL)
{
	return os<<"\nLink parameters:"
			<<"\n Starting Point = ["<<LL.start.i<<" "<<LL.start.j<<" "<<LL.start.k
			<<"]\n Ending Point = ["<<LL.end.i<<" "<<LL.end.j<<" "<<LL.end.k
			<<"]\n l = "<<LL.l
			<<"\n Potential drop = "<<LL.deltaV
			<<"\n E-field = "<<LL.efield<<endl;
}
/**************************************************************************************/