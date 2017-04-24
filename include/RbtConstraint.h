/***********************************************************************
* The rDock program was developed from 1998 - 2006 by the software team 
* at RiboTargets (subsequently Vernalis (R&D) Ltd).
* In 2006, the software was licensed to the University of York for 
* maintenance and distribution.
* In 2012, Vernalis and the University of York agreed to release the 
* program as Open Source software.
* This version is licensed under GNU-LGPL version 3.0 with support from
* the University of Barcelona.
* http://rdock.sourceforge.net/
***********************************************************************/

//Base class for pharmacophore constraints

#ifndef _RBT_CONSTRAINT_H_
#define _RBT_CONSTRAINT_H_

#include "RbtTypes.h"
#include "RbtCoord.h"
#include "RbtAtom.h"
#include "RbtModel.h"
#include "RbtError.h"

class RbtConstraint 
{
public:
  static RbtString _CT;
    ///////////////////
    // Constructors
    ///////////////////
  RbtConstraint(RbtCoord, RbtDouble);
  RbtConstraint(RbtCoord, RbtDouble, RbtDouble, RbtDouble);
  RbtConstraint(const RbtConstraint&);
  RbtConstraint();
 
    ///////////////////
    // Destructor
    //////////////////
  virtual ~RbtConstraint();
  
  void copy(const RbtConstraint &);
  ostream& Print(ostream&) const;
  friend ostream& operator<<(ostream&, const RbtConstraint &);

  RbtCoord GetCoords() const { return coord; };
  RbtDouble GetTolerance() const { return tolerance; };
  RbtDouble GetAParm() const { return a;};
  RbtDouble GetBParm() const { return b;};
  virtual void AddAtomList(RbtModelPtr, RbtBool bCheck=true) = 0;
  RbtDouble Score() const;

    ////////////////////
    // Private methods
    ////////////////////
protected:
  RbtCoord coord;
  RbtDouble tolerance;
  RbtDouble a;
  RbtDouble b;
  RbtBool isexp;
  RbtAtomList m_atomList;
  
private:
    /////////////////////
    // Private data
    /////////////////////  
};

//Useful typedefs
typedef SmartPtr<RbtConstraint> RbtConstraintPtr;  //Smart pointer
typedef vector<RbtConstraintPtr> RbtConstraintList;//Vector of smart pointers
typedef RbtConstraintList::iterator RbtConstraintListIter;
typedef RbtConstraintList::const_iterator RbtConstraintListConstIter;

namespace Rbt {
  RbtConstraintPtr CreateConstraint(RbtCoord& c, RbtDouble& t, RbtString& n, RbtBool bCount=true);
  RbtConstraintPtr CreateConstraint(RbtCoord& c, RbtDouble& t, RbtDouble& a, RbtDouble& b, RbtString& n, RbtBool bCount=true);

//  RbtConstraintPtr CreateConstraintEXP(RbtCoord& c, RbtDouble& t, RbtDouble& a, RbtDouble& b, RbtString& n, RbtBool bCount=true);
  void ZeroCounters();
  void ReadConstraint(istream& ifile, RbtConstraintPtr& cnt, RbtBool bCount=true);
  void ReadConstraintEXP(istream& ifile, RbtConstraintPtr& cnt, RbtBool bCount);
  void ReadConstraintFromMoe(istream& ifile, RbtConstraintPtr& cnt, RbtBool bCount=true);
  void ReadStartMoe(istream& ifile);
  void ReadConstraints(istream& ii, RbtConstraintList& cl, RbtBool bCount=true);
}

//7 Feb 2005 (DM, Enspiral Discovery)
//New constraint type - any heavy atom
//Likely to be inefficient for large ligands, but has the advantage
//of being simple to code and run
class RbtHeavyConstraint : public RbtConstraint
{
    public:
  RbtHeavyConstraint(RbtCoord c, RbtDouble t) : RbtConstraint(c,t) {};
  RbtHeavyConstraint(RbtCoord c, RbtDouble t, RbtDouble a, RbtDouble b) : RbtConstraint(c,t,a,b) {};
  void AddAtomList(RbtModelPtr, RbtBool bCheck=true) throw(RbtError);
  static RbtInt counter;
};

class RbtHBAConstraint : public RbtConstraint
{
    public:
  RbtHBAConstraint(RbtCoord c, RbtDouble t) : RbtConstraint(c,t) {};
  RbtHBAConstraint(RbtCoord c, RbtDouble t, RbtDouble a, RbtDouble b) : RbtConstraint(c,t,a,b) {};
  void AddAtomList(RbtModelPtr, RbtBool bCheck=true) throw(RbtError);
  static RbtInt counter;
};

class RbtHBDConstraint : public RbtConstraint
{
    public:
  RbtHBDConstraint(RbtCoord c, RbtDouble t) : RbtConstraint(c,t){};
  RbtHBDConstraint(RbtCoord c, RbtDouble t, RbtDouble a, RbtDouble b) : RbtConstraint(c,t,a,b) {};
  void AddAtomList(RbtModelPtr, RbtBool bCheck=true) throw(RbtError);
  static RbtInt counter;
};

class RbtHydroConstraint : public RbtConstraint
{
    public:
  RbtHydroConstraint(RbtCoord c, RbtDouble t) : RbtConstraint(c,t){};
  RbtHydroConstraint(RbtCoord c, RbtDouble t, RbtDouble a, RbtDouble b) : RbtConstraint(c,t,a,b) {};
  void AddAtomList(RbtModelPtr, RbtBool bCheck=true) throw(RbtError);
  static RbtInt counter;
};

class RbtHydroAliphaticConstraint : public RbtConstraint
{
    public:
  RbtHydroAliphaticConstraint(RbtCoord c, RbtDouble t) : RbtConstraint(c,t){};
  RbtHydroAliphaticConstraint(RbtCoord c, RbtDouble t, RbtDouble a, RbtDouble b) : RbtConstraint(c,t,a,b) {};
  void AddAtomList(RbtModelPtr, RbtBool bCheck=true) throw(RbtError);
  static RbtInt counter;
};

class RbtHydroAromaticConstraint : public RbtConstraint
{
    public:
  RbtHydroAromaticConstraint(RbtCoord c, RbtDouble t) : RbtConstraint(c,t){};
  RbtHydroAromaticConstraint(RbtCoord c, RbtDouble t, RbtDouble a, RbtDouble b) : RbtConstraint(c,t,a,b) {};
  void AddAtomList(RbtModelPtr, RbtBool bCheck=true) throw(RbtError);
  static RbtInt counter;
};


class RbtNegChargeConstraint : public RbtConstraint
{
    public:
  RbtNegChargeConstraint(RbtCoord c, RbtDouble t) : RbtConstraint(c,t){};
  RbtNegChargeConstraint(RbtCoord c, RbtDouble t, RbtDouble a, RbtDouble b) : RbtConstraint(c,t,a,b) {};
  void AddAtomList(RbtModelPtr, RbtBool bCheck=true) throw(RbtError);
  static RbtInt counter;
};

class RbtPosChargeConstraint : public RbtConstraint
{
    public:
  RbtPosChargeConstraint(RbtCoord c, RbtDouble t) : RbtConstraint(c,t){};
  RbtPosChargeConstraint(RbtCoord c, RbtDouble t, RbtDouble a, RbtDouble b) : RbtConstraint(c,t,a,b) {};
  void AddAtomList(RbtModelPtr, RbtBool bCheck=true) throw(RbtError);
  static RbtInt counter;
};

class RbtRingAromaticConstraint : public RbtConstraint
{
    public:
  RbtRingAromaticConstraint(RbtCoord c, RbtDouble t) : RbtConstraint(c,t){};
  RbtRingAromaticConstraint(RbtCoord c, RbtDouble t, RbtDouble a, RbtDouble b) : RbtConstraint(c,t,a,b) {};
  void AddAtomList(RbtModelPtr, RbtBool bCheck=true) throw(RbtError);
  static RbtInt counter;
};
//
//// DAG Types
//class RbtRingAromaticConstraintEXP : public RbtConstraint
//{
//    public:
//  RbtRingAromaticConstraintEXP(RbtCoord c, RbtDouble t, RbtDouble a, RbtDouble b) : RbtConstraint(c,t,a,b){};
//  void AddAtomList(RbtModelPtr, RbtBool bCheck=true) throw(RbtError);
//  static RbtInt counter;
//};
//
//class RbtPosChargeConstraintEXP : public RbtConstraint
//{
//    public:
//  RbtPosChargeConstraintEXP(RbtCoord c, RbtDouble t, RbtDouble a, RbtDouble b) : RbtConstraint(c,t,a,b){};
//  void AddAtomList(RbtModelPtr, RbtBool bCheck=true) throw(RbtError);
//  static RbtInt counter;
//};
//
//class RbtNegChargeConstraintEXP : public RbtConstraint
//{
//    public:
//  RbtNegChargeConstraintEXP(RbtCoord c, RbtDouble t, RbtDouble a, RbtDouble b) : RbtConstraint(c,t,a,b){};
//  void AddAtomList(RbtModelPtr, RbtBool bCheck=true) throw(RbtError);
//  static RbtInt counter;
//};
//
//class RbtHydroAromaticConstraintEXP : public RbtConstraint
//{
//    public:
//  RbtHydroAromaticConstraintEXP(RbtCoord c, RbtDouble t, RbtDouble a, RbtDouble b) : RbtConstraint(c,t,a,b){};
//  void AddAtomList(RbtModelPtr, RbtBool bCheck=true) throw(RbtError);
//  static RbtInt counter;
//};
//
//class RbtHydroAliphaticConstraintEXP : public RbtConstraint
//{
//    public:
//  RbtHydroAliphaticConstraintEXP(RbtCoord c, RbtDouble t, RbtDouble a, RbtDouble b) : RbtConstraint(c,t,a,b){};
//  void AddAtomList(RbtModelPtr, RbtBool bCheck=true) throw(RbtError);
//  static RbtInt counter;
//};
//
//class RbtHydroConstraintEXP : public RbtConstraint
//{
//    public:
//  RbtHydroConstraintEXP(RbtCoord c, RbtDouble t, RbtDouble a, RbtDouble b) : RbtConstraint(c,t,a,b){};
//  void AddAtomList(RbtModelPtr, RbtBool bCheck=true) throw(RbtError);
//  static RbtInt counter;
//};
//
//class RbtHBAConstraintEXP : public RbtConstraint
//{
//    public:
//  RbtHBAConstraintEXP(RbtCoord c, RbtDouble t, RbtDouble a, RbtDouble b) : RbtConstraint(c,t,a,b){};
//  void AddAtomList(RbtModelPtr, RbtBool bCheck=true) throw(RbtError);
//  static RbtInt counter;
//};
//
//class RbtHBDConstraintEXP : public RbtConstraint
//{
//    public:
//  RbtHBDConstraintEXP(RbtCoord c, RbtDouble t, RbtDouble a, RbtDouble b) : RbtConstraint(c,t,a,b){};
//  void AddAtomList(RbtModelPtr, RbtBool bCheck=true) throw(RbtError);
//  static RbtInt counter;
//};
//
//class RbtHeavyConstraintEXP : public RbtConstraint
//{
//    public:
//  RbtHeavyConstraintEXP(RbtCoord c, RbtDouble t, RbtDouble a, RbtDouble b) : RbtConstraint(c,t,a,b){};
//  void AddAtomList(RbtModelPtr, RbtBool bCheck=true) throw(RbtError);
//  static RbtInt counter;
//};
#endif //_RbtConstraint_H_
