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

#include <sstream>
#include <math.h>
#include "RbtConstraint.h"
#include "RbtLigandError.h"

using std::ostringstream;
using std::ends;

// initialization of the static data of RbtConstraint
RbtString RbtConstraint::_CT("RbtConstraint");
RbtInt RbtHeavyConstraint::counter = 0;//7 Feb 2005 (DM) new constraint type
RbtInt RbtHBAConstraint::counter = 0;
RbtInt RbtHBDConstraint::counter = 0;
RbtInt RbtHydroConstraint::counter = 0;
RbtInt RbtHydroAliphaticConstraint::counter = 0;
RbtInt RbtHydroAromaticConstraint::counter = 0;
RbtInt RbtNegChargeConstraint::counter = 0;
RbtInt RbtPosChargeConstraint::counter = 0;
RbtInt RbtRingAromaticConstraint::counter = 0;
//June 2014 (DAG) new constraint types
//RbtInt RbtHeavyConstraintEXP::counter = 0;
//RbtInt RbtHBAConstraintEXP::counter = 0;
//RbtInt RbtHBDConstraintEXP::counter = 0;
//RbtInt RbtHydroConstraintEXP::counter = 0;
//RbtInt RbtHydroAliphaticConstraintEXP::counter = 0;
//RbtInt RbtHydroAromaticConstraintEXP::counter = 0;
//RbtInt RbtNegChargeConstraintEXP::counter = 0;
//RbtInt RbtPosChargeConstraintEXP::counter = 0;
//RbtInt RbtRingAromaticConstraintEXP::counter = 0;


//DAG added a and b to all constructors
RbtConstraint::RbtConstraint(RbtCoord c, RbtDouble t)
{
    m_atomList.clear();
    coord = c;
    tolerance = t;
    a = 0.0;
    b = 0.0;
    isexp = false;
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtConstraint::RbtConstraint(RbtCoord c, RbtDouble t, RbtDouble _a, RbtDouble _b)
{
    m_atomList.clear();
    coord = c;
    tolerance = t;
    a = _a;
    b = _b;
    isexp = true;
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}

RbtConstraint::RbtConstraint(const RbtConstraint& c)
{
    coord = c.coord;
    tolerance = c.tolerance;
    a = 0.0;
    b = 0.0;
    isexp = false;
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}
 
RbtConstraint::RbtConstraint()
{
    coord = RbtCoord(0,0,0);
    tolerance = 0.0;
    a = 0.0;
    b = 0.0;
    isexp = false;
    _RBTOBJECTCOUNTER_CONSTR_(_CT);
}
 
RbtConstraint::~RbtConstraint()
{
    _RBTOBJECTCOUNTER_DESTR_(_CT);
}

// DAG
void RbtConstraint::copy(const RbtConstraint & c)
{
    coord = c.coord;
    tolerance = c.tolerance;
    a = c.a;
    b = c.b;
    isexp = c.isexp;
}

ostream& RbtConstraint::Print(ostream& o) const
{
    o << coord << "\t" << tolerance << "\t" << a << "\t" << b;
    return o;
}

ostream& operator<<(ostream& s, const RbtConstraint &c)
{
    return (c.Print(s));
}

//DM 17 Oct 2003 - rewrite to avoid step increase in score at the tolerance distance
//Allow for fact that atomlist may be empty
RbtDouble RbtConstraint::Score() const
{
  if (m_atomList.empty()) {
    //Return artifically high score if the ligand has none of the desired features
    //This will only apply to the optional constraints (mandatory constraints are prefiltered)
    //SRC add check of isexp to avoid 99 penalty when restraint not present, from Juan
    if (isexp) {
	return 4.18*a;
    } else {
    	return 99.0;
    }
  }

  //Compile a list of all the distance**2 to the constraint center
  RbtDoubleList d2list;
  for (RbtAtomListConstIter iter = m_atomList.begin(); iter != m_atomList.end() ; iter++) {
    d2list.push_back(Rbt::Length2((*iter)->GetCoords(), coord));
  }
  //Min distance**2
  RbtDouble d2min = *(std::min_element(d2list.begin(),d2list.end()));

  // DAG return score differentially
  // if a and b are != 0.0 return exponential function score
  // else, return quadratic score
  RbtDouble dr = sqrt(d2min) - tolerance;
  
  if (isexp){
      RbtDouble r;
      r = (dr > 0.0) ? 4.18*a*(1-exp(-pow(dr,2)/b)) : 0.0;
      //if (r != 0.0) cout << "Scoring with exp function: " << r << " dr: " << dr << endl;
      return r;
  }else{ //Score is zero up to the tolerance range, then quadratic beyond.
  return (dr > 0.0) ? dr*dr : 0.0;
  }
}


RbtConstraintPtr Rbt::CreateConstraint(RbtCoord& c, RbtDouble& t, RbtString& n, RbtBool bCount)
{
  if (n == "Any") {
    if (bCount) RbtHeavyConstraint::counter++;
    return new RbtHeavyConstraint(c,t);
  }
  else if (n == "Don") {
    if (bCount) RbtHBDConstraint::counter++;
    return new RbtHBDConstraint(c,t);
  }
  else if (n == "Acc") {
    if (bCount) RbtHBAConstraint::counter++;
    return new RbtHBAConstraint(c,t);
  }
  else if (n == "Hyd") {
    if (bCount) RbtHydroConstraint::counter++;
    return new RbtHydroConstraint(c,t);
  }
  else if (n == "Ani") {
    if (bCount) RbtNegChargeConstraint::counter++;
    return new RbtNegChargeConstraint(c,t);
  }
  else if (n == "Cat") {
    if (bCount) RbtPosChargeConstraint::counter++;
    return new RbtPosChargeConstraint(c,t);
  }
  else if (n == "Aro") {
    if (bCount) RbtRingAromaticConstraint::counter++;
    return new RbtRingAromaticConstraint(c,t);
  }
  else if (n == "Hal") {
    if (bCount) RbtHydroAliphaticConstraint::counter++;
    return new RbtHydroAliphaticConstraint(c,t);
  }
  else if (n == "Har") {
    if (bCount) RbtHydroAromaticConstraint::counter++;
    return new RbtHydroAromaticConstraint(c,t);
  }
  else
    throw RbtError(_WHERE_,"Constraint " + n + " not recognized");
}
RbtConstraintPtr Rbt::CreateConstraint(RbtCoord& c, RbtDouble& t, RbtDouble& a, RbtDouble& b, RbtString& n, RbtBool bCount)
{
    cout << "Creating ph4 for " << n << "type EXP" << endl;
  if (n == "Any") {
    if (bCount) RbtHeavyConstraint::counter++;
    return new RbtHeavyConstraint(c,t,a,b);
  }
  else if (n == "Don") {
    if (bCount) RbtHBDConstraint::counter++;
    return new RbtHBDConstraint(c,t,a,b);
  }
  else if (n == "Acc") {
    if (bCount) RbtHBAConstraint::counter++;
    return new RbtHBAConstraint(c,t,a,b);
  }
  else if (n == "Hyd") {
    if (bCount) RbtHydroConstraint::counter++;
    return new RbtHydroConstraint(c,t,a,b);
  }
  else if (n == "Ani") {
    if (bCount) RbtNegChargeConstraint::counter++;
    return new RbtNegChargeConstraint(c,t,a,b);
  }
  else if (n == "Cat") {
    if (bCount) RbtPosChargeConstraint::counter++;
    return new RbtPosChargeConstraint(c,t,a,b);
  }
  else if (n == "Aro") {
    if (bCount) RbtRingAromaticConstraint::counter++;
    return new RbtRingAromaticConstraint(c,t,a,b);
  }
  else if (n == "Hal") {
    if (bCount) RbtHydroAliphaticConstraint::counter++;
    return new RbtHydroAliphaticConstraint(c,t,a,b);
  }
  else if (n == "Har") {
    if (bCount) RbtHydroAromaticConstraint::counter++;
    return new RbtHydroAromaticConstraint(c,t,a,b);
  }
  else
    throw RbtError(_WHERE_,"Constraint " + n + " not recognized");
}

//RbtConstraintPtr Rbt::CreateConstraintEXP(RbtCoord& c, RbtDouble& t, RbtDouble& a, RbtDouble& b, RbtString& n, RbtBool bCount)
//{
//  if (n == "Any") {
//    if (bCount) RbtHeavyConstraintEXP::counter++;
//    return new RbtHeavyConstraintEXP(c,t,a,b);
//  }
//  else if (n == "Don") {
//    if (bCount) RbtHBDConstraintEXP::counter++;
//    return new RbtHBDConstraintEXP(c,t,a,b);
//  }
//  else if (n == "Acc") {
//    if (bCount) RbtHBAConstraintEXP::counter++;
//    return new RbtHBAConstraintEXP(c,t,a,b);
//  }
//  else if (n == "Hyd") {
//    if (bCount) RbtHydroConstraintEXP::counter++;
//    return new RbtHydroConstraintEXP(c,t,a,b);
//  }
//  else if (n == "Ani") {
//    if (bCount) RbtNegChargeConstraintEXP::counter++;
//    return new RbtNegChargeConstraintEXP(c,t,a,b);
//  }
//  else if (n == "Cat") {
//    if (bCount) RbtPosChargeConstraintEXP::counter++;
//    return new RbtPosChargeConstraintEXP(c,t,a,b);
//  }
//  else if (n == "Aro") {
//    if (bCount) RbtRingAromaticConstraintEXP::counter++;
//    return new RbtRingAromaticConstraintEXP(c,t,a,b);
//  }
//  else if (n == "Hal") {
//    if (bCount) RbtHydroAliphaticConstraintEXP::counter++;
//    return new RbtHydroAliphaticConstraintEXP(c,t,a,b);
//  }
//  else if (n == "Har") {
//    if (bCount) RbtHydroAromaticConstraintEXP::counter++;
//    return new RbtHydroAromaticConstraintEXP(c,t,a,b);
//  }
//  else
//    throw RbtError(_WHERE_,"Constraint " + n + " not recognized");
//}



void Rbt::ZeroCounters()
{
  RbtHeavyConstraint::counter = 0;
  RbtHBAConstraint::counter = 0;
  RbtHBDConstraint::counter = 0;
  RbtHydroConstraint::counter = 0;
  RbtHydroAliphaticConstraint::counter = 0;
  RbtHydroAromaticConstraint::counter = 0;
  RbtNegChargeConstraint::counter = 0;
  RbtPosChargeConstraint::counter = 0;
  RbtRingAromaticConstraint::counter = 0;
//  RbtHeavyConstraintEXP::counter = 0;
//  RbtHBAConstraintEXP::counter = 0;
//  RbtHBDConstraintEXP::counter = 0;
//  RbtHydroConstraintEXP::counter = 0;
//  RbtHydroAliphaticConstraintEXP::counter = 0;
//  RbtHydroAromaticConstraintEXP::counter = 0;
//  RbtNegChargeConstraintEXP::counter = 0;
//  RbtPosChargeConstraintEXP::counter = 0;
//  RbtRingAromaticConstraintEXP::counter = 0;
}

void Rbt::ReadConstraint(istream& ifile, RbtConstraintPtr& cnt, RbtBool bCount)
{
  RbtString n;
  RbtCoord c;
  RbtDouble t;
  ifile >> c;
  ifile >> t;
  ifile >> n;
  cout << c << "\t" << t << "\t" << n << endl;
  if (ifile.fail())
    throw RbtError(_WHERE_,"Problems reading constraint " + n );
  cnt = Rbt::CreateConstraint(c,t,n,bCount);
}

void Rbt::ReadConstraintEXP(istream& ifile, RbtConstraintPtr& cnt, RbtBool bCount)
{
  RbtString n;
  RbtCoord c;
  RbtDouble t;
  RbtDouble a;
  RbtDouble b;
  ifile >> c;
  ifile >> t;
  ifile >> a;
  ifile >> b;
  ifile >> n;
  cout << c << "\t" << t << "\t" << a << "\t" << b << "\t" << n << endl;
  if (ifile.fail())
    throw RbtError(_WHERE_,"Problems reading constraint " + n );
  cnt = Rbt::CreateConstraint(c,t,a,b,n,bCount);
}

void Rbt::ReadConstraintFromMoe(istream& ifile, RbtConstraintPtr& cnt, RbtBool bCount)
{
    string n, extra;
    RbtDouble x, y, z, t;
    ifile >> n;
    cout << n << endl;
    ifile >> extra;
    ifile >> x >> y >> z >> t;
    ifile >> extra >> extra;
    if (ifile.fail())
        throw RbtError(_WHERE_,"Problems reading constraint " + n );
    RbtCoord c(x,y,z);
    cout << c << "\t" << t << "\t" << n << endl;
    cnt = Rbt::CreateConstraint(c,t,n,bCount);
}

void Rbt::ReadStartMoe(istream& ifile)
{
    const int MAXLINELENGTH = 100; // I'm quite sure the max length for any 
                                   // line in moe is 72, but
                                   // prefer to be on the safe side 
    char* line = new char[MAXLINELENGTH+1];
    RbtConstraintPtr cons;
    char ch;
    ifile.getline(line, MAXLINELENGTH);  // get rid of first four lines
    ifile.getline(line, MAXLINELENGTH);
    ifile.getline(line, MAXLINELENGTH);  
        // This line has some data that could be neccesary later on, like
        // partial match information also, it contains the
        // values for the fields comment and title, that can be as
        // long as necessary. Need to make sure I read all the lines
        // before the next commentted line (starting with a #)
    while ((ch = ifile.peek()) != '#')
        ifile.getline(line, MAXLINELENGTH);
    ifile.getline(line, MAXLINELENGTH);  // next comment line
}

void Rbt::ReadConstraints(istream& ifile, RbtConstraintList& cl, RbtBool bCount)
{
    char ch;
    RbtBool moe = false;
    RbtBool exp = false;
    cl.clear();
    RbtConstraintPtr cons;
    ch = ifile.peek();
    if (ch == '#')  // moe file
    {
        moe = true;
        ReadStartMoe(ifile);
    }else if (ch == '!')
    {
	exp = true;
        char* linescape = new char[100+1];
        ifile.getline(linescape, 100);
    }
    while ((!moe && !ifile.eof()) || 
           (moe  && ((ch = ifile.peek()) != '#')))
    {
        if (moe){
            ReadConstraintFromMoe(ifile, cons, bCount);
	}else {
            if (exp){
                cout << "Reading EXP constraint file" << endl;
                ReadConstraintEXP(ifile, cons, bCount);
            }else{
                ReadConstraint(ifile, cons, bCount);
            }
        }
        
        cl.push_back(cons);
        while (ifile.get(ch))
        {
            if (!isspace(ch))
            {
                ifile.putback(ch);
                break;
            }
        }
    }
    if (cl.size() == 0)
        cout << RbtConstraint::_CT 
             << "** WARNING: number of constraints read is 0\n";
}

//07 Feb 2005 (DM) - new constraint type, any heavy atom
void RbtHeavyConstraint::AddAtomList(RbtModelPtr lig, RbtBool bCheck) throw (RbtError)
{
  m_atomList = Rbt::GetAtomList(lig->GetAtomList(), std::not1(Rbt::isAtomicNo_eq(1)));
  if (bCheck && (m_atomList.size() < counter) && !isexp) {
    ostringstream ostr;
    ostr << "The ligand has only " << m_atomList.size()
         << " heavy atom(s) (" << counter << " required)" << ends;
    throw RbtLigandError(_WHERE_, ostr.str());
  }
}

//16 May 2003 (DM) - limit to neutral H-bond donor hydrogens
void RbtHBDConstraint::AddAtomList(RbtModelPtr lig, RbtBool bCheck) throw (RbtError)
{
  m_atomList = Rbt::GetAtomList(lig->GetAtomList(), Rbt::isAtomHBondDonor());
  m_atomList = Rbt::GetAtomList(m_atomList, std::not1(Rbt::isAtomCationic()));
  if (bCheck && (m_atomList.size() < counter) && !isexp) {
    ostringstream ostr;
    ostr << "The ligand has only " << m_atomList.size() 
	 << " neutral H-bond donor hydrogens(s) (" << counter << " required)" << ends;
    throw RbtLigandError(_WHERE_, ostr.str());
  }
}

//16 May 2003 (DM) - limit to neutral H-bond acceptors
void RbtHBAConstraint::AddAtomList(RbtModelPtr lig, RbtBool bCheck) throw (RbtError)
{
  m_atomList = Rbt::GetAtomList(lig->GetAtomList(), Rbt::isAtomHBondAcceptor());
  m_atomList = Rbt::GetAtomList(m_atomList, std::not1(Rbt::isAtomAnionic()));
  if (bCheck && (m_atomList.size() < counter) && !isexp) {
    ostringstream ostr;
    ostr << "The ligand has only " << m_atomList.size() 
	 << " neutral H-bond acceptor(s) (" << counter << " required)" << ends;
    throw RbtLigandError(_WHERE_, ostr.str());
  }
}

void RbtHydroConstraint::AddAtomList(RbtModelPtr lig, RbtBool bCheck) throw (RbtError)
{
  m_atomList = Rbt::GetAtomList(lig->GetAtomList(), Rbt::isAtomLipophilic());
  if (bCheck && (m_atomList.size() < counter) && !isexp) {
    ostringstream ostr;
    ostr << "The ligand has only " << m_atomList.size() 
	 << " hydrophobic atom(s) (" << counter << " required)" << ends;
    throw RbtLigandError(_WHERE_, ostr.str());
  }
}

void RbtHydroAliphaticConstraint::AddAtomList(RbtModelPtr lig, RbtBool bCheck) throw (RbtError)
{
  m_atomList = Rbt::GetAtomList(lig->GetAtomList(), Rbt::isAtomLipophilic());
  m_atomList = Rbt::GetAtomList(m_atomList, 
				Rbt::isHybridState_eq(RbtAtom::SP3));
  if (bCheck && (m_atomList.size() < counter) && !isexp) {
    ostringstream ostr;
    ostr << "The ligand has only " << m_atomList.size() 
	 << " hydrophobic aliphatic atom(s) (" << counter << " required)" << ends;
    throw RbtLigandError(_WHERE_, ostr.str());
  }
}

void RbtHydroAromaticConstraint::AddAtomList(RbtModelPtr lig, RbtBool bCheck) throw (RbtError)
{
  m_atomList = Rbt::GetAtomList(lig->GetAtomList(), Rbt::isAtomLipophilic());
  m_atomList = Rbt::GetAtomList(m_atomList, 
				Rbt::isHybridState_eq(RbtAtom::AROM));
  if (bCheck && (m_atomList.size() < counter) && !isexp) {
    ostringstream ostr;
    ostr << "The ligand has only " << m_atomList.size() 
	 << " hydrophobic aromatic atom(s) (" << counter << " required)" << ends;
    throw RbtLigandError(_WHERE_, ostr.str());
  }
}

void RbtNegChargeConstraint::AddAtomList(RbtModelPtr lig, RbtBool bCheck) throw (RbtError)
{
  m_atomList = Rbt::GetAtomList(lig->GetAtomList(), Rbt::isAtomAnionic());
  if (bCheck && (m_atomList.size() < counter) && !isexp) {
    ostringstream ostr;
    ostr << "The ligand has only " << m_atomList.size() 
	 << " negatively charged atom(s) (" << counter << " required)" << ends;
    throw RbtLigandError(_WHERE_, ostr.str());
  }
}

void RbtPosChargeConstraint::AddAtomList(RbtModelPtr lig, RbtBool bCheck) throw (RbtError)
{
  m_atomList = Rbt::GetAtomList(lig->GetAtomList(), Rbt::isAtomCationic());
  if (bCheck && (m_atomList.size() < counter) && !isexp) {
    ostringstream ostr;
    ostr << "The ligand has only " << m_atomList.size() 
	 << " positively charged atom(s) (" << counter << " required)" << ends;
    throw RbtLigandError(_WHERE_, ostr.str());
  }
}


void RbtRingAromaticConstraint::AddAtomList(RbtModelPtr lig, RbtBool bCheck) throw (RbtError)
{
    m_atomList.clear();
    RbtAtomList at = lig->GetAtomList();
    RbtAtomListList ligRingLists = lig->GetRingAtomLists();
    for (RbtAtomListListConstIter rIter = ligRingLists.begin() ;
         rIter != ligRingLists.end(); rIter++)
    {
        if (Rbt::GetNumAtoms(*rIter,Rbt::isPiAtom()) == 
                 (*rIter).size())
        {
            RbtPseudoAtomPtr spPseudoAtom = 
                   lig->AddPseudoAtom(*rIter);
            m_atomList.push_back(spPseudoAtom);
        }
    }
    if (bCheck && (m_atomList.size() < counter) && !isexp)
      {
	ostringstream ostr;
	ostr << "The ligand has only " << m_atomList.size() 
	     << " aromatic ring(s) (" << counter << " required)" << ends;
        throw RbtLigandError(_WHERE_, ostr.str());
      }
}
