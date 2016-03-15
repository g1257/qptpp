/*
Copyright (c) 2015, UT-Battelle, LLC

QuantumPerturbation, Version 0.1

This file is part of QuantumPerturbation.
QuantumPerturbation is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
QuantumPerturbation is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with QuantumPerturbation. If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef BRACKET_H
#define BRACKET_H

#include "EquationTerm.h"

namespace QuantumPerturbation {

template<typename RealType>
class Bracket {

	typedef Bracket<RealType> ThisType;

public:

	typedef EquationTerm<RealType> EquationTermType;
	typedef typename EquationTermType::StateType StateType;
	typedef typename EquationTermType::EnergyType EnergyType;


	Bracket(const EquationTermType& ket1,const EquationTermType& bra1)
	    : ket(ket1),
	      bra(bra1),
	      factor_(ket.factor() * bra.factor()),
	      a_(EnergyType('p',-1)),
	      op_(EquationTermType::OP_IDENTITY),
	      energy_(""),
	      bracket_("")
	{
		setEnergy();
		setBracket();
	}

	void operator+=(const ThisType& other)
	{
		assert(sameType(other));
		assert(energy_ == other.energy_);
		factor_ += other.factor();
	}

	RealType factor() const
	{
		return factor_;
	}

	PsimagLite::String a() const
	{
		if (a_.second <= 0) return "";
		PsimagLite::String str("a");
		str += a_.first;
		str += "(" + ttos(a_.second) + ")";
		return str;
	}

	const PsimagLite::String& energy() const
	{
		return energy_;
	}

	const PsimagLite::String& op() const
	{
		return bracket_;
	}

	bool sameType(const ThisType& other) const
	{
		return (bracket_ == other.bracket_ &&
		        energy_ == other.energy_ &&
		        op_ == other.op_);
	}

	bool isNull() const
	{
		if (bracket_ == "<m(0)|n(0)> ") return true;

		return false;
	}

	void replace()
	{
		if (ket.stateChar() != 'n' || a_.second >= 0)
			return;

		if (ket.stateIndex()!= 0) {
			a_ = StateType('p',ket.stateIndex());
			ket.state(StateType('p',0));
			setBracket();
		}

		if (bra.stateIndex() !=0) return;

		if (op_ == EquationTermType::OP_IDENTITY) {
			bracket_ = "";
			a_.first = bra.stateChar();
		} else if (op_ == EquationTermType::OP_V) {
			bracket_ = "V(" + ttos(bra.stateChar());
			bracket_ += "," + ttos(ket.stateChar()) + ")";
		}
	}

private:

	void setEnergy()
	{
		energy_= ket.energyString() + bra.energyString();
	}

	void setBracket()
	{
		op_ =  bra.op() + ket.op();
		switch (op_) {
		case 0:
			setBracketIdentity();
			break;
		case 1:
			setBracketH0();
			break;
		case 2:
			setBracketV();
			break;
		default:
			bracket_="<" + bra.stateString() +
			        "|" + opToString(op_) +
			        "|" + ket.stateString() +
			        ">";
		}
		bracket_ += " ";
	}

	PsimagLite::String opToString(SizeType op) const
	{
		switch(op) {
		case 0:
			assert(false);
			return "";
		case 1:
			return "H0";
		case 2:
			return "V";
		default:
			return "UNIMPLEMENTED";
		}
	}

	void setBracketIdentity()
	{
		if (bra.stateChar() != ket.stateChar()) {
			bracket_="<" + bra.stateString() + "|" + ket.stateString() + ">";
			return;
		}

		if (bra.stateIndex() == 0 && ket.stateIndex() == 0) {
			return;
		}
		if (bra.stateIndex() == 0 && ket.stateIndex() == 1) {
			factor_=0;
			return;
		}
		if (bra.stateIndex() == 1 && ket.stateIndex() == 0) {
			factor_=0;
			return;
		}
		bracket_="<" + bra.stateString() + "|" + ket.stateString() + ">";
	}

	void setBracketH0()
	{
		if (ket.stateIndex() == 0) {
			EnergyType energy1(ket.stateChar(),0);
			StateType state1(ket.stateChar(),0);
			EquationTermType dummy(EquationTermType::OP_IDENTITY,state1,energy1,0);
			energy_ += dummy.energyString();
			assert(op_ >= EquationTermType::OP_H0);
			op_ -= EquationTermType::OP_H0;
			setBracketIdentity();
			return;
		}

		if (bra.stateIndex() == 0) {
			EnergyType energy1(ket.stateChar(),0);
			StateType state1(ket.stateChar(),0);
			EquationTermType dummy(EquationTermType::OP_IDENTITY,state1,energy1,0);
			energy_ += dummy.energyString();
			assert(op_ >= EquationTermType::OP_H0);
			op_ -= EquationTermType::OP_H0;
			setBracketIdentity();
			return;
		}

		bracket_="<" + bra.stateString() + "|H0|" + ket.stateString() + ">";
	}

	void setBracketV()
	{
		bracket_="<" + bra.stateString() + "|V|" + ket.stateString() + ">";
	}

	EquationTermType ket;
	EquationTermType bra;
	RealType factor_;
	EnergyType a_;
	SizeType op_;
	PsimagLite::String energy_;
	PsimagLite::String bracket_;

}; // class Bracket

template<typename RealType>
std::ostream& operator<<(std::ostream& os,const Bracket<RealType>& eqterm)
{
	if (eqterm.factor()==0)
		return os;
	else if (eqterm.factor()==1)
		os<<"";
	else if (eqterm.factor()== -1)
		os<<"- ";
	else
		os<<eqterm.factor();

	os<<eqterm.a();

//	switch (eqterm.order()) {
//	case 0:
//		break;
//	case 1:
//		os<<"l*";
//		break;
//	default:
//		os<<"l^"<<eqterm.order()<<"*";
//	}

	os<<eqterm.energy();

	os<<eqterm.op();

	return os;
}

template<typename RealType>
bool operator==(const Bracket<RealType>& a,const Bracket<RealType>& b)
{
	return (a.sameType(b));
}

} // namespace QuantumPerturbation

#endif // BRACKET_H

