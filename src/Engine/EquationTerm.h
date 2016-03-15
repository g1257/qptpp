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
#ifndef EQUATIONTERM_H
#define EQUATIONTERM_H

#include "Vector.h"

namespace QuantumPerturbation {

template<typename RealType>
class EquationTerm {

public:

	typedef std::pair<char,SizeType> StateType;
	typedef std::pair<char,int> EnergyType;

	enum {OP_IDENTITY=0, OP_H0=1, OP_V=2};

	static const int ENERGY_NONE = -1;

	EquationTerm(SizeType op1,
	             const std::pair<char,SizeType>& state1,
	             const std::pair<char,int>& energy1,
	             SizeType order1)
	: op_(op1),
	  state_(state1),
	  energy_(energy1),
	  order_(order1),
	  factor_(1.0),
	  hasNegativeSign_(factor_ < 0)
	{}

	EquationTerm& operator*=(const RealType& val)
	{
		factor_ *=val;
		hasNegativeSign_ = (factor_ < 0);
		return *this;
	}

	void state(const StateType& state1)
	{
		state_ = state1;
	}

	void print(std::ostream& os) const
	{
		if (factor_ == 1)
			os<<"";
		else if (factor_ == -1)
			os<<" - ";
		else
			os<< factor_<<"*";

	//	switch (order) {
	//	case 0:
	//		break;
	//	case 1:
	//		os<<"l*";
	//		break;
	//	default:
	//		os<<"l^"<<order<<"*";
	//	}

		if (energy_.second >= 0)
			os<<"E"<<energy_.first<<"("<<energy_.second<<")*";
		switch (op_) {
		case OP_IDENTITY:
			break;
		case OP_H0:
			os<<"H0";
			break;
		case OP_V:
			os<<"V";
			break;
		default:
			assert(false);
		}
		os<<"|"<<state_.first<<"("<<state_.second<<")> ";
	}

	const SizeType& order() const { return order_; }

	const RealType& factor() const { return factor_; }

	bool hasNegativeSign() const { return hasNegativeSign_; }

	const int& energyIndex() const { return energy_.second; }

	PsimagLite::String energyString() const
	{
		if (energy_.second < 0) return "";

		return PsimagLite::String("E") +
		        energy_.first + "(" + ttos(energy_.second) + ")";
	}

	const SizeType& op() const { return op_; }

	PsimagLite::String stateString() const
	{
		PsimagLite::String str1(" ");
		str1[0] = state_.first;
		str1 += "(";
		str1 += ttos(state_.second);
		str1 += ")";
		return str1;
	}

	const SizeType& stateIndex() const
	{
		return state_.second;
	}

	const char& stateChar() const { return state_.first; }

private:

	SizeType op_;
	std::pair<char,SizeType> state_;
	std::pair<char,int> energy_;
	SizeType order_;
	RealType factor_;
	bool hasNegativeSign_;

}; // EquationTerm

template<typename T>
struct IsEquationTerm {
	enum {True = false};
};

template<typename RealType>
struct IsEquationTerm<EquationTerm<RealType> > {
	enum {True = true};
};

template<typename RealType>
std::ostream& operator<<(std::ostream& os,
                         const EquationTerm<RealType>& eqterm)
{
	eqterm.print(os);

	return os;
}

template<typename VectorEquationTermType>
typename PsimagLite::EnableIf<
PsimagLite::IsVectorLike<VectorEquationTermType>::True &
IsEquationTerm<typename VectorEquationTermType::value_type>::True,
void>::Type print(std::ostream& os,const VectorEquationTermType& eqterm)
{
	//os<<"Vector.size="<<eqterm.size()<<"\n";
	for (SizeType i=0;i<eqterm.size();i++) {
		os<<eqterm[i];
		if (i+1<eqterm.size() && eqterm[i+1].factor()>0)
			os<<"+";
	}
}
} //namespace QuantumPerturbation

#endif // EQUATIONTERM_H

