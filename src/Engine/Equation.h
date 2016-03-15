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
#ifndef EQUATION_H
#define EQUATION_H

#include "Vector.h"
#include "BracketCollection.h"

namespace QuantumPerturbation {

template<typename RealType>
class Equation {

public:

	typedef BracketCollection<RealType> BracketCollectionType;
	typedef Equation<RealType> ThisType;
	typedef typename BracketCollectionType::EquationTermType EquationTermType;
	typedef typename PsimagLite::Vector<EquationTermType>::Type FullTermType;
	typedef typename BracketCollectionType::BracketType BracketType;
	typedef typename EquationTermType::StateType StateType;
	typedef typename EquationTermType::EnergyType EnergyType;

	Equation()
	: orderMax_(0)
	{}

	Equation(const FullTermType& left,const FullTermType& right)
	: orderMax_(0)
	{
		add(left,1.0);
		add(right,-1.0);
		orderMax_=maxOrderOf();
	}

	void getEquationAtOrder(ThisType& equation,SizeType order) const
	{
		for (SizeType i=0;i<data_.size();i++) {
			if (data_[i].order()!=order) continue;
			equation.push_back(data_[i]);
		}
	}

	BracketCollectionType findEnergy(int en) const
	{
		FullTermType energy;
		arrangeForEnergy(energy,en);
		StateType state('n',0);
		EnergyType energy1('n',-1);
		EquationTermType n0(EquationTermType::OP_IDENTITY,state,energy1,0);

		BracketCollectionType fullBracket;
		for (SizeType i=0;i<energy.size();i++) {
			BracketType bracket(energy[i],n0);
			fullBracket.push_back(bracket);
		}
		fullBracket.simplify();
		fullBracket.replace();
		return fullBracket;
	}

	BracketCollectionType findState(SizeType ind) const
	{
		BracketCollectionType fullBracket;
		StateType state('m',0);
		EnergyType energy1('m',-1);
		EquationTermType n0(EquationTermType::OP_IDENTITY,state,energy1,0);
		n0 *= (-1.0);
		for (SizeType i=0;i<data_.size();i++) {
			if (data_[i].op() == EquationTermType::OP_H0) continue;
			PsimagLite::String str("n(");
			str += ttos(data_[i].stateIndex());
			str += ")";
			if (data_[i].op() == n0.op() &&
			    data_[i].stateString() == str &&
			    data_[i].energyString() == "En(0)") continue;
			BracketType bracket(data_[i],n0);
			fullBracket.push_back(bracket);
		}

		fullBracket.simplify();
		fullBracket.replace();
		return fullBracket;
	}

	void push_back(const EquationTermType& t)
	{
		data_.push_back(t);
	}

	const EquationTermType& operator[](SizeType i) const
	{
		assert(i<data_.size());
		return data_[i];
	}

	void print(std::ostream& os) const
	{
		for (SizeType i=0;i<data_.size();++i) {
			os<<data_[i];
			if (i+1 < data_.size() && !data_[i+1].hasNegativeSign())
				os<<" + ";
		}
		os<<"= 0\n";
	}

private:

	void add(const FullTermType& x,const RealType& factor)
	{
		for (SizeType i=0;i<x.size();i++) {
			SizeType len = data_.size();
			data_.push_back(x[i]);
			data_[len] *= factor;
		}
	}

	void arrangeForEnergy(FullTermType& newTerm,int en) const
	{
		SizeType count = 0;
		SizeType start = 0;
		for (SizeType i=0;i<data_.size();i++) {
			if (data_[i].energyIndex() == en) {
				if (count>0) {
					PsimagLite::String str(__FILE__);
					str += " " + ttos(__LINE__) + "\n";
					str += "More than one term with energy " + ttos(en) + " found\n";
					throw std::runtime_error(str.c_str());
				}
				count++;
			} else {
				start = newTerm.size();
				newTerm.push_back(data_[i]);
			}
		}
		if (count==0) return;
		for (SizeType i=start;i<newTerm.size();i++) {
			newTerm[i] *= (-1.0);
		}
	}

	SizeType maxOrderOf() const
	{
		SizeType order = 0;
		for (SizeType i=0;i<data_.size();i++) {
			if (data_[i].order()>order)
				order=data_[i].order();
		}
		return order;
	}

	FullTermType data_;
	SizeType orderMax_;

}; // Equation

template<typename RealType>
std::ostream& operator<<(std::ostream& os,Equation<RealType>& eq)
{
	eq.print(os);
	return os;
}

} // QuantumPerturbation

#endif // EQUATION_H

