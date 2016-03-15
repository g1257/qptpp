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
#ifndef EXPANSION_H
#define EXPANSION_H

#include "Equation.h"
#include "ExpansionHelper.h"

namespace QuantumPerturbation {

template<typename ModelType>
class Expansion {

public:

	typedef typename ModelType::RealType RealType;
	typedef Equation<RealType> EquationType;
	typedef typename EquationType::BracketCollectionType BracketCollectionType;
	typedef ExpansionHelper<BracketCollectionType,ModelType> ExpansionHelperType;
	typedef typename EquationType::FullTermType FullTermType;
	typedef typename EquationType::EquationTermType EquationTermType;
	typedef typename EquationTermType::StateType StateType;
	typedef typename EquationTermType::EnergyType EnergyType;

	Expansion(ModelType& model,SizeType mode,SizeType order)
	    : helper_(model,mode),equations_(order+1)
	{}

	void expand()
	{
		FullTermType left;
		FullTermType right;

		for (SizeType term=1;term<3;term++) {
			for (SizeType st=0;st<equations_.size();st++) {
				SizeType order = term-1+st;
				if (order>equations_.size()-1) continue;
				StateType state('n',st);
				EnergyType energy('n',-1);
				EquationTermType eqterm(term,state,energy,order);
				left.push_back(eqterm);
			}
		}
		for (SizeType en=0;en<equations_.size();en++) {
			for (SizeType st=0;st<equations_.size();st++) {
				SizeType order = en+st;
				if (order>equations_.size()-1) continue;
				StateType state('n',st);
				EnergyType energy('n',en);
				EquationTermType eqterm(EquationTermType::OP_IDENTITY,state,energy,order);
				right.push_back(eqterm);
			}
		}
		EquationType equation(left,right);

		for (SizeType i=0;i<equations_.size();i++) {
			equation.getEquationAtOrder(equations_[i],i);
			std::cout<<equations_[i];
		}

		helper_.printHeader(std::cout);

		for (SizeType i=1;i<equations_.size();i++) {
			BracketCollectionType sol = equations_[i].findEnergy(i);
			helper_.push("En(" + ttos(i) + ") = ",sol);

			BracketCollectionType sol2 = equations_[i].findState(i);
			helper_.push("E(m,n) am(" + ttos(i) + ") = ",sol2);
		}
	}

private:

	ExpansionHelperType helper_;
	typename PsimagLite::Vector<EquationType>::Type equations_;

}; // class Expansion

} // namespace QuantumPerturbation

#endif // EXPANSION_H

