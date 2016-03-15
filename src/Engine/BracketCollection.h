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
#ifndef BRACKETCOLLECTION_H
#define BRACKETCOLLECTION_H

#include "Bracket.h"
#include "Vector.h"

namespace QuantumPerturbation {

template<typename RealType_>
class BracketCollection {

	typedef typename PsimagLite::Vector<Bracket<RealType_> >::Type VectorType;
	typedef typename VectorType::iterator IteratorType;

public:

	typedef RealType_ RealType;
	typedef Bracket<RealType> BracketType;
	typedef typename BracketType::EquationTermType EquationTermType;

	void push_back(const BracketType& item)
	{
		data_.push_back(item);
	}

	void simplify()
	{
		VectorType bag;
		for (SizeType i=0;i<data_.size();i++) {
			if (data_[i].isNull()) continue;
			IteratorType it = find(bag.begin(),bag.end(),data_[i]);
			if (it!=bag.end()) {
				*it += data_[i];
			} else {
				bag.push_back(data_[i]);
			}
		}
		data_=bag;
	}

	void replace()
	{
		for (SizeType i=0;i<data_.size();i++)
			data_[i].replace();
	}

	void print(std::ostream& os) const
	{
		for (SizeType i=0;i<size();i++) {
			os<<data_[i];
			if (i+1<size() && data_[i+1].factor()>0)
				os<<" + ";
		}
	}

	const BracketType& operator()(SizeType i) const
	{
		assert(i<data_.size());
		return data_[i];
	}

	const SizeType size() const { return data_.size(); }

private:

	VectorType data_;

}; // BracketCollection

template<typename RealType>
std::ostream& operator<<(std::ostream& os,const BracketCollection<RealType>& eqterm)
{
	eqterm.print(os);
	return os;
}

} // namespace QuantumPerturbation

#endif // BRACKETCOLLECTION_H

