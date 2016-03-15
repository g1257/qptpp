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
#ifndef VECTOR_OF_ENERGIES_H
#define VECTOR_OF_ENERGIES_H

#include "Vector.h"
#include <utility>
#include <cstdlib>

namespace QuantumPerturbation {

template<typename RealType>
class VectorOfEnergies {

	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef std::pair<SizeType,SizeType> PairSizeType;

public:

	typedef RealType value_type;

	VectorOfEnergies(const VectorRealType& energies,
	               const PsimagLite::Vector<PairSizeType>::Type& iperm)
	    : energies_(energies),iperm_(iperm),n_(energies.size()),dummy_(0)
	{}

	const RealType& operator[](SizeType ind) const
	{
		const PairSizeType& p = getPermOf(ind);

		dummy_ = energies_[p.first] + energies_[p.second];
		return dummy_;
	}

	RealType& operator[](SizeType ind)
	{
		return dummy_;
	}

	SizeType size() const
	{
		return n_*n_;
	}

private:

	PairSizeType getPermOf(SizeType ind) const
	{
		if (iperm_.size()>0) return iperm_[ind];

		div_t q = div(ind,n_);
		return PairSizeType(q.quot,q.rem);
	}

	const VectorRealType& energies_;
	const PsimagLite::Vector<PairSizeType>::Type& iperm_;
	SizeType n_;
	mutable RealType dummy_;
};

} // namespace QuantumPerturbation

#endif // VECTOR_OF_ENERGIES_H

