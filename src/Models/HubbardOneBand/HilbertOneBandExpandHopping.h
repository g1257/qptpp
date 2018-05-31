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
#ifndef HILBERT_ONE_BAND_EXPAND_HOPPING
#define HILBERT_ONE_BAND_EXPAND_HOPPING

#include <iostream>
#include "Vector.h"
#include "TypeToString.h"
#include "BitManip.h"

namespace QuantumPerturbation {

template<typename GeometryParamsType>
class HilbertOneBandExpandHopping {

	typedef typename GeometryParamsType::RealType RealType;

public:

	typedef std::pair<SizeType,SizeType> PairSizeType;
	typedef PsimagLite::Vector<PairSizeType>::Type VectorPairSizeType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef PsimagLite::Vector<VectorSizeType>::Type VectorVectorSizeType;

	struct Sset {

		Sset(SizeType s1, SizeType offset1, const VectorSizeType& length1)
		    : s(s1),offset(offset1),length(length1)
		{}

		SizeType s;
		SizeType offset;
		VectorSizeType length;
	}; // struct Sset

	HilbertOneBandExpandHopping(const GeometryParamsType& geometryParams,
	                            const PsimagLite::Vector<size_t>::Type& ne)
	    : geometryParams_(geometryParams),
	      ne_(ne),
	      kets_(geometryParams_.sites + 1, geometryParams_.sites + 1)
	{
		SizeType n = geometryParams_.sites;

		for (SizeType i = 1; i <= n; ++i) {
			for (SizeType j = 0; j <= i; ++j) {
				VectorSizeType* v = new VectorSizeType;
				findKet(*v,i,j);
				kets_(i,j) = v;
			}
		}
	}

	~HilbertOneBandExpandHopping()
	{
		SizeType n = geometryParams_.sites;
		for (SizeType i = 1; i <= n; ++i) {
			for (SizeType j = 0; j <= i; ++j) {
				delete kets_(i,j);
				kets_(i,j) = 0;
			}
		}
	}

	PairSizeType getKet(SizeType ind, const Sset& sSet) const
	{
		SizeType s = sSet.s;

		if (s >= cachedKets_.size()) {
			cacheKets(sSet);
		}

		assert(s < cachedKets_.size());
		assert(ind < cachedKets_[s].size());
		return cachedKets_[s][ind];
	}

	SizeType calcCombinations(SizeType n, SizeType npart) const
	{
		SizeType fullHilbertSize = 1;

		for (SizeType m=1; m <= npart; n--,m++)
			fullHilbertSize=fullHilbertSize*n/m;

		return fullHilbertSize;
	}

private:

	void cacheKets(const Sset& sSet) const
	{
		SizeType s = sSet.s;
		SizeType total = sSet.length[0] * sSet.length[1] * sSet.length[2];
		VectorPairSizeType v(total);

		for (SizeType ind = 0; ind < total; ++ind) {
			const VectorSizeType& coordinates = findCoordinates(ind,sSet);
			PairSizeType ket;
			for (SizeType i = 0; i < coordinates.size(); i++) {
				PairSizeType prevKet = ket;
				findKet(ket,i,prevKet,sSet,coordinates[i]);
			}

			assert(PsimagLite::BitManip::count(ket.first & ket.second) == int(s));
			assert(PsimagLite::BitManip::count(ket.first) == int(ne_[0]));
			assert(PsimagLite::BitManip::count(ket.second) == int(ne_[1]));

			v[ind] = ket;
		}

		if (cachedKets_.size() != s)
			throw PsimagLite::RuntimeError("cacheKets");

		cachedKets_.push_back(v);

	}

	void cacheCoordinates(const Sset& sSet) const
	{
		SizeType s = sSet.s;
		SizeType total = sSet.length[0] * sSet.length[1] * sSet.length[2];
		std::cerr<<"cacheCoordinates s="<<s<<"\n";

		VectorVectorSizeType v;
		for (SizeType ind = 0; ind < total; ++ind) {
			VectorSizeType coordinates(3);
			findCoordinates_(coordinates,ind,sSet);
			v.push_back(coordinates);
		}

		if (cachedCoordinates_.size() != s)
			throw PsimagLite::RuntimeError("cacheCoordinates");

		cachedCoordinates_.push_back(v);
	}

	const VectorSizeType& findCoordinates(SizeType ind, const Sset& sSet) const
	{
		SizeType s = sSet.s;
		if (s >= cachedCoordinates_.size()) {
			cacheCoordinates(sSet);
		}

		assert(ind < cachedCoordinates_[s].size());
		return cachedCoordinates_[s][ind];
	}

	void findCoordinates_(VectorSizeType& coordinates,
	                      SizeType ind,
	                      const Sset& sSet) const
	{
		assert(coordinates.size() == 3);

		assert(ind < sSet.length[1] * sSet.length[0] * sSet.length[2]);
		SizeType n = sSet.length[1] * sSet.length[0];
		div_t q = div(ind,n);
		ind = q.rem;
		coordinates[2] = q.quot;

		n = sSet.length[0];
		q = div(ind,n);

		coordinates[0] = q.rem;
		coordinates[1] = q.quot;

		for (SizeType i = 0; i < coordinates.size(); ++i)
			assert(coordinates[i] < sSet.length[i]);
	}

	void findKet(PairSizeType& ket,
	             SizeType coordinateIndex,
	             const PairSizeType& prevKet,
	             const Sset& sSet,
	             SizeType coordinate) const
	{
		SizeType total = 0;
		SizeType take = 0;
		SizeType linSize = geometryParams_.sites;

		if (coordinateIndex == 0) {
			total = linSize;
			take = sSet.s;
		} else if (coordinateIndex == 1) {
			assert(sSet.s < linSize);
			total = linSize - sSet.s;
			take = ne_[0] - sSet.s;
		} else if (coordinateIndex == 2) {
			total = linSize - ne_[0];
			take = ne_[1] - sSet.s;
		} else {
			assert(false);
		}

		findKet(ket.first,total,take,coordinate);
		collateKet(ket,coordinateIndex,prevKet);
	}

	void findKet(SizeType& ket,
	             SizeType total,
	             SizeType npart,
	             SizeType coordinate) const
	{
		assert(coordinate < kets_(total,npart)->size());
		ket = kets_(total,npart)->operator[](coordinate);
	}

	void findKet(VectorSizeType& v,
	             SizeType total,
	             SizeType npart) const
	{
		SizeType ket = (1ul<<npart)-1;

		if (npart == 0) {
			v.push_back(ket);
			return;
		}

		SizeType fullHilbertSize = calcCombinations(total,npart);

		for (SizeType i=0;i<fullHilbertSize;i++) {
			v.push_back(ket);
			SizeType n=0;
			SizeType m = 0;
			for (;(ket&3)!=1;n++,ket>>=1) {
				m += ket&1;
			}
			ket = ((ket+1)<<n) ^ ((1<<m)-1);
		}
	}

	void collateKet(PairSizeType& ket,
	                SizeType coordinateIndex,
	                const PairSizeType& prevKet) const
	{
		if (coordinateIndex == 0) {
			ket.second = ket.first;
			return;

		} else if (coordinateIndex == 1) {
			assert(prevKet.first == prevKet.second);

			SizeType result = prevKet.first;
			addToKet(result,prevKet.first,ket.first);
			ket.first = result;
			ket.second = prevKet.second;

		} else if (coordinateIndex == 2) {
			SizeType result = prevKet.second;
			addToKet(result,(prevKet.first | prevKet.second),ket.first);
			ket.first = prevKet.first;
			ket.second = result;

		} else {
			assert(false);
		}
	}

	void addToKet(SizeType& dest,
	              SizeType prev,
	              SizeType whatToAdd) const
	{
		SizeType offset = 0;
		SizeType level = 0;

		while (whatToAdd > 0) {
			if (prev > 0) {
				while (prev & 1) {
					offset++;
					prev >>= 1;
				}

				prev >>= 1;
			}

			if (whatToAdd & 1) {
				SizeType mask = (1 << (offset + level));
				assert((dest & mask) == 0);
				dest |= mask;
			}

			level++;
			whatToAdd >>= 1;
		}
	}

	const GeometryParamsType& geometryParams_;
	const PsimagLite::Vector<size_t>::Type& ne_; // n. of up (= n. of  down electrons)
	PsimagLite::Matrix<VectorSizeType*> kets_;
	mutable PsimagLite::Vector<VectorVectorSizeType>::Type cachedCoordinates_;
	mutable PsimagLite::Vector<VectorPairSizeType>::Type cachedKets_;
};

} // namespace QuantumPerturbation

#endif // HILBERT_ONE_BAND_EXPAND_HOPPING

