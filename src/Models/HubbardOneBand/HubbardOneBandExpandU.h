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
#ifndef HUBBARD_ONE_BAND_EXPAND_U
#define HUBBARD_ONE_BAND_EXPAND_U

#include <iostream>
#include "Vector.h"
#include "../../FreeFermions/src/GeometryParameters.h"
#include "../../FreeFermions/src/Engine.h"
#include "../../FreeFermions/src/GeometryLibrary.h"
#include "TypeToString.h"
#include "../../FreeFermions/src/CreationOrDestructionOp.h"
#include "../../FreeFermions/src/HilbertState.h"
#include "VectorOfEnergies.h"
#include "Sort.h"
#include "BitManip.h"
#include "InputNg.h"
#include "../../FreeFermions/src/InputCheck.h"

namespace QuantumPerturbation {

template<typename RealType_>
class HubbardOneBandExpandU {

	typedef PsimagLite::Matrix<RealType_> MatrixType;
	typedef PsimagLite::InputNg<FreeFermions::InputCheck> InputNgType;
	typedef FreeFermions::GeometryParameters<RealType_,InputNgType::Readable> GeometryParamsType;
	typedef FreeFermions::GeometryLibrary<MatrixType,GeometryParamsType> GeometryLibraryType;
	typedef FreeFermions::Engine<RealType_,RealType_> EngineType;
	typedef FreeFermions::CreationOrDestructionOp<EngineType> OperatorType;
	typedef FreeFermions::HilbertState<OperatorType> HilbertStateType;
	typedef std::pair<SizeType,SizeType> PairSizeType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<RealType_>::Type VectorRealType;

public:

	typedef RealType_ RealType;

	HubbardOneBandExpandU(InputNgType::Readable& io)
	    : geometryParams_(io),
	      geometry_(geometryParams_),
	      engine_(geometry_.matrix(),1,true),
	      ne_(1,GeometryParamsType::readElectrons(io,geometryParams_.sites)),
	      gs_(engine_,ne_)
	{
		/* compute size of basis */
		SizeType fullHilbertSize = 1;
		int n = engine_.size();
		SizeType npart = ne_[0];
		SizeType m=1;
		for (;m<=npart;n--,m++)
			fullHilbertSize=fullHilbertSize*n/m;

		VectorSizeType data(fullHilbertSize);
		energies_.resize(fullHilbertSize);

		/* define basis states */
		SizeType ket = (1ul<<npart)-1;
		for (SizeType i=0;i<data.size();i++) {
			data[i] = ket;
			energies_[i] = energy(ket);
			n=m=0;
			for (;(ket&3)!=1;n++,ket>>=1) {
				m += ket&1;
			}
			ket = ((ket+1)<<n) ^ ((1<<m)-1);
		}

		/* sort one sector */
		PsimagLite::Sort<VectorRealType> sort;
		VectorSizeType iperm(data.size());
		sort.sort(energies_,iperm);

		cacheValues(data,iperm);
		data.clear();

		/* both spin sectors */
		VectorOfEnergies<RealType> energies2(energies_,iperm_);

		/* order by energy */
		PsimagLite::Sort<VectorOfEnergies<RealType> > sort2;
		iperm.resize(energies2.size());
		sort2.sort(energies2,iperm);

		cachePerm(iperm);
		iperm.clear();

		std::cerr<<geometry_;

		RealType sum = 0;
		for (SizeType i=0;i<ne_[0];i++)
			sum += engine_.eigenvalue(i);
		std::cerr<<"Energy="<<sum<<"\n";
//		std::cerr<<"E2 (gues) = "<<guessE2()<<"\n";
	}

	SizeType total() const { return iperm_.size(); }

	RealType v(SizeType i,SizeType j) const
	{
		const PairSizeType& ip = iperm_[i];
		const PairSizeType& jp = iperm_[j];

		RealType sum = 0;
		for (SizeType site = 0; site < engine_.size(); ++site)
			sum += value(data_[ip.first],data_[jp.first],site) *
			       value(data_[ip.second],data_[jp.second],site);

		return sum;
	}

	RealType e(SizeType i,SizeType j) const
	{
		assert(i!=j);
		const PairSizeType& ip = iperm_[i];
		const PairSizeType& jp = iperm_[j];
		RealType result =  energies_[ip.first] - energies_[jp.first] +
		        energies_[ip.second] - energies_[jp.second];
		if (result == 0)
			err("Ground state is degenerate? (not supported yet)\n");

		return result;
	}

private:

	void cacheValues(const VectorSizeType& data,const VectorSizeType& iperm)
	{
		SizeType n = data.size();
		data_.resize(n);
		for (SizeType i = 0; i < n; ++i)
			data_[i] = data[iperm[i]];
	}

	void cachePerm(const VectorSizeType& iperm)
	{
		SizeType n = energies_.size();
		SizeType bigN  = iperm.size();
		assert(n*n == bigN);
		iperm_.resize(bigN);
		for (SizeType i = 0; i < iperm.size(); ++i) {
			SizeType ind = iperm[i];
			div_t q = div(ind,n);
			iperm_[i] = PairSizeType(q.quot,q.rem);
		}
	}

	RealType energy(SizeType ind) const
	{
		SizeType level = 0;
		RealType sum = 0.0;
		while (ind > 0) {
			if (ind & 1) sum += engine_.eigenvalue(level);
			level++;
			ind >>= 1;
		}

		return sum;
	}

	RealType value(SizeType ind, SizeType jnd, SizeType site) const
	{
		PairSizeType levels;
		SizeType c = findDifferences(levels,ind,jnd);

		if (c != 0 && c != 2) return 0.0;

		if (c == 0) {
			assert(ind == jnd);
			return valueEqual(ind,site);
		}

		assert(c == 2);

		SizeType ind2 = ind;
		ind2 ^= (1 << levels.first);
		ind2 ^= (1 << levels.second);
		assert(ind2 == jnd);
		return valueNonEqual(ind,levels,site);
	}

	RealType valueEqual(SizeType ind, SizeType site) const
	{
		RealType sum = 0.0;
		SizeType level1 = 0;

		SizeType ind2 = ind;
		while (ind2 > 0) {
			if (ind2 & 1)
				sum += PsimagLite::conj(eigenvector(level1,site))*eigenvector(level1,site);

			level1++;
			ind2 >>= 1;
		}

		return sum;
	}

	RealType valueNonEqual(SizeType ind,
	                       const PairSizeType& levels,
	                       SizeType site) const
	{
		RealType sum = PsimagLite::conj(eigenvector(levels.first,site))*
		        eigenvector(levels.second,site);

		return sum;
	}

	SizeType findDifferences(PairSizeType& levels,SizeType ind,SizeType jnd) const
	{
		SizeType x = (ind ^ jnd);
		SizeType c = PsimagLite::BitManip::count(x);

		if (c != 0 && c != 2) return c;

		if (c == 0) return c;

		levels.first = findLevel(x & ind);
		levels.second = findLevel(x & jnd);

		assert(levels.first != levels.second);

		return c;
	}

	SizeType findLevel(SizeType ind) const
	{
		assert(PsimagLite::BitManip::count(ind) == 1);

		SizeType y = ind;
		SizeType level = 0;
		while (y > 0) {
			y >>= 1;
			level++;
		}

		assert(level>0);
		return level-1;
	}

	RealType guessE2()
	{
		RealType sum = 0;
		SizeType occupied = ne_[0];
		for (SizeType site = 0; site < engine_.size(); ++site) {
			for (SizeType s = 0; s < occupied; s++) {
				for (SizeType sp = occupied; sp< engine_.size(); ++sp) {
					RealType tmp1 = std::conj(eigenvector(s,site))*
					        eigenvector(sp,site);
					RealType num1 = engine_.eigenvalue(s)-engine_.eigenvalue(sp);
					for (SizeType l = 0; l < occupied; l++) {
						for (SizeType lp = occupied; lp< engine_.size(); ++lp) {

							RealType tmp2 = std::conj(eigenvector(l,site))*
							        eigenvector(lp,site);
							RealType num2 = engine_.eigenvalue(l)-engine_.eigenvalue(lp);
							RealType tmp = tmp1 * tmp2;
							sum += (tmp * tmp) / (num1 + num2);
						}
					}
				}
			}
		}

		return sum;
	}

	const RealType& eigenvector(SizeType i, SizeType j) const
	{
		return engine_.eigenvector(i,j);
	}

	GeometryParamsType geometryParams_;
	GeometryLibraryType geometry_;
	EngineType engine_;
	PsimagLite::Vector<SizeType>::Type ne_; // n. of up (= n. of  down electrons)
	HilbertStateType gs_;
	VectorRealType energies_;
//	typename PsimagLite::Vector<PsimagLite::Matrix<RealType> >::Type values_;
	VectorSizeType data_;
	PsimagLite::Vector<PairSizeType>::Type iperm_;
};

} // namespace QuantumPerturbation

#endif // HUBBARD_ONE_BAND_EXPAND_U

