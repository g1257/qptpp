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
#ifndef HUBBARD_ONE_BAND_EXPAND_HOPPING
#define HUBBARD_ONE_BAND_EXPAND_HOPPING

#include <iostream>
#include "Vector.h"
#include "../../FreeFermions/src/GeometryParameters.h"
#include "../../FreeFermions/src/GeometryLibrary.h"
#include "TypeToString.h"
#include "../../FreeFermions/src/HilbertState.h"
#include "Sort.h"
#include "BitManip.h"
#include "HilbertOneBandExpandHopping.h"

namespace QuantumPerturbation {

template<typename RealType_>
class HubbardOneBandExpandHopping {

	typedef PsimagLite::Matrix<RealType_> MatrixType;
	typedef PsimagLite::InputNg<FreeFermions::InputCheck> InputNgType;
	typedef FreeFermions::GeometryParameters<RealType_,InputNgType::Readable> GeometryParamsType;
	typedef HilbertOneBandExpandHopping<GeometryParamsType> HilbertSpaceType;
	typedef FreeFermions::GeometryLibrary<MatrixType,GeometryParamsType> GeometryLibraryType;
	typedef typename HilbertSpaceType::VectorSizeType VectorSizeType;
	typedef typename PsimagLite::Vector<RealType_>::Type VectorRealType;
	typedef typename HilbertSpaceType::PairSizeType PairSizeType;
	typedef typename PsimagLite::Vector<PairSizeType>::Type VectorPairSizeType;
	typedef typename HilbertSpaceType::Sset SsetType;

public:

	typedef RealType_ RealType;

	HubbardOneBandExpandHopping(InputNgType::Readable& io)
	    : geometryParams_(io),
	      geometry_(geometryParams_),
	      ne_(2,GeometryParamsType::readElectrons(io,geometryParams_.sites)),
	      hilbertSpace_(geometryParams_,ne_),
	      prevSectorEmpty_(false)

	{
		io.read(u_,"hubbardU");

		if (!isConstantVector(u_)) {
			throw PsimagLite::RuntimeError("All hubbardU values must be the same\n");
		}

		/* compute size of basis */
		SizeType n = geometryParams_.sites;
		SizeType npart = ne_[0];
		SizeType fullHilbertSize = hilbertSpace_.calcCombinations(n,npart);

		VectorSizeType data(fullHilbertSize);

		/* define basis states */
		SizeType ket = (1ul<<npart)-1;
		for (SizeType i=0;i<fullHilbertSize;i++) {
			data[i] = ket;
			n=0;
			SizeType m = 0;
			for (;(ket&3)!=1;n++,ket>>=1) {
				m += ket&1;
			}
			ket = ((ket+1)<<n) ^ ((1<<m)-1);
		}

		/* both spin sectors */
		VectorRealType energies2(data.size()*data.size());
		computeEnergies(energies2,data);

		/* order by energy */
		PsimagLite::Sort<VectorRealType> sort2;
		VectorSizeType iperm(energies2.size());
		sort2.sort(energies2,iperm);

		iperm.clear();

		cacheOffsets(energies2);

		setSsetTypes();

		std::cerr<<geometry_;
		std::cerr<<"Energy=0\n";
	}

	SizeType total() const
	{
		assert(offsets_.size() >= 1);
		return offsets_.size() - 1;
	}

	RealType v(SizeType i,SizeType j) const
	{
		SizeType diff = (i >= j) ? i - j : j - i;

		if (diff > 1) return 0.0;

		if (i >= eigenvector_.size()) {
			computeEigenvector(i);
			if (i >= eigenvector_.size()) return 0.0;
		}

		if (j >= eigenvector_.size()) {
			computeEigenvector(j);
			if (j >= eigenvector_.size()) return 0.0;
		}

		const VectorRealType& vx = eigenvector_[i];
		const VectorRealType& vy = eigenvector_[j];
		SizeType totali = vx.size();
		SizeType totalj = vy.size();

		RealType sum = 0.0;
		for (SizeType x = 0; x < totali; ++x) {
			PairSizeType bra = hilbertSpace_.getKet(x,sSetOf(i));
			for (SizeType y = 0; y < totalj; ++y) {
				PairSizeType ket = hilbertSpace_.getKet(y,sSetOf(j));
				sum += PsimagLite::conj(vx[x]) * vy[y] * value(ket,bra);
			}
		}

		return sum;
	}

	RealType e(SizeType i,SizeType j) const
	{
		assert(i != j);
		return (i - j) * u_[0];
	}

private:

	void computeMatrix(MatrixType& m, SizeType ind) const
	{
		if (ind == 0)
			return computeMatrix00(m);

		if (prevSectorEmpty_)
			return computeRectangularMatrix(m,ind);

		SizeType qsize = 2 * ind;
		assert(qsize > 0);
		VectorRealType q(qsize,0);
		SizeType totalqs = q.size() - 1;
		VectorRealType max(q.size());

		setMax(max);
		VectorRealType r(max[totalqs]);

		bool flag = true;
		while (flag) {
			assert(q[totalqs] < r.size());
			r[q[totalqs]] += getR(q,ind);
			flag = advanceVector(q,max);
		}

		m.resize(r.size(),r.size());

		int sign = (ind & 1) ? -1 : 1;
		for (SizeType i = 0; i < r.size(); ++i) {
			for (SizeType j = 0; j < r.size(); ++j) {
				m(i,j) = r[i] * r[j] * sign;
			}
		}
	}

	RealType getR(const VectorRealType& q, SizeType sector) const
	{
		SizeType totalqs = static_cast<SizeType>(q.size()*0.5);
		RealType prod = 1.0;
		for (SizeType nq = 0; nq < totalqs; ++nq) {
			assert(2*nq + 1 < q.size());
			PairSizeType bra = hilbertSpace_.getKet(q[2*nq],sSetOf(nq));
			PairSizeType ket = hilbertSpace_.getKet(q[2*nq+1],sSetOf(nq+1));
			assert(nq < eigenvector_.size());
			assert(q[nq] < eigenvector_[nq].size());
			prod *= value(bra,ket) * eigenvector_[nq][q[nq]];
		}

		return prod;
	}

	void computeMatrix00(MatrixType& m) const
	{
		VectorPairSizeType basis;
		setBasis(basis,0);

		SizeType n = basis.size();
		m.resize(n,n);

		for (SizeType i = 0; i < n; ++i) {
			PairSizeType ket = basis[i];
			for (SizeType j = 0; j < n; ++j) {
				PairSizeType bra = basis[j];

				m(i,j) = value(ket,bra);

			}
		}

		assert(isHermitian(m));
	}

	bool advanceVector(VectorRealType& q, const VectorRealType& max) const
	{
		for (SizeType i = 0; i < q.size(); ++i) {
			q[i]++;
			if (q[i] < max[i]) return true;
			q[i] = 0;
		}

		return false;
	}

	void setMax(VectorRealType& max) const
	{
		SizeType j = 0;
		for (SizeType i = 0; i < max.size(); ++i) {
			assert(j+1 < offsets_.size());
			max[i] = offsets_[j+1] - offsets_[j];
			if (i & 1) continue;
			j++;
		}
	}

	bool isAlmostZero(const MatrixType& m) const
	{
		RealType sum = 0.0;
		for (SizeType i = 0; i < m.n_row(); ++i)
			for (SizeType j = 0; j < m.n_row(); ++j)
				sum += PsimagLite::conj(m(i,j)) * m(i,j);

		return (sum < 1e-6);
	}

	void computeEigenvector(SizeType thisSector) const
	{
		MatrixType m;
		computeMatrix(m,thisSector);

		if (isAlmostZero(m)) {
			if (prevSectorEmpty_) {
				throw PsimagLite::RuntimeError("computeEigenvector\n");
			}

			prevSectorEmpty_ = true;
			return;
		}

		if (m.n_row() == m.n_col()) {
			SizeType n = m.n_row();
			VectorRealType eigs(n);
			diag(m,eigs,'V');

			VectorRealType v(m.n_row());
			for (SizeType i = 0; i < m.n_row(); ++i)
				v[i] = m(i,0);
			assert(eigenvector_.size() == thisSector);
			eigenvector_.push_back(v);
		} else {
			MatrixType vt;
			VectorRealType s;
			svd('A',m,s,vt);

			assert(thisSector > 0);
			SizeType prevSector = thisSector -1;
			SizeType thisSectorSize = offsets_[thisSector + 1] - offsets_[thisSector];
			SizeType prevSectorSize = offsets_[thisSector] - offsets_[prevSector];
			assert(m.n_row() == prevSectorSize);
			VectorRealType v(prevSectorSize);
			for (SizeType i = 0; i < m.n_row(); ++i)
				v[i] = m(i,0);
			assert(eigenvector_.size() == prevSector);
			eigenvector_.push_back(v);

			v.resize(thisSectorSize);
			assert(vt.n_col() == thisSectorSize);
			for (SizeType i = 0; i < vt.n_col(); ++i)
				v[i] = vt(0,i);
			assert(eigenvector_.size() == thisSector);
			eigenvector_.push_back(v);
		}

		prevSectorEmpty_ = false;
	}

	void computeRectangularMatrix(MatrixType& m, SizeType thisSector) const
	{
		if (thisSector != 1)
			throw PsimagLite::RuntimeError("computeRectangularMatrix");

		SizeType prevSector = thisSector - 1;
		SizeType thisSize = offsets_[thisSector + 1] - offsets_[thisSector];
		SizeType prevSize = offsets_[thisSector] - offsets_[prevSector];

		m.resize(prevSize,thisSize);

		int sign = (thisSector & 1) ? -1 : 1;
		for (SizeType i = 0; i < prevSize; ++i) {
			PairSizeType bra = hilbertSpace_.getKet(i,sSetOf(prevSector));
			for (SizeType j = 0; j < thisSize; ++j) {
				PairSizeType ket = hilbertSpace_.getKet(j,sSetOf(thisSector));
				m(i,j) = sign * value(bra,ket);
			}
		}
	}

	void computeEnergies(VectorRealType& energies2,const VectorSizeType& data) const
	{
		SizeType n = data.size();
		for (SizeType i = 0; i < n; ++i) {
			for (SizeType j = 0; j < n; ++j) {
				energies2[i+j*n] = calcEnergy(data[i],data[j]);
			}
		}
	}

	RealType calcEnergy(SizeType up, SizeType down) const
	{
		SizeType mask = (up & down);
		SizeType site = 0;
		RealType sum = 0;

		while (mask > 0) {
			if (mask & 1) sum += u_[site];
			site++;
			mask >>= 1;
		}

		return sum;
	}

	void cacheOffsets(const VectorRealType& energies2)
	{
		RealType previousEnergy = energies2[0] - 100.;
		for (SizeType i = 0; i < energies2.size(); ++i) {
			if (energies2[i] == previousEnergy) continue;

			offsets_.push_back(i);
			previousEnergy = energies2[i];
		}
		offsets_.push_back(energies2.size());
	}

	RealType value(PairSizeType ind, PairSizeType jnd) const
	{
		if (ind.first == jnd.first)
			return value(ind.second,jnd.second);

		if (ind.second == jnd.second)
			return value(ind.first,jnd.first);

		return 0.0;
	}

	RealType value(SizeType ind, SizeType jnd) const
	{
		SizeType mask (ind ^ jnd);
		SizeType c = PsimagLite::BitManip::count(mask);

		if (c != 2) return 0.0;

		SizeType i = 0;
		SizeType j = 0;
		SizeType counter = 0;
		SizeType level = 0;

		while (mask > 0) {

			if (mask & 1) {
				if (counter == 0) {
					i = level;
					counter++;
				} else if (counter == 1) {
					j = level;
					counter++;
					break;
				}
			}

			mask >>= 1;
			level++;
		}

		return geometry_.matrix()(i,j);
	}

	void setSsetTypes()
	{
		SizeType doublyOccupiedPlus1 = *std::min_element(ne_.begin(),ne_.end());
		doublyOccupiedPlus1++;

		for (SizeType s = 0; s < doublyOccupiedPlus1; ++s)
			setSsetType(s);
	}

	void setSsetType(SizeType s)
	{
		VectorSizeType length(3);
		SizeType linSize = geometryParams_.sites;
		length[0] = hilbertSpace_.calcCombinations(linSize,s);
		SizeType prod = length[0];
		length[1] = hilbertSpace_.calcCombinations(linSize - s,ne_[0] -s);
		prod *= length[1];
		length[2] = hilbertSpace_.calcCombinations(linSize - ne_[0],ne_[1] -s);
		prod *= length[2];

		SizeType offset = (s == 0) ? 0 : sSet_[s-1].offset + prod;

		SsetType sSet(s,offset,length);
		sSet_.push_back(sSet);
	}

	void setBasis(VectorPairSizeType& basis,SizeType indexOfSet) const
	{
		SizeType total = offsets_[indexOfSet+1] - offsets_[indexOfSet];

		basis.resize(total);

		for (SizeType i = 0; i < total; ++i)
			basis[i] = hilbertSpace_.getKet(i,sSetOf(indexOfSet));
	}

	SizeType whatPartition(SizeType ind) const
	{
		assert(offsets_.size()>0);
		SizeType last = offsets_.size() - 1;
		for (int i = last; i>=0; --i)
			if (ind >= offsets_[i]) return i;

		throw PsimagLite::RuntimeError("whatPartition");
	}

	bool isConstantVector(const VectorRealType& v) const
	{
		if (v.size() == 0) return true;

		for (SizeType i = 0; i < v.size(); ++i)
			if (v[i] != v[0]) return false;

		return true;
	}

	const SsetType& sSetOf(SizeType ind) const
	{
		assert(ind < sSet_.size());
		return sSet_[ind];
	}

	GeometryParamsType geometryParams_;
	GeometryLibraryType geometry_;
	PsimagLite::Vector<size_t>::Type ne_; // n. of up (= n. of  down electrons)
	HilbertSpaceType hilbertSpace_;
	mutable bool prevSectorEmpty_;
	VectorRealType u_;
	VectorSizeType offsets_;
	mutable typename PsimagLite::Vector<VectorRealType>::Type eigenvector_;
	typename PsimagLite::Vector<SsetType>::Type sSet_;
};

} // namespace QuantumPerturbation

#endif // HUBBARD_ONE_BAND_EXPAND_HOPPING
