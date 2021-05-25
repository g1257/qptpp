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
#ifndef EXPANSION_HELPER_H
#define EXPANSION_HELPER_H

#include <iostream>
#include "Vector.h"

namespace QuantumPerturbation {

template<typename BracketCollectionType,typename ModelType>
class ExpansionHelper {

	typedef typename BracketCollectionType::RealType RealType_;
	typedef std::pair<SizeType,SizeType> PairSizeType;
	typedef typename PsimagLite::Vector<RealType_>::Type VectorRealType;

	enum {TYPE_E, TYPE_A};

	enum {BRACKET_V,BRACKET_SIMPLE};

	enum {MODE_INFO = 1, MODE_ACTION = 2};

public:

	typedef RealType_ RealType;
	typedef typename BracketCollectionType::BracketType BracketType;
	typedef typename BracketType::EquationTermType EquationTermType;

	ExpansionHelper(const ModelType& model,SizeType mode)
	    : model_(model),mode_(mode),energies_(10,0),a_(10)
	{
		std::cout.precision(10);
	}

	void printHeader(std::ostream& os) const
	{
		if (!(mode_ & MODE_INFO)) return;
		os<<"-------------------\n";
		os<<"Notation V(m,n) == <m(0)|V|n(0)>\n";
		os<<"Notation E(m,n) == Em(0) -En(0)\n";
		os<<"-------------------\n";
		os<<"#FORMULAS START\n";
	}

	void push(const PsimagLite::String& str,
	                  const BracketCollectionType& bracketCollection)
	{
		if (mode_ & MODE_INFO) std::cout<<str<<bracketCollection<<"\n";
		if (!(mode_ & MODE_ACTION)) return;

		PairSizeType typeAndNumber = getTypeAndNumber(str);
		if (typeAndNumber.first == TYPE_E) {
			RealType value = getEnergy(bracketCollection);
			if (typeAndNumber.second>=energies_.size()) {
				energies_.resize(2*typeAndNumber.second+1);
			}
			energies_[typeAndNumber.second] = value;
			std::cout<<"energies["<<typeAndNumber.second<<"]="<<value<<"\n";
		} else {
			if (typeAndNumber.second>=a_.size()) {
				a_.resize(2*typeAndNumber.second+1);
			}
			VectorRealType& v = a_[typeAndNumber.second];
			v.resize(model_.total(),0);
			getState(v,bracketCollection);
//			std::cout<<"a["<<typeAndNumber.second<<"]= ";
//			const VectorRealType& vtmp = a_[typeAndNumber.second];
//			for (SizeType i = 0; i < vtmp.size(); ++i)
//				std::cout<<vtmp[i]<<" ";
//			std::cout<<"\n";
		}
	}

private:

	const RealType& a(SizeType x,SizeType y) const
	{
		assert(x < a_.size());
		assert(y < a_[x].size());
		assert(x > 0);
		return a_[x][y];
	}

	RealType getEnergy(const BracketCollectionType& bracketCollection) const
	{
		RealType value = 0;
		for (SizeType i = 0; i < bracketCollection.size(); ++i) {
			const BracketType& bracket = bracketCollection(i);
			if (bracket.factor() == 0) continue;
			value += getValue(bracket);
		}
		return value;
	}

	RealType getValue(const BracketType& bracket) const
	{
		SizeType whatBracketType = findBracketType(bracket.op());

		if (whatBracketType == BRACKET_V) return getVvalue(bracket);

		return simpleValue(bracket);
	}

	SizeType findBracketType(const PsimagLite::String& str) const
	{
		if (str.length() == 0) return BRACKET_SIMPLE;
		assert(str.length()>=6);

		if (str[0] == 'V' && str[1] == '(' && str[3] == ',' && str[5] == ')') {
			return BRACKET_V;
		}

		throw PsimagLite::RuntimeError("findBracketType");
	}

	RealType getVvalue(const BracketType& bracket) const
	{
		const PsimagLite::String& str = bracket.op();
		assert(str[0] == 'V' &&
		       str[1] == '(' &&
		       str[2] == 'n' &&
		       str[3] == ',' &&
		       str[5] == ')');

		if (str[4] == 'n') return model_.v(0,0);
		if (str[4] != 'p')
			throw PsimagLite::RuntimeError("getVvalue");

		SizeType indexOfA = getIndexOf("ap",bracket.a());

		//ap(indexOfA)V(n,p)
		RealType sum = 0.0;
		for (SizeType p = 1; p < model_.total(); p++) {
			sum += a(indexOfA,p) * model_.v(0,p) * bracket.factor();
		}
		return sum;
	}

	SizeType getIndexOf(const PsimagLite::String& match,
	                    const PsimagLite::String& str) const
	{
		if (str.length()<3 ||
		    match.length() !=2 ||
		    str[0] != match[0] ||
		    str[1] != match[1] ||
		    str[2] != '(') throw PsimagLite::RuntimeError("getIndexOf");

		return getNumberBeforeParens(str,3);
	}

	RealType simpleValue(const BracketType& bracket) const
	{
		// an(4)En(2)
		SizeType indexOfA = getIndexOf("an",bracket.a());
		SizeType indexOfE = getIndexOf("En",bracket.energy());

		return a(indexOfA,0)*energies_[indexOfE]*bracket.factor();
	}

	PairSizeType getTypeAndNumber(const PsimagLite::String& str) const
	{
		if (str.length()<3 ||
		    str[0] != 'E' ||
		    (str[1] != 'n' && str[1] != '(') ||
		    (str[2] != '(' && str[2] != 'm')) {
			throw PsimagLite::RuntimeError("getTypeAndNumber");
		}

		PsimagLite::String buffer = "";
		for (SizeType i = 3; i < str.length(); ++i) {
			if (str[i]==',') return getTypeAndNumberAux(str,i);
			if (str[i]==')') break;
			buffer += str[i];
		}

		return PairSizeType(TYPE_E,stringToNumber(buffer));

	}

	PairSizeType getTypeAndNumberAux(const PsimagLite::String& str,
	                                 SizeType ind) const
	{
		assert(str[ind] == ',');
		while (ind < str.length() && str[ind] != ')') ind++;

		if (str.length()<ind+5 ||
		    str[ind+1] != ' ' ||
		    str[ind+2] != 'a' ||
		    str[ind+3] != 'm' ||
		    str[ind+4] != '(') {
			throw PsimagLite::RuntimeError("getTypeAndNumberAux");
		}

		ind += 5;

		SizeType number = getNumberBeforeParens(str,ind);

		return PairSizeType(TYPE_A,number);
	}

	SizeType getNumberBeforeParens(const PsimagLite::String& str,
	                               SizeType ind) const
	{
		PsimagLite::String buffer = "";
		for (SizeType i = ind; i < str.length(); ++i) {
			if (str[i]==')') break;
			buffer += str[i];
		}

		return stringToNumber(buffer);
	}

	void getState(VectorRealType& v,
	              const BracketCollectionType& bracketCollection) const
	{
		for (SizeType i = 0; i < bracketCollection.size(); ++i) {
			const BracketType& bracket = bracketCollection(i);
			if (bracket.factor() == 0) continue;
			accState(v,bracket);
		}
	}

	void accState(VectorRealType& v,const BracketType& bracket) const
	{
		SizeType whatBracketType = findBracketType(bracket.op());

		return (whatBracketType == BRACKET_V) ? accVstate(v,bracket) :
		                                        accSimpleState(v,bracket);
	}

	void accVstate(VectorRealType& v, const BracketType& bracket) const
	{
		const PsimagLite::String& str = bracket.op();
		assert(str[0] == 'V' &&
		       str[1] == '(' &&
		       str[2] == 'm' &&
		       str[3] == ',' &&
		       str[5] == ')');

		if (str[4] == 'n') {
			for (SizeType m = 1; m < v.size(); m++)
				v[m] += model_.v(m,0)*bracket.factor()/model_.e(m,0);

			return;
		}
		if (str[4] != 'p')
			throw PsimagLite::RuntimeError("getVstate");

		SizeType indexOfA = getIndexOf("ap",bracket.a());

		//ap(indexOfA)V(m,p)

		for (SizeType m = 1; m < v.size(); m++) {
			RealType sum = 0.0;
			RealType tmp = bracket.factor()/model_.e(m,0);
			for (SizeType p = 0; p < model_.total(); p++) {
				sum += a(indexOfA,p) * model_.v(m,p) * tmp;
			}
			v[m] += sum;
		}
	}

	void accSimpleState(VectorRealType& v, const BracketType& bracket) const
	{
		// am(5)En(1)
		SizeType indexOfA = getIndexOf("am",bracket.a());
		SizeType indexOfE = getIndexOf("En",bracket.energy());

		RealType tmp = energies_[indexOfE]*bracket.factor();
		for (SizeType m = 1; m < v.size(); m++) {
			v[m] += a(indexOfA,m)*tmp/model_.e(m,0);
		}
	}

	SizeType stringToNumber(const PsimagLite::String& str) const
	{
		std::stringstream o;
		o<<str;
		SizeType n = 0;
		o>>n;
		return n;
	}

	const ModelType& model_;
	SizeType mode_;
	VectorRealType energies_;
	typename PsimagLite::Vector<VectorRealType>::Type a_;
};

} // namespace QuantumPerturbation

#endif // EXPANSION_HELPER_H

