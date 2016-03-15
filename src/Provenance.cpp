#include "Provenance.h"

std::ostream& operator<<(std::ostream& os,const Provenance&)
{
	os<<"qpt++ version "<<QPTPP_VERSION<<"\n";
	os<<"PsimagLite version "<<PSIMAGLITE_VERSION<<"\n";
	return os;
}

