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
#include <iostream>
#include <unistd.h>
#include <vector>
#include <utility>
#include <cassert>
#include <cstdlib>
#include <algorithm>
#include <stdexcept>
#include "TypeToString.h"
#include "Models/HubbardOneBand/HubbardOneBandExpandU.h"
#include "Models/HubbardOneBand/HubbardOneBandExpandHopping.h"
#include "Expansion.h"
#include "../../FreeFermions/src/GeometryParameters.h"
#include "Concurrency.h"
#include "Provenance.h"

const PsimagLite::String license=
"Copyright (c) 2015-2016, UT-Battelle, LLC\n"
"All rights reserved\n"
"\n"
"[qpt++, Version 0.]\n"
"\n"
"---------------------------------------------------------\n"
"THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND\n"
"CONTRIBUTORS \"AS IS\" AND ANY EXPRESS OR IMPLIED\n"
"WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED\n"
"WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A\n"
"PARTICULAR PURPOSE ARE DISCLAIMED.\n"
"\n"
"Please see full open source license included in file LICENSE.\n"
"---------------------------------------------------------\n"
"\n";

using namespace QuantumPerturbation;

typedef double RealType;
typedef PsimagLite::InputNg<FreeFermions::InputCheck> InputNgType;
typedef FreeFermions::GeometryParameters<RealType,InputNgType::Readable> GeometryParamsType;

template<typename ModelType>
void mainLoop(const PsimagLite::String& file,SizeType n, SizeType mode)
{
	FreeFermions::InputCheck inputCheck;
	InputNgType::Writeable ioWriteable(file,inputCheck);
	InputNgType::Readable io(ioWriteable);

	ModelType model(io);

	Expansion<ModelType> expansion(model,mode,n);

	expansion.expand();
}

int main(int argc,char *argv[])
{
	SizeType n = 0;
	int opt = 0;
	PsimagLite::String modelName = "dummy";
	PsimagLite::String file("");
	SizeType mode = 0;
	PsimagLite::String expand("");
	bool versionOnly = false;

	while ((opt = getopt(argc, argv, "n:f:m:e:V")) != -1) {
		switch (opt) {
		case 'n':
			n = atoi(optarg);
			break;
		case 'f':
			file=optarg;
			break;
		case 'm':
			mode=atoi(optarg);
			break;
		case 'e':
			expand = optarg;
			break;
		case 'V':
			versionOnly = true;
			break;
		default:
			std::cerr<<"USAGE is "<<argv[0]<<" -f file -e expansion -n n [-m mode] | -V\n";
			return 1;
		}
	}

	bool b = (n == 0 || file == "" || expand == "");
	if (b && !versionOnly) {
		std::cerr<<"USAGE is "<<argv[0]<<" -f file -e expansion -n n [-m mode]\n";
		return 1;
	}

	size_t npthreads = 1;
	PsimagLite::Concurrency concurrency(&argc,&argv,npthreads);

	// print license
	if (PsimagLite::Concurrency::root()) {
		std::cout<<license;
		Provenance provenance;
		std::cout<<provenance;
	}

	if (versionOnly) return 0;

	GeometryParamsType::readLabel(modelName,file,"Model=");

	modelName += "Expand";
	modelName += expand;

	if (modelName == "HubbardOneBandExpandU") {
		mainLoop<HubbardOneBandExpandU<RealType> >(file,n,mode);
	} else if (modelName == "HubbardOneBandExpandHopping") {
		mainLoop<HubbardOneBandExpandHopping<RealType> >(file,n,mode);
	} else {
		throw PsimagLite::RuntimeError("Unknown model");
	}
}

