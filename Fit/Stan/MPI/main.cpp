#include <cmdstan/command.hpp>
#include <stan/services/error_codes.hpp>
#include <boost/exception/diagnostic_information.hpp>
#include <boost/exception_ptr.hpp>
#include "mpi.h"
#include "model.hpp"
#include "Pred.h"
#include <stdio.h>

int main(int argc, const char* argv[]) {
	int initflag;
	MPI_Initialized(&initflag);
	if(!initflag) MPI_Init(NULL, NULL);
	char* Argv = new char[40];
	strcpy(Argv, argv[0]);
	detect.initialize(argc, &Argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if(rank==0) { //root process
		try {
			return cmdstan::command<stan_model>(argc,argv);
		} catch (const std::exception& e) {
			std::cout << e.what() << std::endl;
			return stan::services::error_codes::SOFTWARE;
		}
	}
	else {
		while(true) {
			detect.AssociateRun();
		}
	}
	return 0;
}

