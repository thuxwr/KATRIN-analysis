#include <cmdstan/command.hpp>
#include <stan/services/error_codes.hpp>
#include <boost/exception/diagnostic_information.hpp>
#include <boost/exception_ptr.hpp>
#include "mpi.h"
#include "model.hpp"
#include "Pred.h"
#include <stdio.h>

int main(int argc, const char* argv[]) {
	MPI_Init(NULL, NULL);
	char* Argv = new char[40];
	strcpy(Argv, argv[0]);
	detect.initialize(argc, &Argv);
	int rank, size;
	double* pars = new double[10];
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
			MPI_Recv(pars, 10, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			detect.SetScatParams(pars[0], pars[1], pars[2], pars[3], pars[4], pars[5], pars[6]);
			detect.SetupResponse(pars[7], pars[8], pars[9]);
		}
	}
	return 0;
}

