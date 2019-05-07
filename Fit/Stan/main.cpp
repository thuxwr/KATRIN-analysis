#include <cmdstan/command.hpp>
#include <stan/services/error_codes.hpp>
#include <boost/exception/diagnostic_information.hpp> 
#include <boost/exception_ptr.hpp>
#include "model.hpp"
#include "Pred.h"

int main(int argc, const char* argv[]) {
	char* Argv = new char[40];
	strcpy(Argv, argv[0]);
	detect.initialize(argc, &Argv);
  try {
    return cmdstan::command<stan_model>(argc,argv);
  } catch (const std::exception& e) {
    std::cout << e.what() << std::endl;
    return stan::services::error_codes::SOFTWARE;
  }
}
