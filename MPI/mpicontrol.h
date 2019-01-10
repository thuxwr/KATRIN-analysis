/*
	 This program aims to control MPI multiprocessing.

	 Weiran, Jan.8, 2019.
*/

#ifndef mpicontrol_h
#define mpicontrol_h

#include "mpi.h"
#include <iostream>

using namespace std;

class mpicontrol
{
	public:
		mpicontrol() {
			int IsInitialized;
			MPI_Initialized(&IsInitialized);
			if(!IsInitialized) MPI_Init(NULL, NULL);
			SleepFlag = 0;
			root = 0;
		}
		~mpicontrol() {}

		void Activate(int rank) {
			int currentrank;
			MPI_Comm_rank(MPI_COMM_WORLD, &currentrank);
			if(currentrank!=root) {
				cout << "Only root process can activate other processes!" << endl;
				return;
			}

			MPI_Send(&SleepFlag, 1, MPI_INT, rank, 0, MPI_COMM_WORLD);
		}

		void ActivateAll() {
			int currentrank;
			MPI_Comm_rank(MPI_COMM_WORLD, &currentrank);
			if(currentrank!=root) {
				cout << "Only root process can activate other processes!" << endl;
				return;
			}

			int totalnum;
			MPI_Comm_size(MPI_COMM_WORLD, &totalnum);
			for(int i=0; i<totalnum; i++) {
				if(i==root) continue;
				Activate(i);
			}
		}

		void SetRoot(int rank) { root = rank; }

		void Sleep() {
			int currentrank;
			MPI_Comm_rank(MPI_COMM_WORLD, &currentrank);
			if(currentrank==root) {
				cout << "Root process cannot be shut up!" << endl;
				return;
			}

			MPI_Recv(&SleepFlag, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}



	private:
		int SleepFlag; 
		int root;

};

#endif

