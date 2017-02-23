from mpi4py import MPI
import models
import cobra.test
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
batch_size = 1

