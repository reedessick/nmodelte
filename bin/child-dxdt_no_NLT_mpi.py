#!/usr/local/bin/python

from mpi4py import MPI
import numpy as np
from network_flow import __dxmdt_no_NLT_mpi

### figure out your MPI world
comm = MPI.Comm.Get_parent()
size = comm.Get_size()
rank = comm.Get_rank()

### collect basic info from parent
dimension = np.array(0, dtype="i")
comm.Bcast([dimension, MPI.INT], root=0)
#print "child", rank, dimension

Oorb = np.array(0, dtype="d")
comm.Bcast([Oorb, MPI.INT], root=0)
#print "child", rank, Oorb

nlm = np.empty((dimension/2, 3), dtype="i")
comm.Bcast([nlm, MPI.INT], root=0)
#print "child", rank, nlm
#print nlm[0][-1]

wyU = comm.bcast(root=0) # python lists
#print "child", rank, wyU

K   = comm.bcast(root=0)
#print "child", rank, K

### figure out best division of labor ==> compute Msets
num_k_par = int(np.ceil(1.0*sum([len(k) for k in K]) / size))
setNo = 0
n_k = 0
Mset = []
for modeNo, k in enumerate(K):
  n_k += len(k)
  if setNo == rank:
    Mset.append(modeNo) # we only care about THIS instance
  if n_k >= num_k_par:
    n_k = 0
    setNo += 1
    if setNo > rank: # we only care about THIS instance
      break

### define vectors that we share with parent
t = np.array(0.0, dtype="f")
xvec = np.empty(dimension, dtype="f")
dxmdt = np.zeros(dimension, dtype="f")

### define persistent communication channels with parent
rcvt = comm.Recv_init([t, MPI.FLOAT], source=0)#, tag=rank)
rcvx = comm.Recv_init([xvec, MPI.FLOAT], source=0)#, tag=size+rank)
snd = comm.Send_init([dxmdt, MPI.FLOAT], dest=0)#, tag=rank)

### always listen for new data and process it as it arrives
### this loop is terminated when the process dies (caused by parent's termination)
while True:
  rcvt.Start() # request t from parent
  rcvt.Wait() # wait for t
  rcvx.Start() # request xvec from parent
  rcvx.Wait() # wait for xvec

#  print "child", rank, xvec
  __dxmdt_no_NLT_mpi(dxmdt, Mset, t, xvec, Oorb, wyU, nlm, K) ### must modify array in place!
#  print "child", rank, dxmdt

  snd.Start() # send dxmdt to parent
              # we don't need to wait because nothing will be sent until this is received

