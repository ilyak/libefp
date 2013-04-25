TOP= ..

include $(TOP)/config.inc

check:
	@EFPMD=../efpmd/src/efpmd ./run.sh

check-omp:
	@for thr in 1 2 3; do \
		OMP_NUM_THREADS=$$thr EFPMD=../efpmd/src/efpmd ./run.sh; \
	done

check-mpi:
	@for prc in 1 2 3; do \
		EFPMD="mpirun -np $$prc ../efpmd/src/efpmd" ./run.sh; \
	done

clean:
	$(RM) *.out

.PHONY: check check-omp check-mpi clean
