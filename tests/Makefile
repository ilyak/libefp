TOP= ..

include $(TOP)/config.inc

check:
	@EFPMD=../efpmd/src/efpmd ./run.sh

checkomp:
	@for thr in 1 2 3; do \
		echo "Testing for $$thr thread(s)..."; \
		OMP_NUM_THREADS=$$thr EFPMD=../efpmd/src/efpmd ./run.sh; \
	done

checkmpi:
	@for prc in 1 2 3; do \
		echo "Testing for $$prc process(es)..."; \
		OMP_NUM_THREADS=1 EFPMD="mpirun -np $$prc ../efpmd/src/efpmd" ./run.sh; \
	done

clean:
	$(RM) *.out

.PHONY: check checkomp checkmpi clean
