#!/bin/sh

for TEST in `basename -s .in *.in`; do
	$EFPMD $TEST.in > $TEST.out

	if grep -q "DOES NOT MATCH" $TEST.out; then
		echo -e "\033[33;31mFAILURE: $TEST\033[0m"
	else
		echo -e "\033[33;32mSUCCESS: $TEST\033[0m"
	fi
done
