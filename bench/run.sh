#!bin/bash

../build/bench/workload Bloom mixed 50 0 randint point 0 zipfian

../build/bench/workload SuRF mixed 50 0 randint point 0 zipfian

../build/bench/workload SuRFHash mixed 50 0 randint point 0 zipfian

../build/bench/workload SuRFReal mixed 50 0 randint point 0 zipfian


../build/bench/workload Bloom mixed 50 0 timestamp point 0 zipfian

../build/bench/workload SuRF mixed 50 0 timestamp point 0 zipfian

../build/bench/workload SuRFHash mixed 50 0 timestamp point 0 zipfian

../build/bench/workload SuRFReal mixed 50 0 timestamp point 0 zipfian


../build/bench/workload Bloom mixed 50 0 email point 0 zipfian

../build/bench/workload SuRF mixed 50 0 email point 0 zipfian

../build/bench/workload SuRFHash mixed 50 0 email point 0 zipfian

../build/bench/workload SuRFReal mixed 50 0 email point 0 zipfian


#../build/bench/workload SuRF mixed 50 0 randint range 352187318272 zipfian

#../build/bench/workload SuRFHash mixed 50 0 randint range 352187318272 zipfian

#../build/bench/workload SuRFReal mixed 50 0 randint range 352187318272 zipfian


#../build/bench/workload SuRF mixed 50 0 timestamp range 1000000000 zipfian

#../build/bench/workload SuRFHash mixed 50 0 timestamp range 1000000000 zipfian

#../build/bench/workload SuRFReal mixed 50 0 timestamp range 1000000000 zipfian
