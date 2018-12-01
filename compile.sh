#g++ -g
#g++ -O3 -march=native -flto -fopenmp -funroll-loops\
icpc -flto -xHost -qopenmp -ipo -O3 -Wall\
    src/executable.cpp\
    src/big_int.cpp\
    src/big_int_math.cpp\
    src/big_int_io.cpp\
    src/big_int_private.cpp\
    src/sieve.cpp\
    src/sieve_startup.cpp\
    src/sieve_numtheo.cpp\
    src/sieve_gaussian.cpp -o jans

