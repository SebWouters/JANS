#g++ -g
#g++ -O3 -march=native -flto -fopenmp -funroll-loops\
#icpc -flto -xHost -qopenmp -ipo -O3 -Wall\
#g++ -O3 -march=native -flto -fopenmp -funroll-loops\
g++ -g\
    src/executable.cpp\
    src/big_int.cpp\
    src/big_int_math.cpp\
    src/big_int_io.cpp\
    src/big_int_private.cpp\
    src/sieve.cpp\
    src/sieve_startup.cpp\
    src/sieve_numtheo.cpp\
    src/solver/clean.cpp\
    src/solver/power_contributions.cpp\
    src/solver/gaussian.cpp -o jans

