#g++ -g
g++ -O3 -march=native -flto\
    src/main.cpp\
    src/big_int.cpp\
    src/big_int_math.cpp\
    src/big_int_io.cpp\
    src/big_int_private.cpp\
    src/sieve.cpp\
    src/sieve_startup.cpp\
    src/sieve_numtheo.cpp -o jans

