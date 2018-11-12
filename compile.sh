#g++ -O3 -march=native -flto\
g++ -g\
    src/main.cpp\
    src/big_int.cpp\
    src/big_int_math.cpp\
    src/big_int_io.cpp\
    src/big_int_private.cpp\
    src/sieve.cpp -o jans

