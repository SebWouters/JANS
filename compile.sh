#g++ -O3 -march=native -flto\
g++ -g    main.cpp\
    big_int.cpp\
    big_int_math.cpp\
    big_int_io.cpp\
    big_int_private.cpp\
    sieve.cpp -o jans
