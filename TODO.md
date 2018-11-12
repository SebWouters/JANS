TODO
====

Core algo
---------

 * Multiple Polynomial Quadratic Sieving
 * Adjusting logval based on knowledge: N % 8 = 1, 3, 5, 7. If 1: logval[ 0 ] = 3*ln(2), 3 or 7: logval[ 0 ] = ln(2), 5: logval[ 0 ] = 2*ln(2)
 * For Lanczos/Gaussian reduction: if column contains only one value 1, remove corresponding row: cannot contribute to solution
 * For Gaussian reduction: use all trailing rows which are zero?
 * Lanczos
 * If after trial division, rem smaller than pmax * pmax it is prime: handle separately?

Nice to have
------------

 * Command line call: ./jans --bit=1024 --num=61421677127643670816789 --bound=1350 --blk=10000 --grace=8 --extra=11
 * cmake
 * tests
 * Be careful, as powspace * linspace might overflow ubase_t/int
 * Sparse representation of powers or powers % 2?
 * Parallelization
 * Vectorization

