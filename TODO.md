TODO
====

 * Sieve negative offsets Q(x) = x*x - n with x < sqrt(n), with extra leading bit in powers
 * MPQS
 * Adjusting logval based on knowledge: N % 8 = 1, 3, 5, 7. If 1: logval[ 0 ] = 3*ln(2), 3 or 7: logval[ 0 ] = ln(2), 5: logval[ 0 ] = 2*ln(2)
 * For Gaussian reduction: if column contains only one value 1, remove corresponding row: cannot contribute to solution
 * For Gaussian reduction: use all trailing rows which are zero?
 * Lanczos
 * If after trial division, rem < p_max^2 it is prime: handle separately?
 * Be careful, as linspace * linspace might overflow ubase_t/int

