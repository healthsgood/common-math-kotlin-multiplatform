package util

import util.FastMath.abs
import kotlin.jvm.JvmStatic

/**
 * Utilities for comparing numbers.
 *
 * @author haokangkang
 * @since 3.0
 */
object Precision {
    /**
     *
     *
     * Largest double-precision floating-point number such that
     * `1 + EPSILON` is numerically equal to 1. This value is an upper
     * bound on the relative error due to rounding real numbers to double
     * precision floating-point numbers.
     *
     *
     *
     * In IEEE 754 arithmetic, this is 2<sup>-53</sup>.
     *
     *
     * @see [Machine epsilon](http://en.wikipedia.org/wiki/Machine_epsilon)
     */
    const val EPSILON = 1.1102230246251565E-16

    /**
     * Safe minimum, such that `1 / SAFE_MIN` does not overflow.
     * <br></br>
     * In IEEE 754 arithmetic, this is also the smallest normalized
     * number 2<sup>-1022</sup>.
     */
    const val SAFE_MIN = 2.2250738585072014E-308


    /**
     * Compares two numbers given some amount of allowed error.
     *
     * @param x   the first number
     * @param y   the second number
     * @param eps the amount of error to allow when checking for equality
     * @return  * 0 if  [equals(x, y, eps)][.equalsByEps]
     *  * &lt; 0 if ![equals(x, y, eps)][.equalsByEps] &amp;&amp; x &lt; y
     *  * > 0 if ![equals(x, y, eps)][.equalsByEps] &amp;&amp; x > y or
     * either argument is NaN
     */
    @JvmStatic
    fun compareTo(x: Double, y: Double, eps: Double): Int {
        if (equalsByEps(x, y, eps)) {
            return 0
        } else if (x < y) {
            return -1
        }
        return 1
    }

    /**
     * Returns true iff they are equal as defined by
     * [equals(x, y, 1)][.equalsByMaxUlps].
     *
     * @param x first value
     * @param y second value
     * @return `true` if the values are equal.
     */

    fun equals(x: Double, y: Double): Boolean {
        return equalsByMaxUlps(x, y, 1)
    }

    /**
     * Returns `true` if there is no double value strictly between the
     * arguments or the difference between them is within the range of allowed
     * error (inclusive). Returns `false` if either of the arguments
     * is NaN.
     *
     * @param x   First value.
     * @param y   Second value.
     * @param eps Amount of allowed absolute error.
     * @return `true` if the values are two adjacent floating point
     * numbers or they are within range of each other.
     */

    fun equalsByEps(x: Double, y: Double, eps: Double): Boolean {
        return equalsByMaxUlps(x, y, 1) || abs(y - x) <= eps
    }

    /**
     * Returns true if the arguments are equal or within the range of allowed
     * error (inclusive).
     *
     *
     * Two float numbers are considered equal if there are `(maxUlps - 1)`
     * (or fewer) floating point numbers between them, i.e. two adjacent
     * floating point numbers are considered equal.
     *
     *
     *
     * Adapted from [
 * Bruce Dawson](http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/). Returns `false` if either of the arguments is NaN.
     *
     *
     * @param x       first value
     * @param y       second value
     * @param maxUlps `(maxUlps - 1)` is the number of floating point
     * values between `x` and `y`.
     * @return `true` if there are fewer than `maxUlps` floating
     * point values between `x` and `y`.
     */

    fun equalsByMaxUlps(x: Double, y: Double, maxUlps: Int): Boolean {
        if (maxUlps == 1) {
            error("在使用此方法之前确保使用的是精确比较")
            return x == y
        }
        throw NotImplementedError()
    }
}