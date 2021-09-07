package complex

import exception.NullArgumentException
import util.FastMath.abs
import util.FastMath.sqrt
import util.MathUtils.checkNotNull
import util.MathUtils.equals
import util.MathUtils.hash
import util.Precision
import util.Precision.equalsByEps
import util.Precision.equalsByMaxUlps

/**
 * Representation of a Complex number, i.e. a number which has both a
 * real and imaginary part.
 *
 *
 * Implementations of arithmetic operations handle `NaN` and
 * infinite values according to the rules for [java.lang.Double], i.e.
 * [.equals] is an equivalence relation for all instances that have
 * a `NaN` in either real or imaginary part, e.g. the following are
 * considered equal:
 *
 *  * `1 + NaNi`
 *  * `NaN + i`
 *  * `NaN + NaNi`
 *
 *
 * Note that this contradicts the IEEE-754 standard for floating
 * point numbers (according to which the test `x == x` must fail if
 * `x` is `NaN`). The method
 * [ equals for primitive double][Precision.equalsByMaxUlps] in [Precision]
 * conforms with IEEE-754 while this class conforms with the standard behavior
 * for Java object types.
 */
class Complex constructor(
    /**
     * The real part.
     */
    val real: Double,
    /**
     * The imaginary part.
     */
    val imaginary: Double = 0.0
) {
    /**
     * Access the imaginary part.
     *
     * @return the imaginary part.
     */
    /**
     * Access the real part.
     *
     * @return the real part.
     */
    /**
     * Checks whether either or both parts of this complex number is
     * `NaN`.
     *
     * @return true if either or both parts of this complex number is
     * `NaN`; false otherwise.
     */
    /**
     * Record whether this complex number is equal to NaN.
     */

    val isNaN: Boolean
    /**
     * Checks whether either the real or imaginary part of this complex number
     * takes an infinite value (either `Double.POSITIVE_INFINITY` or
     * `Double.NEGATIVE_INFINITY`) and neither part
     * is `NaN`.
     *
     * @return true if one or both parts of this complex number are infinite
     * and neither part is `NaN`.
     */
    /**
     * Record whether this complex number is infinite.
     */

    val isInfinite: Boolean

    /**
     * Return the absolute value of this complex number.
     * Returns `NaN` if either real or imaginary part is `NaN`
     * and `Double.POSITIVE_INFINITY` if neither part is `NaN`,
     * but at least one part is infinite.
     *
     * @return the absolute value.
     */
    fun abs(): Double {
        if (isNaN) {
            return Double.NaN
        }
        if (isInfinite) {
            return Double.POSITIVE_INFINITY
        }
        return if (abs(real) < abs(imaginary)) {
            if (imaginary == 0.0) {
                return abs(real)
            }
            val q = real / imaginary
            abs(imaginary) * sqrt(1 + q * q)
        } else {
            if (real == 0.0) {
                return abs(imaginary)
            }
            val q = imaginary / real
            abs(real) * sqrt(1 + q * q)
        }
    }

    /**
     * Returns a `Complex` whose value is
     * `(this + addend)`.
     * Uses the definitional formula
     *
     *
     * `(a + bi) + (c + di) = (a+c) + (b+d)i`
     *
     * If either `this` or `addend` has a `NaN` value in
     * either part, [.NaN] is returned; otherwise `Infinite`
     * and `NaN` values are returned in the parts of the result
     * according to the rules for [java.lang.Double] arithmetic.
     *
     * @param addend Value to be added to this `Complex`.
     * @return `this + addend`.
     * @throws NullArgumentException if `addend` is `null`.
     */
    @Throws(NullArgumentException::class)
    fun add(addend: Complex): Complex {
        checkNotNull(addend)
        return if (isNaN || addend.isNaN) {
            NaN
        } else createComplex(
            real + addend.real,
            imaginary + addend.imaginary
        )
    }

    /**
     * Returns a `Complex` whose value is `(this + addend)`,
     * with `addend` interpreted as a real number.
     *
     * @param addend Value to be added to this `Complex`.
     * @return `this + addend`.
     * @see .add
     */
    fun add(addend: Double): Complex {
        return if (isNaN || addend.isNaN()) {
            NaN
        } else createComplex(real + addend, imaginary)
    }

    /**
     * Returns the conjugate of this complex number.
     * The conjugate of `a + bi` is `a - bi`.
     *
     *
     * [.NaN] is returned if either the real or imaginary
     * part of this Complex number equals `Double.NaN`.
     *
     *
     * If the imaginary part is infinite, and the real part is not
     * `NaN`, the returned value has infinite imaginary part
     * of the opposite sign, e.g. the conjugate of
     * `1 + POSITIVE_INFINITY i` is `1 - NEGATIVE_INFINITY i`.
     *
     *
     * @return the conjugate of this Complex object.
     */
    fun conjugate(): Complex {
        return if (isNaN) {
            NaN
        } else createComplex(real, -imaginary)
    }

    /**
     * Returns a `Complex` whose value is
     * `(this / divisor)`.
     * Implements the definitional formula
     * <pre>
     * `
     * a + bi          ac + bd + (bc - ad)i
     * ----------- = -------------------------
     * c + di         c<sup>2</sup> + d<sup>2</sup>
    ` *
    </pre> *
     * but uses
     * [
 * prescaling of operands](http://doi.acm.org/10.1145/1039813.1039814) to limit the effects of overflows and
     * underflows in the computation.
     *
     *
     * `Infinite` and `NaN` values are handled according to the
     * following rules, applied in the order presented:
     *
     *  * If either `this` or `divisor` has a `NaN` value
     * in either part, [.NaN] is returned.
     *
     *  * If `divisor` equals [.ZERO], [.NaN] is returned.
     *
     *  * If `this` and `divisor` are both infinite,
     * [.NaN] is returned.
     *
     *  * If `this` is finite (i.e., has no `Infinite` or
     * `NaN` parts) and `divisor` is infinite (one or both parts
     * infinite), [.ZERO] is returned.
     *
     *  * If `this` is infinite and `divisor` is finite,
     * `NaN` values are returned in the parts of the result if the
     * [java.lang.Double] rules applied to the definitional formula
     * force `NaN` results.
     *
     *
     *
     * @param divisor Value by which this `Complex` is to be divided.
     * @return `this / divisor`.
     * @throws NullArgumentException if `divisor` is `null`.
     */
    @Throws(NullArgumentException::class)
    fun divide(divisor: Complex): Complex {
        checkNotNull(divisor)
        if (isNaN || divisor.isNaN) {
            return NaN
        }
        val c = divisor.real
        val d = divisor.imaginary
        if (c == 0.0 && d == 0.0) {
            return NaN
        }
        if (divisor.isInfinite && !isInfinite) {
            return ZERO
        }
        return if (abs(c) < abs(d)) {
            val q = c / d
            val denominator = c * q + d
            createComplex(
                (real * q + imaginary) / denominator,
                (imaginary * q - real) / denominator
            )
        } else {
            val q = d / c
            val denominator = d * q + c
            createComplex(
                (imaginary * q + real) / denominator,
                (imaginary - real * q) / denominator
            )
        }
    }

    /**
     * Test for equality with another object.
     * If both the real and imaginary parts of two complex numbers
     * are exactly the same, and neither is `Double.NaN`, the two
     * Complex objects are considered to be equal.
     * The behavior is the same as for JDK's [ Double][Double.equals]:
     *
     *  * All `NaN` values are considered to be equal,
     * i.e, if either (or both) real and imaginary parts of the complex
     * number are equal to `Double.NaN`, the complex number is equal
     * to `NaN`.
     *
     *  *
     * Instances constructed with different representations of zero (i.e.
     * either "0" or "-0") are *not* considered to be equal.
     *
     *
     *
     * @param other Object to test for equality with this instance.
     * @return `true` if the objects are equal, `false` if object
     * is `null`, not an instance of `Complex`, or not equal to
     * this instance.
     */
    override fun equals(other: Any?): Boolean {
        if (this === other) {
            return true
        }
        if (other is Complex) {
            val c = other
            return if (c.isNaN) {
                isNaN
            } else {
                equals(real, c.real) &&
                        equals(imaginary, c.imaginary)
            }
        }
        return false
    }

    /**
     * Get a hashCode for the complex number.
     * Any `Double.NaN` value in real or imaginary part produces
     * the same hash code `7`.
     *
     * @return a hash code value for this object.
     */
    override fun hashCode(): Int {
        return if (isNaN) {
            7
        } else 37 * (17 * hash(imaginary) +
                hash(real))
    }

    /**
     * Returns a `Complex` whose value is `(-this)`.
     * Returns `NaN` if either real or imaginary
     * part of this Complex number is `Double.NaN`.
     *
     * @return `-this`.
     */
    fun negate(): Complex {
        return if (isNaN) {
            NaN
        } else createComplex(-real, -imaginary)
    }

    /**
     * Returns a `Complex` whose value is
     * `(this - subtrahend)`.
     * Uses the definitional formula
     *
     *
     * `(a + bi) - (c + di) = (a-c) + (b-d)i`
     *
     * If either `this` or `subtrahend` has a `NaN]` value in either part,
     * [.NaN] is returned; otherwise infinite and `NaN` values are
     * returned in the parts of the result according to the rules for
     * [java.lang.Double] arithmetic.
     *
     * @param subtrahend value to be subtracted from this `Complex`.
     * @return `this - subtrahend`.
     * @throws NullArgumentException if `subtrahend` is `null`.
     */
    @Throws(NullArgumentException::class)
    fun subtract(subtrahend: Complex): Complex {
        checkNotNull(subtrahend)
        return if (isNaN || subtrahend.isNaN) {
            NaN
        } else createComplex(
            real - subtrahend.real,
            imaginary - subtrahend.imaginary
        )
    }

    /**
     * Returns a `Complex` whose value is
     * `(this - subtrahend)`.
     *
     * @param subtrahend value to be subtracted from this `Complex`.
     * @return `this - subtrahend`.
     * @see .subtract
     */
    fun subtract(subtrahend: Double): Complex {
        return if (isNaN || subtrahend.isNaN()) {
            NaN
        } else createComplex(real - subtrahend, imaginary)
    }

    /**
     * Create a complex number given the real and imaginary parts.
     *
     * @param realPart      Real part.
     * @param imaginaryPart Imaginary part.
     * @return a new complex number instance.
     * @see .valueOf
     * @since 1.2
     */
    protected fun createComplex(
        realPart: Double,
        imaginaryPart: Double
    ): Complex {
        return Complex(realPart, imaginaryPart)
    }

    /**
     * {@inheritDoc}
     */
    override fun toString(): String {
        return "($real, $imaginary)"
    }

    companion object {
        /**
         * The square root of -1. A number representing "0.0 + 1.0i"
         */
        val I = Complex(0.0, 1.0)
        // CHECKSTYLE: stop ConstantName
        /**
         * A complex number representing "NaN + NaNi"
         */
        val NaN = Complex(Double.NaN, Double.NaN)
        // CHECKSTYLE: resume ConstantName
        /**
         * A complex number representing "+INF + INFi"
         */
        val INF = Complex(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY)

        /**
         * A complex number representing "1.0 + 0.0i"
         */
        val ONE = Complex(1.0, 0.0)

        /**
         * A complex number representing "0.0 + 0.0i"
         */
        val ZERO = Complex(0.0, 0.0)
        /**
         * Test for the floating-point equality between Complex objects.
         * It returns `true` if both arguments are equal or within the
         * range of allowed error (inclusive).
         *
         * @param x       First value (cannot be `null`).
         * @param y       Second value (cannot be `null`).
         * @param maxUlps `(maxUlps - 1)` is the number of floating point
         * values between the real (resp. imaginary) parts of `x` and
         * `y`.
         * @return `true` if there are fewer than `maxUlps` floating
         * point values between the real (resp. imaginary) parts of `x`
         * and `y`.
         * @see Precision.equalsByMaxUlps
         * @since 3.3
         */
        /**
         * Returns `true` iff the values are equal as defined by
         * [equals(x, y, 1)][.equals].
         *
         * @param x First value (cannot be `null`).
         * @param y Second value (cannot be `null`).
         * @return `true` if the values are equal.
         * @since 3.3
         */

        fun equals(x: Complex, y: Complex, maxUlps: Int = 1): Boolean {
            return equalsByMaxUlps(x.real, y.real, maxUlps) &&
                    equalsByMaxUlps(x.imaginary, y.imaginary, maxUlps)
        }

        /**
         * Returns `true` if, both for the real part and for the imaginary
         * part, there is no double value strictly between the arguments or the
         * difference between them is within the range of allowed error
         * (inclusive).  Returns `false` if either of the arguments is NaN.
         *
         * @param x   First value (cannot be `null`).
         * @param y   Second value (cannot be `null`).
         * @param eps Amount of allowed absolute error.
         * @return `true` if the values are two adjacent floating point
         * numbers or they are within range of each other.
         * @see Precision.equalsByEps
         * @since 3.3
         */
        fun equals(x: Complex, y: Complex, eps: Double): Boolean {
            return equalsByEps(x.real, y.real, eps) &&
                    equalsByEps(x.imaginary, y.imaginary, eps)
        }

        /**
         * Create a complex number given the real and imaginary parts.
         *
         * @param realPart      Real part.
         * @param imaginaryPart Imaginary part.
         * @return a Complex instance.
         */
        fun valueOf(
            realPart: Double,
            imaginaryPart: Double
        ): Complex {
            return if (realPart.isNaN() ||
                imaginaryPart.isNaN()
            ) {
                NaN
            } else Complex(realPart, imaginaryPart)
        }

        /**
         * Create a complex number given only the real part.
         *
         * @param realPart Real part.
         * @return a Complex instance.
         */
        fun valueOf(realPart: Double): Complex {
            return if (realPart.isNaN()) {
                NaN
            } else Complex(realPart)
        }
    }
    /**
     * Create a complex number given the real and imaginary parts.
     *
     * @param real      Real part.
     * @param imaginary Imaginary part.
     */
    /**
     * Create a complex number given only the real part.
     *
     * @param real Real part.
     */
    init {
        isNaN = real.isNaN() || imaginary.isNaN()
        isInfinite = !isNaN &&
                real.isInfinite() || imaginary.isInfinite()
    }
}