package util

import exception.NullArgumentException

/**
 * Miscellaneous utility functions.
 *
 * @author haokangkang
 * @see Precision
 */
object MathUtils {
    /**
     * Returns an integer hash code representing the given double value.
     *
     * @param value the value to be hashed
     * @return the hash code
     */

    fun hash(value: Double): Int {
        return value.hashCode()
    }

    /**
     * Returns `true` if the values are equal according to semantics of
     * [Double.equals].
     *
     * @param x Value
     * @param y Value
     * @return `new Double(x).equals(new Double(y))`
     */

    fun equals(x: Double, y: Double): Boolean {
        return x == y
    }

    /**
     * Returns an integer hash code representing the given double array.
     *
     * @param value the value to be hashed (may be null)
     * @return the hash code
     * @since 1.2
     */

    fun hash(value: DoubleArray?): Int {
        return value.hashCode()
    }

    /**
     * Checks that an object is not null.
     *
     * @param o Object to be checked.
     * @throws NullArgumentException if `o` is `null`.
     */

    @Throws(NullArgumentException::class)
    fun checkNotNull(o: Any?) {
        if (o == null) {
            throw NullArgumentException()
        }
    }
}