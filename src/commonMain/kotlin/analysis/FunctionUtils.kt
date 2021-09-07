package analysis


/**
 * Utilities for manipulating function objects.
 *
 * @author haokangkang
 * @since 3.0
 */
object FunctionUtils {
    /**
     * Creates a unary function by fixing the second argument of a binary function.
     *
     * @param f     Binary function.
     * @param fixed value to which the second argument of `f` is set.
     * @return the unary function h(x) = f(x, fixed)
     */

    fun fix2ndArgument(
        f: BivariateFunction,
        fixed: Double
    ): UnivariateFunction {
        return object : UnivariateFunction {
            /** {@inheritDoc}  */
            override fun value(x: Double): Double {
                return f.value(x, fixed)
            }
        }
    }
}