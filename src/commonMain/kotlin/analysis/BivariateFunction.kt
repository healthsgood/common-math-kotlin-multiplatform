package analysis

/**
 * An interface representing a bivariate real function.
 *
 * @author haokangkang
 * @since 2.1
 */
interface BivariateFunction {
    /**
     * Compute the value for the function.
     *
     * @param x Abscissa for which the function value should be computed.
     * @param y Ordinate for which the function value should be computed.
     * @return the value.
     */
    fun value(x: Double, y: Double): Double
}