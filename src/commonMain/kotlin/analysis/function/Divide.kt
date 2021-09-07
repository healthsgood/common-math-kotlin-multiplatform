package analysis.function

import analysis.BivariateFunction

/**
 * Divide the first operand by the second.
 *
 * @author haokangkang
 * @since 3.0
 */
class Divide : BivariateFunction {
    /**
     * {@inheritDoc}
     */
    override fun value(x: Double, y: Double): Double {
        return x / y
    }
}