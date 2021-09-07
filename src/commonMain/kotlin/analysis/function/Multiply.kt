package analysis.function

import analysis.BivariateFunction

/**
 * Multiply the two operands.
 *
 * @author haokangkang
 * @since 3.0
 */
class Multiply : BivariateFunction {
    /**
     * {@inheritDoc}
     */
    override fun value(x: Double, y: Double): Double {
        return x * y
    }
}