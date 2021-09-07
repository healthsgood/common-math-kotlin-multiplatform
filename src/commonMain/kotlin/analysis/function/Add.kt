package analysis.function

import analysis.BivariateFunction

/**
 * Add the two operands.
 *
 * @author haokangkang
 * @since 3.0
 */
class Add : BivariateFunction {
    /**
     * {@inheritDoc}
     */
    override fun value(x: Double, y: Double): Double {
        return x + y
    }
}