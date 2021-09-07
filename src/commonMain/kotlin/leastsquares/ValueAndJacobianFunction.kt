package leastsquares

import linear.RealMatrix
import linear.RealVector

/**
 * A interface for functions that compute a vector of values and can compute their
 * derivatives (Jacobian).
 *
 * @author haokangkang
 * @since 3.4
 */
interface ValueAndJacobianFunction : MultivariateJacobianFunction {
    /**
     * Compute the value.
     *
     * @param params Point.
     * @return the value at the given point.
     */
    fun computeValue(params: DoubleArray?): RealVector?

    /**
     * Compute the Jacobian.
     *
     * @param params Point.
     * @return the Jacobian at the given point.
     */
    fun computeJacobian(params: DoubleArray?): RealMatrix?
}