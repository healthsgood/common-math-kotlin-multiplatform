package leastsquares

import linear.RealMatrix
import linear.RealVector
import util.Pair

/**
 * A interface for functions that compute a vector of values and can compute their
 * derivatives (Jacobian).
 *
 * @author haokangkang
 * @since 3.3
 */
interface MultivariateJacobianFunction {
    /**
     * Compute the function value and its Jacobian.
     *
     * @param point the abscissae
     * @return the values and their Jacobian of this vector valued function.
     */
    fun value(point: RealVector?): Pair<RealVector?, RealMatrix?>?
}