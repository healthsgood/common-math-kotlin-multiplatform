package leastsquares

import leastsquares.LeastSquaresOptimizer.Optimum
import leastsquares.LeastSquaresProblem.Evaluation
import linear.RealMatrix
import linear.RealVector

/**
 * A pedantic implementation of [Optimum].
 *
 * @since 3.3
 */
internal class OptimumImpl
/**
 * Construct an optimum from an evaluation and the values of the counters.
 *
 * @param value       the function value
 * @param evaluations number of times the function was evaluated
 * @param iterations  number of iterations of the algorithm
 */(
    /**
     * abscissa and ordinate
     */
    private val value: Evaluation,
    /**
     * number of evaluations to compute this optimum
     */
    private val evaluations: Int,
    /**
     * number of iterations to compute this optimum
     */
    private val iterations: Int
) : Optimum {
    /* auto-generated implementations */
    /**
     * {@inheritDoc}
     */
    override fun getEvaluations(): Int {
        return evaluations
    }

    /**
     * {@inheritDoc}
     */
    override fun getIterations(): Int {
        return iterations
    }

    /**
     * {@inheritDoc}
     */
    override fun getCovariances(threshold: Double): RealMatrix? {
        return value.getCovariances(threshold)
    }

    /**
     * {@inheritDoc}
     */
    override fun getSigma(covarianceSingularityThreshold: Double): RealVector? {
        return value.getSigma(covarianceSingularityThreshold)
    }

    /**
     * {@inheritDoc}
     */
    override fun getJacobian(): RealMatrix? {
        return value.getJacobian()
    }

    /**
     * {@inheritDoc}
     */
    override fun getCost(): Double {
        return value.getCost()
    }

    /**
     * {@inheritDoc}
     */
    override fun getResiduals(): RealVector? {
        return value.getResiduals()
    }

    /**
     * {@inheritDoc}
     */
    override fun getPoint(): RealVector? {
        return value.getPoint()
    }
}