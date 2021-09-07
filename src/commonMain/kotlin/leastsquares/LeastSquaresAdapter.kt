package leastsquares

import leastsquares.LeastSquaresProblem.Evaluation
import linear.RealVector
import optim.ConvergenceChecker
import util.Incrementor

/**
 * An adapter that delegates to another implementation of [LeastSquaresProblem].
 * //todo 此类没有意义，考虑删除
 *
 * @author haokangkang
 * @since 3.3
 */
open class LeastSquaresAdapter
/**
 * Delegate the [LeastSquaresProblem] interface to the given implementation.
 *
 * @param problem the delegate
 */(
    /**
     * the delegate problem
     */
    private val problem: LeastSquaresProblem
) : LeastSquaresProblem {
    /**
     * {@inheritDoc}
     */
    override fun getStart(): RealVector? {
        return problem.getStart()
    }

    /**
     * {@inheritDoc}
     */
    override fun getObservationSize(): Int {
        return problem.getObservationSize()
    }

    /**
     * {@inheritDoc}
     */
    override fun getParameterSize(): Int {
        return problem.getParameterSize()
    }

    /**
     * {@inheritDoc}
     *
     * @param point
     */
    override fun evaluate(point: RealVector?): Evaluation? {
        return problem.evaluate(point)
    }

    /**
     * {@inheritDoc}
     */
    override val evaluationCounter: Incrementor
        get() = problem.evaluationCounter

    /**
     * {@inheritDoc}
     */
    override val iterationCounter: Incrementor
        get() = problem.iterationCounter

    /**
     * {@inheritDoc}
     */
    override val convergenceChecker: ConvergenceChecker<Evaluation?>?
        get() = problem.convergenceChecker
}