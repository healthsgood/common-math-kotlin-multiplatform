package optim

import exception.TooManyEvaluationsException
import exception.TooManyIterationsException
import util.Incrementor
import util.Incrementor.MaxCountExceededCallback

/**
 * Base class for implementing optimization problems. It contains the boiler-plate code
 * for counting the number of evaluations of the objective function and the number of
 * iterations of the algorithm, and storing the convergence checker.
 *
 * @param <PAIR> Type of the point/value pair returned by the optimization algorithm.
 * @author haokangkang
 * @since 3.3
</PAIR> */
abstract class AbstractOptimizationProblem<PAIR>
/**
 * Create an [AbstractOptimizationProblem] from the given data.
 *
 * @param maxEvaluations the number of allowed model function evaluations.
 * @param maxIterations  the number of allowed iterations.
 * @param checker        the convergence checker.
 */ protected constructor(
    /**
     * max evaluations
     */
    private val maxEvaluations: Int,
    /**
     * max iterations
     */
    private val maxIterations: Int,
    /**
     * Convergence checker.
     */
    override val convergenceChecker: ConvergenceChecker<PAIR>?
) : OptimizationProblem<PAIR> {
    /**
     * {@inheritDoc}
     */

    /**
     * {@inheritDoc}
     */
    override val evaluationCounter: Incrementor
        get() = Incrementor(maxEvaluations, MAX_EVAL_CALLBACK)

    /**
     * {@inheritDoc}
     */
    override val iterationCounter: Incrementor
        get() = Incrementor(maxIterations, MAX_ITER_CALLBACK)

    /**
     * Defines the action to perform when reaching the maximum number of evaluations.
     */
    private class MaxEvalCallback : MaxCountExceededCallback {
        /**
         * {@inheritDoc}
         *
         * @throws TooManyEvaluationsException
         */
        override fun trigger(max: Int) {
            throw TooManyEvaluationsException(max)
        }
    }

    /**
     * Defines the action to perform when reaching the maximum number of evaluations.
     */
    private class MaxIterCallback : MaxCountExceededCallback {
        /**
         * {@inheritDoc}
         *
         * @throws TooManyIterationsException
         */
        override fun trigger(max: Int) {
            throw TooManyIterationsException(max)
        }
    }

    companion object {
        /**
         * Callback to use for the evaluation counter.
         */
        private val MAX_EVAL_CALLBACK = MaxEvalCallback()

        /**
         * Callback to use for the iteration counter.
         */
        private val MAX_ITER_CALLBACK = MaxIterCallback()
    }
}