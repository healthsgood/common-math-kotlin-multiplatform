package optim

import util.Incrementor

/**
 * Common settings for all optimization problems. Includes divergence and convergence
 * criteria.
 *
 * @param <PAIR> The type of value the [convergence][.getConvergenceChecker] will operate on. It should include the value of the model
 * function and point where it was evaluated.
 * @author haokangkang
 * @since 3.3
</PAIR> */
interface OptimizationProblem<PAIR> {
    /**
     * Get a independent Incrementor that counts up to the maximum number of evaluations
     * and then throws an exception.
     *
     * @return a counter for the evaluations.
     */
    val evaluationCounter: Incrementor

    /**
     * Get a independent Incrementor that counts up to the maximum number of iterations
     * and then throws an exception.
     *
     * @return a counter for the evaluations.
     */
    val iterationCounter: Incrementor

    /**
     * Gets the convergence checker.
     *
     * @return the object used to check for convergence.
     */
    val convergenceChecker: ConvergenceChecker<PAIR>?
}