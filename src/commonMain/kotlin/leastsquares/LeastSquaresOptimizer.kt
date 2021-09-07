package leastsquares

import leastsquares.LeastSquaresProblem.Evaluation

/**
 * An algorithm that can be applied to a non-linear least squares problem.
 *
 * @author haokangkang
 * @since 3.3
 */
interface LeastSquaresOptimizer {
    /**
     * Solve the non-linear least squares problem.
     *
     * @param leastSquaresProblem the problem definition, including model function and
     * convergence criteria.
     * @return The optimum.
     */
    fun optimize(leastSquaresProblem: LeastSquaresProblem?): Optimum?

    /**
     * The optimum found by the optimizer. This object contains the point, its value, and
     * some metadata.
     */
    interface Optimum : Evaluation {
        /**
         * Get the number of times the model was evaluated in order to produce this
         * optimum.
         *
         * @return the number of model (objective) function evaluations
         */
        fun getEvaluations(): Int

        /**
         * Get the number of times the algorithm iterated in order to produce this
         * optimum. In general least squares it is common to have one [ ][.getEvaluations] per iterations.
         *
         * @return the number of iterations
         */
        fun getIterations(): Int
    }
}