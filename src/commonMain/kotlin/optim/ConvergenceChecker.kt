package optim

/**
 * This interface specifies how to check if an optimization algorithm has
 * converged.
 * <br></br>
 * Deciding if convergence has been reached is a problem-dependent issue. The
 * user should provide a class implementing this interface to allow the
 * optimization algorithm to stop its search according to the problem at hand.
 * <br></br>
 *
 * @author haokangkang
 * @since 3.0
 */
interface ConvergenceChecker<PAIR> {
    /**
     * Check if the optimization algorithm has converged.
     *
     * @param iteration Current iteration.
     * @param previous  Best point in the previous iteration.
     * @param current   Best point in the current iteration.
     * @return `true` if the algorithm is considered to have converged.
     */
    fun converged(iteration: Int, previous: PAIR, current: PAIR): Boolean
}