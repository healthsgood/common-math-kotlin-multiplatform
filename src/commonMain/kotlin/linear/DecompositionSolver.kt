package linear

import linear.exception.SingularMatrixException

/**
 * Interface handling decomposition algorithms that can solve A  X = B.
 *
 *
 * Decomposition algorithms decompose an A matrix has a product of several specific
 * matrices from which they can solve A  X = B in least squares sense: they find X
 * such that ||A  X - B|| is minimal.
 *
 *
 * Some solvers like [LUDecomposition] can only find the solution for
 * square matrices and when the solution is an exact linear solution, i.e. when
 * ||A  X - B|| is exactly 0. Other solvers can also find solutions
 * with non-square matrix A and with non-null minimal norm. If an exact linear
 * solution exists it is also the minimal norm solution.
 *
 * @author haokangkang
 * @since 2.0
 */
interface DecompositionSolver {
    /**
     * Solve the linear equation A  X = B for matrices A.
     *
     *
     * The A matrix is implicit, it is provided by the underlying
     * decomposition algorithm.
     *
     * @param b right-hand side of the equation A  X = B
     * @return a vector X that minimizes the two norm of A  X - B
     * @throws exception.DimensionMismatchException if the matrices dimensions do not match.
     * @throws SingularMatrixException                                       if the decomposed matrix is singular.
     */
    @Throws(SingularMatrixException::class)
    fun solve(b: RealVector?): RealVector?

    /**
     * Solve the linear equation A  X = B for matrices A.
     *
     *
     * The A matrix is implicit, it is provided by the underlying
     * decomposition algorithm.
     *
     * @param b right-hand side of the equation A  X = B
     * @return a matrix X that minimizes the two norm of A  X - B
     * @throws exception.DimensionMismatchException if the matrices dimensions do not match.
     * @throws SingularMatrixException                                       if the decomposed matrix is singular.
     */
    @Throws(SingularMatrixException::class)
    fun solve(b: RealMatrix?): RealMatrix?

    /**
     * Check if the decomposed matrix is non-singular.
     *
     * @return true if the decomposed matrix is non-singular.
     */
    val isNonSingular: Boolean

    /**
     * Get the [pseudo-inverse](http://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_pseudoinverse)
     * of the decomposed matrix.
     *
     *
     * *This is equal to the inverse  of the decomposed matrix, if such an inverse exists.*
     *
     *
     * If no such inverse exists, then the result has properties that resemble that of an inverse.
     *
     *
     * In particular, in this case, if the decomposed matrix is A, then the system of equations
     * \( A x = b \) may have no solutions, or many. If it has no solutions, then the pseudo-inverse
     * \( A^+ \) gives the "closest" solution \( z = A^+ b \), meaning \( \left \| A z - b \right \|_2 \)
     * is minimized. If there are many solutions, then \( z = A^+ b \) is the smallest solution,
     * meaning \( \left \| z \right \|_2 \) is minimized.
     *
     *
     * Note however that some decompositions cannot compute a pseudo-inverse for all matrices.
     * For example, the [LUDecomposition] is not defined for non-square matrices to begin
     * with. The [QRDecomposition] can operate on non-square matrices, but will throw
     * [SingularMatrixException] if the decomposed matrix is singular. Refer to the javadoc
     * of specific decomposition implementations for more details.
     *
     * @return pseudo-inverse matrix (which is the inverse, if it exists),
     * if the decomposition can pseudo-invert the decomposed matrix
     * @throws SingularMatrixException if the decomposed matrix is singular and the decomposition
     * can not compute a pseudo-inverse
     */
    @Throws(SingularMatrixException::class)
    fun getInverse(): RealMatrix?
}