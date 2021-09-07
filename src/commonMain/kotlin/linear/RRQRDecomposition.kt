package linear

import linear.MatrixUtils.createRealIdentityMatrix
import linear.MatrixUtils.createRealMatrix
import util.FastMath.min

/**
 * Calculates the rank-revealing QR-decomposition of a matrix, with column pivoting.
 *
 * The rank-revealing QR-decomposition of a matrix A consists of three matrices Q,
 * R and P such that AP=QR.  Q is orthogonal (Q<sup>T</sup>Q = I), and R is upper triangular.
 * If A is mn, Q is mm and R is mn and P is nn.
 *
 * QR decomposition with column pivoting produces a rank-revealing QR
 * decomposition and the [.getRank] method may be used to return the rank of the
 * input matrix A.
 *
 * This class compute the decomposition using Householder reflectors.
 *
 * For efficiency purposes, the decomposition in packed form is transposed.
 * This allows inner loop to iterate inside rows, which is much more cache-efficient
 * in Java.
 *
 * This class is based on the class with similar name from the
 * [JAMA](http://math.nist.gov/javanumerics/jama/) library, with the
 * following changes:
 *
 *  * a [getQT][.getQT] method has been added,
 *  * the `solve` and `isFullRank` methods have been replaced
 * by a [getSolver][.getSolver] method and the equivalent methods
 * provided by the returned [DecompositionSolver].
 *
 *
 * @see [MathWorld](http://mathworld.wolfram.com/QRDecomposition.html)
 *
 * @see [Wikipedia](http://en.wikipedia.org/wiki/QR_decomposition)
 *
 * @since 3.2
 */
class RRQRDecomposition
/**
 * Calculates the QR-decomposition of the given matrix.
 *
 * @param matrix    The matrix to decompose.
 * @param threshold Singularity threshold.
 * @see .RRQRDecomposition
 */
/**
 * Calculates the QR-decomposition of the given matrix.
 * The singularity threshold defaults to zero.
 *
 * @param matrix The matrix to decompose.
 * @see .RRQRDecomposition
 */
constructor(matrix: RealMatrix?, threshold: Double = 0.0) : QRDecomposition(matrix!!, threshold) {
    /**
     * An array to record the column pivoting for later creation of P.
     */
    private lateinit var p: IntArray

    /**
     * Cached value of P.
     */
    private var cachedP: RealMatrix? = null

    /**
     * Decompose matrix.
     *
     * @param qrt transposed matrix
     */
    override fun decompose(qrt: Array<DoubleArray?>?) {
        p = IntArray(qrt!!.size)
        for (i in p.indices) {
            p[i] = i
        }
        super.decompose(qrt)
    }

    /**
     * Perform Householder reflection for a minor A(minor, minor) of A.
     *
     * @param minor minor index
     * @param qrt   transposed matrix
     */
    override fun performHouseholderReflection(minor: Int, qrt: Array<DoubleArray?>?) {
        var l2NormSquaredMax = 0.0
        // Find the unreduced column with the greatest L2-Norm
        var l2NormSquaredMaxIndex = minor
        for (i in minor until qrt!!.size) {
            var l2NormSquared = 0.0
            for (j in 0 until qrt[i]!!.size) {
                l2NormSquared += qrt[i]!![j] * qrt[i]!![j]
            }
            if (l2NormSquared > l2NormSquaredMax) {
                l2NormSquaredMax = l2NormSquared
                l2NormSquaredMaxIndex = i
            }
        }
        // swap the current column with that with the greated L2-Norm and record in p
        if (l2NormSquaredMaxIndex != minor) {
            val tmp1 = qrt[minor]
            qrt[minor] = qrt[l2NormSquaredMaxIndex]
            qrt[l2NormSquaredMaxIndex] = tmp1
            val tmp2 = p[minor]
            p[minor] = p[l2NormSquaredMaxIndex]
            p[l2NormSquaredMaxIndex] = tmp2
        }
        super.performHouseholderReflection(minor, qrt)
    }

    /**
     * Returns the pivot matrix, P, used in the QR Decomposition of matrix A such that AP = QR.
     *
     *
     * If no pivoting is used in this decomposition then P is equal to the identity matrix.
     *
     * @return a permutation matrix.
     */
    fun getP(): RealMatrix {
        if (cachedP == null) {
            val n = p.size
            cachedP = createRealMatrix(n, n)
            for (i in 0 until n) {
                cachedP!!.setEntry(p[i], i, 1.0)
            }
        }
        return cachedP!!
    }

    /**
     * Return the effective numerical matrix rank.
     *
     * The effective numerical rank is the number of non-negligible
     * singular values.
     *
     * This implementation looks at Frobenius norms of the sequence of
     * bottom right submatrices.  When a large fall in norm is seen,
     * the rank is returned. The drop is computed as:
     * <pre>
     * (thisNorm/lastNorm) * rNorm < dropThreshold
    </pre> *
     *
     *
     * where thisNorm is the Frobenius norm of the current submatrix,
     * lastNorm is the Frobenius norm of the previous submatrix,
     * rNorm is is the Frobenius norm of the complete matrix
     *
     *
     * @param dropThreshold threshold triggering rank computation
     * @return effective numerical matrix rank
     */
    fun getRank(dropThreshold: Double): Int {
        val r = r
        val rows = r.getRowDimension()
        val columns = r.getColumnDimension()
        var rank = 1
        var lastNorm = r.getFrobeniusNorm()
        val rNorm = lastNorm
        while (rank < min(rows, columns)) {
            val thisNorm = r.getSubMatrix(rank, rows - 1, rank, columns - 1)!!.getFrobeniusNorm()
            if (thisNorm == 0.0 || thisNorm / lastNorm * rNorm < dropThreshold) {
                break
            }
            lastNorm = thisNorm
            rank++
        }
        return rank
    }

    /**
     * Get a solver for finding the A  X = B solution in least square sense.
     *
     *
     * Least Square sense means a solver can be computed for an overdetermined system,
     * (i.e. a system with more equations than unknowns, which corresponds to a tall A
     * matrix with more rows than columns). In any case, if the matrix is singular
     * within the tolerance set at [construction][RRQRDecomposition.RRQRDecomposition], an error will be triggered when
     * the [solve][DecompositionSolver.solve] method will be called.
     *
     *
     * @return a solver
     */
    override val solver: DecompositionSolver
        get() = Solver(super.solver, getP())

    /**
     * Specialized solver.
     */
    private class Solver
    /**
     * Build a solver from decomposed matrix.
     *
     * @param upper upper level solver.
     * @param p     permutation matrix
     */(
        /**
         * Upper level solver.
         */
        private val upper: DecompositionSolver?,
        /**
         * A permutation matrix for the pivots used in the QR decomposition
         */
        private val p: RealMatrix
    ) : DecompositionSolver {
        /**
         * {@inheritDoc}
         */
        override val isNonSingular: Boolean
            get() = upper!!.isNonSingular

        /**
         * {@inheritDoc}
         */
        override fun solve(b: RealVector?): RealVector? {
            return p.operate(upper!!.solve(b))
        }

        /**
         * {@inheritDoc}
         */
        override fun solve(b: RealMatrix?): RealMatrix? {
            return p.multiply(upper!!.solve(b))
        }

        /**
         * {@inheritDoc}
         *
         * @throws SingularMatrixException if the decomposed matrix is singular.
         */
        override fun getInverse(): RealMatrix? {
            return solve(createRealIdentityMatrix(p.getRowDimension()))
        }
    }
}