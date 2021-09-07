package linear

import exception.MaxCountExceededException
import exception.util.LocalizedFormats
import linear.MatrixUtils.createRealMatrix
import linear.exception.NonSquareMatrixException
import util.FastMath.abs
import util.FastMath.max
import util.FastMath.min
import util.FastMath.sqrt
import util.Precision
import util.Precision.equalsByEps

/**
 * Class transforming a general real matrix to Schur form.
 *
 * A m  m matrix A can be written as the product of three matrices: A = P
 *  T  P<sup>T</sup> with P an orthogonal matrix and T an quasi-triangular
 * matrix. Both P and T are m  m matrices.
 *
 * Transformation to Schur form is often not a goal by itself, but it is an
 * intermediate step in more general decomposition algorithms like
 * [eigen decomposition][EigenDecomposition]. This class is therefore
 * intended for internal use by the library and is not public. As a consequence
 * of this explicitly limited scope, many methods directly returns references to
 * internal arrays, not copies.
 *
 * This class is based on the method hqr2 in class EigenvalueDecomposition
 * from the [JAMA](http://math.nist.gov/javanumerics/jama/) library.
 *
 * @see [Schur Decomposition - MathWorld](http://mathworld.wolfram.com/SchurDecomposition.html)
 *
 * @see [Schur Decomposition - Wikipedia](http://en.wikipedia.org/wiki/Schur_decomposition)
 *
 * @see [Householder Transformations](http://en.wikipedia.org/wiki/Householder_transformation)
 *
 * @since 3.1
 */
internal class SchurTransformer(matrix: RealMatrix) {
    /**
     * P matrix.
     */
    private val matrixP: Array<DoubleArray?>?

    /**
     * T matrix.
     */
    private val matrixT: Array<DoubleArray?>?

    /**
     * Epsilon criteria taken from JAMA code (originally was 2^-52).
     */
    private val epsilon = Precision.EPSILON

    /**
     * Cached value of P.
     */
    private var cachedP: RealMatrix?

    /**
     * Cached value of T.
     */
    private var cachedT: RealMatrix?

    /**
     * Cached value of PT.
     */
    private var cachedPt: RealMatrix?

    /**
     * Returns the matrix P of the transform.
     *
     * P is an orthogonal matrix, i.e. its inverse is also its transpose.
     *
     * @return the P matrix
     */
    val p: RealMatrix
        get() {
            if (cachedP == null) {
                cachedP = createRealMatrix(matrixP)
            }
            return cachedP!!
        }// return the cached matrix

    /**
     * Returns the transpose of the matrix P of the transform.
     *
     * P is an orthogonal matrix, i.e. its inverse is also its transpose.
     *
     * @return the transpose of the P matrix
     */
    val pT: RealMatrix?
        get() {
            if (cachedPt == null) {
                cachedPt = p.transpose()
            }

            // return the cached matrix
            return cachedPt
        }// return the cached matrix

    /**
     * Returns the quasi-triangular Schur matrix T of the transform.
     *
     * @return the T matrix
     */
    val t: RealMatrix
        get() {
            if (cachedT == null) {
                cachedT = createRealMatrix(matrixT)
            }

            // return the cached matrix
            return cachedT!!
        }

    /**
     * Transform original matrix to Schur form.
     *
     * @throws MaxCountExceededException if the transformation does not converge
     */
    private fun transform() {
        val n = matrixT!!.size

        // compute matrix norm
        val norm = norm

        // shift information
        val shift = ShiftInfo()

        // Outer loop over eigenvalue index
        var iteration = 0
        var iu = n - 1
        while (iu >= 0) {

            // Look for single small sub-diagonal element
            val il = findSmallSubDiagonalElement(iu, norm)

            // Check for convergence
            if (il == iu) {
                // One root found
                matrixT[iu]!![iu] += shift.exShift
                iu--
                iteration = 0
            } else if (il == iu - 1) {
                // Two roots found
                var p = (matrixT[iu - 1]!![iu - 1] - matrixT[iu]!![iu]) / 2.0
                var q = p * p + matrixT[iu]!![iu - 1] * matrixT[iu - 1]!![iu]
                matrixT[iu]!![iu] += shift.exShift
                matrixT[iu - 1]!![iu - 1] += shift.exShift
                if (q >= 0) {
                    var z = sqrt(abs(q))
                    z = if (p >= 0) {
                        p + z
                    } else {
                        p - z
                    }
                    val x = matrixT[iu]!![iu - 1]
                    val s = abs(x) + abs(z)
                    p = x / s
                    q = z / s
                    val r = sqrt(p * p + q * q)
                    p /= r
                    q /= r

                    // Row modification
                    for (j in iu - 1 until n) {
                        z = matrixT[iu - 1]!![j]
                        matrixT[iu - 1]!![j] = q * z + p * matrixT[iu]!![j]
                        matrixT[iu]!![j] = q * matrixT[iu]!![j] - p * z
                    }

                    // Column modification
                    for (i in 0..iu) {
                        z = matrixT[i]!![iu - 1]
                        matrixT[i]!![iu - 1] = q * z + p * matrixT[i]!![iu]
                        matrixT[i]!![iu] = q * matrixT[i]!![iu] - p * z
                    }

                    // Accumulate transformations
                    for (i in 0..n - 1) {
                        z = matrixP!![i]!![iu - 1]
                        matrixP[i]!![iu - 1] = q * z + p * matrixP[i]!![iu]
                        matrixP[i]!![iu] = q * matrixP[i]!![iu] - p * z
                    }
                }
                iu -= 2
                iteration = 0
            } else {
                // No convergence yet
                computeShift(il, iu, iteration, shift)

                // stop transformation after too many iterations
                if (++iteration > MAX_ITERATIONS) {
                    throw MaxCountExceededException(
                        LocalizedFormats.CONVERGENCE_FAILED,
                        MAX_ITERATIONS
                    )
                }

                // the initial houseHolder vector for the QR step
                val hVec = DoubleArray(3)
                val im = initQRStep(il, iu, shift, hVec)
                performDoubleQRStep(il, im, iu, shift, hVec)
            }
        }
    }// as matrix T is (quasi-)triangular, also take the sub-diagonal element into account

    /**
     * Computes the L1 norm of the (quasi-)triangular matrix T.
     *
     * @return the L1 norm of matrix T
     */
    private val norm: Double
        private get() {
            var norm = 0.0
            for (i in matrixT!!.indices) {
                // as matrix T is (quasi-)triangular, also take the sub-diagonal element into account
                for (j in max(i - 1, 0) until matrixT.size) {
                    norm += abs(matrixT[i]!![j])
                }
            }
            return norm
        }

    /**
     * Find the first small sub-diagonal element and returns its index.
     *
     * @param startIdx the starting index for the search
     * @param norm     the L1 norm of the matrix
     * @return the index of the first small sub-diagonal element
     */
    private fun findSmallSubDiagonalElement(startIdx: Int, norm: Double): Int {
        var l = startIdx
        while (l > 0) {
            var s = abs(matrixT!![l - 1]!![l - 1]) + abs(matrixT[l]!![l])
            if (s == 0.0) {
                s = norm
            }
            if (abs(matrixT[l]!![l - 1]) < epsilon * s) {
                break
            }
            l--
        }
        return l
    }

    /**
     * Compute the shift for the current iteration.
     *
     * @param l         the index of the small sub-diagonal element
     * @param idx       the current eigenvalue index
     * @param iteration the current iteration
     * @param shift     holder for shift information
     */
    private fun computeShift(l: Int, idx: Int, iteration: Int, shift: ShiftInfo) {
        // Form shift
        shift.x = matrixT!![idx]!![idx]
        shift.w = 0.0
        shift.y = shift.w
        if (l < idx) {
            shift.y = matrixT[idx - 1]!![idx - 1]
            shift.w = matrixT[idx]!![idx - 1] * matrixT[idx - 1]!![idx]
        }

        // Wilkinson's original ad hoc shift
        if (iteration == 10) {
            shift.exShift += shift.x
            for (i in 0..idx) {
                matrixT[i]!![i] -= shift.x
            }
            val s = abs(
                matrixT[idx]!![idx - 1]
            ) + abs(matrixT[idx - 1]!![idx - 2])
            shift.x = 0.75 * s
            shift.y = 0.75 * s
            shift.w = -0.4375 * s * s
        }

        // MATLAB's new ad hoc shift
        if (iteration == 30) {
            var s = (shift.y - shift.x) / 2.0
            s = s * s + shift.w
            if (s > 0.0) {
                s = sqrt(s)
                if (shift.y < shift.x) {
                    s = -s
                }
                s = shift.x - shift.w / ((shift.y - shift.x) / 2.0 + s)
                for (i in 0..idx) {
                    matrixT[i]!![i] -= s
                }
                shift.exShift += s
                shift.w = 0.964
                shift.y = shift.w
                shift.x = shift.y
            }
        }
    }

    /**
     * Initialize the householder vectors for the QR step.
     *
     * @param il    the index of the small sub-diagonal element
     * @param iu    the current eigenvalue index
     * @param shift shift information holder
     * @param hVec  the initial houseHolder vector
     * @return the start index for the QR step
     */
    private fun initQRStep(il: Int, iu: Int, shift: ShiftInfo, hVec: DoubleArray): Int {
        // Look for two consecutive small sub-diagonal elements
        var im = iu - 2
        while (im >= il) {
            val z = matrixT!![im]!![im]
            val r = shift.x - z
            val s = shift.y - z
            hVec[0] = (r * s - shift.w) / matrixT[im + 1]!![im] + matrixT[im]!![im + 1]
            hVec[1] = matrixT[im + 1]!![im + 1] - z - r - s
            hVec[2] = matrixT[im + 2]!![im + 1]
            if (im == il) {
                break
            }
            val lhs = abs(matrixT[im]!![im - 1]) * (abs(hVec[1]) + abs(hVec[2]))
            val rhs = abs(hVec[0]) * (abs(matrixT[im - 1]!![im - 1]) +
                    abs(z) +
                    abs(matrixT[im + 1]!![im + 1]))
            if (lhs < epsilon * rhs) {
                break
            }
            im--
        }
        return im
    }

    /**
     * Perform a double QR step involving rows l:idx and columns m:n
     *
     * @param il    the index of the small sub-diagonal element
     * @param im    the start index for the QR step
     * @param iu    the current eigenvalue index
     * @param shift shift information holder
     * @param hVec  the initial houseHolder vector
     */
    private fun performDoubleQRStep(
        il: Int, im: Int, iu: Int,
        shift: ShiftInfo, hVec: DoubleArray
    ) {
        val n = matrixT!!.size
        var p = hVec[0]
        var q = hVec[1]
        var r = hVec[2]
        for (k in im..iu - 1) {
            val notlast = k != iu - 1
            if (k != im) {
                p = matrixT[k]!![k - 1]
                q = matrixT[k + 1]!![k - 1]
                r = if (notlast) matrixT[k + 2]!![k - 1] else 0.0
                shift.x = abs(p) + abs(q) + abs(r)
                if (equalsByEps(shift.x, 0.0, epsilon)) {
                    continue
                }
                p /= shift.x
                q /= shift.x
                r /= shift.x
            }
            var s = sqrt(p * p + q * q + r * r)
            if (p < 0.0) {
                s = -s
            }
            if (s != 0.0) {
                if (k != im) {
                    matrixT[k]!![k - 1] = -s * shift.x
                } else if (il != im) {
                    matrixT[k]!![k - 1] = -matrixT[k]!![k - 1]
                }
                p += s
                shift.x = p / s
                shift.y = q / s
                val z = r / s
                q /= p
                r /= p

                // Row modification
                for (j in k until n) {
                    p = matrixT[k]!![j] + q * matrixT[k + 1]!![j]
                    if (notlast) {
                        p += r * matrixT[k + 2]!![j]
                        matrixT[k + 2]!![j] -= p * z
                    }
                    matrixT[k]!![j] -= p * shift.x
                    matrixT[k + 1]!![j] -= p * shift.y
                }

                // Column modification
                for (i in 0..min(iu, k + 3)) {
                    p = shift.x * matrixT[i]!![k] + shift.y * matrixT[i]!![k + 1]
                    if (notlast) {
                        p += z * matrixT[i]!![k + 2]
                        matrixT[i]!![k + 2] -= p * r
                    }
                    matrixT[i]!![k] -= p
                    matrixT[i]!![k + 1] -= p * q
                }

                // Accumulate transformations
                val high = matrixT.size - 1
                for (i in 0..high) {
                    p = shift.x * matrixP!![i]!![k] + shift.y * matrixP[i]!![k + 1]
                    if (notlast) {
                        p += z * matrixP[i]!![k + 2]
                        matrixP[i]!![k + 2] -= p * r
                    }
                    matrixP[i]!![k] -= p
                    matrixP[i]!![k + 1] -= p * q
                }
            } // (s != 0)
        } // k loop

        // clean up pollution due to round-off errors
        for (i in im + 2..iu) {
            matrixT[i]!![i - 2] = 0.0
            if (i > im + 2) {
                matrixT[i]!![i - 3] = 0.0
            }
        }
    }

    /**
     * Internal data structure holding the current shift information.
     * Contains variable names as present in the original JAMA code.
     */
    private class ShiftInfo {
        // CHECKSTYLE: stop all
        /**
         * x shift info
         */
        var x = 0.0

        /**
         * y shift info
         */
        var y = 0.0

        /**
         * w shift info
         */
        var w = 0.0

        /**
         * Indicates an exceptional shift.
         */
        var exShift = 0.0 // CHECKSTYLE: resume all
    }

    companion object {
        /**
         * Maximum allowed iterations for convergence of the transformation.
         */
        private const val MAX_ITERATIONS = 100
    }

    /**
     * Build the transformation to Schur form of a general real matrix.
     *
     * @param matrix matrix to transform
     * @throws NonSquareMatrixException if the matrix is not square
     */
    init {
        if (!matrix.isSquare) {
            throw NonSquareMatrixException(
                matrix.getRowDimension(),
                matrix.getColumnDimension()
            )
        }
        val transformer = HessenbergTransformer(matrix)
        matrixT = transformer.h!!.getData()
        matrixP = transformer.p!!.getData()
        cachedT = null
        cachedP = null
        cachedPt = null

        // transform matrix
        transform()
    }
}