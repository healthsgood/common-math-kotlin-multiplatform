package linear

import crossjvm.ArrayUtil
import linear.MatrixUtils.createRealMatrix
import linear.exception.NonSquareMatrixException
import util.FastMath.sqrt

/**
 * Class transforming a symmetrical matrix to tridiagonal shape.
 *
 * A symmetrical m  m matrix A can be written as the product of three matrices:
 * A = Q  T  Q<sup>T</sup> with Q an orthogonal matrix and T a symmetrical
 * tridiagonal matrix. Both Q and T are m  m matrices.
 *
 * This implementation only uses the upper part of the matrix, the part below the
 * diagonal is not accessed at all.
 *
 * Transformation to tridiagonal shape is often not a goal by itself, but it is
 * an intermediate step in more general decomposition algorithms like [ ]. This class is therefore intended for internal
 * use by the library and is not public. As a consequence of this explicitly limited scope,
 * many methods directly returns references to internal arrays, not copies.
 *
 * @since 2.0
 */
internal class TriDiagonalTransformer(matrix: RealMatrix) {
    /**
     * Get the Householder vectors of the transform.
     *
     * Note that since this class is only intended for internal use,
     * it returns directly a reference to its internal arrays, not a copy.
     *
     * @return the main diagonal elements of the B matrix
     */
    /**
     * Householder vectors.
     */
    val householderVectorsRef: Array<DoubleArray?>?
    /**
     * Get the main diagonal elements of the matrix T of the transform.
     *
     * Note that since this class is only intended for internal use,
     * it returns directly a reference to its internal arrays, not a copy.
     *
     * @return the main diagonal elements of the T matrix
     */
    /**
     * Main diagonal.
     */
    val mainDiagonalRef: DoubleArray
    /**
     * Get the secondary diagonal elements of the matrix T of the transform.
     *
     * Note that since this class is only intended for internal use,
     * it returns directly a reference to its internal arrays, not a copy.
     *
     * @return the secondary diagonal elements of the T matrix
     */
    /**
     * Secondary diagonal.
     */
    val secondaryDiagonalRef: DoubleArray

    /**
     * Cached value of Q.
     */
    private var cachedQ: RealMatrix?

    /**
     * Cached value of Qt.
     */
    private var cachedQt: RealMatrix?

    /**
     * Cached value of T.
     */
    private var cachedT: RealMatrix?

    /**
     * Returns the matrix Q of the transform.
     *
     * Q is an orthogonal matrix, i.e. its transpose is also its inverse.
     *
     * @return the Q matrix
     */
    val q: RealMatrix?
        get() {
            if (cachedQ == null) {
                cachedQ = qT.transpose()
            }
            return cachedQ
        }// build up first part of the matrix by applying Householder transforms

    // return the cached matrix
    /**
     * Returns the transpose of the matrix Q of the transform.
     *
     * Q is an orthogonal matrix, i.e. its transpose is also its inverse.
     *
     * @return the Q matrix
     */
    val qT: RealMatrix
        get() {
            if (cachedQt == null) {
                val m = householderVectorsRef!!.size
                val qta = Array<DoubleArray?>(m) { DoubleArray(m) }

                // build up first part of the matrix by applying Householder transforms
                for (k in m - 1 downTo 1) {
                    val hK = householderVectorsRef[k - 1]
                    qta[k]!![k] = 1.0
                    if (hK!![k] != 0.0) {
                        val inv = 1.0 / (secondaryDiagonalRef[k - 1] * hK[k])
                        var beta = 1.0 / secondaryDiagonalRef[k - 1]
                        qta[k]!![k] = 1 + beta * hK[k]
                        for (i in k + 1 until m) {
                            qta[k]!![i] = beta * hK[i]
                        }
                        for (j in k + 1 until m) {
                            beta = 0.0
                            for (i in k + 1 until m) {
                                beta += qta[j]!![i] * hK[i]
                            }
                            beta *= inv
                            qta[j]!![k] = beta * hK[k]
                            for (i in k + 1 until m) {
                                qta[j]!![i] += beta * hK[i]
                            }
                        }
                    }
                }
                qta[0]!![0] = 1.0
                cachedQt = createRealMatrix(qta)
            }

            // return the cached matrix
            return cachedQt!!
        }// return the cached matrix

    /**
     * Returns the tridiagonal matrix T of the transform.
     *
     * @return the T matrix
     */
    val t: RealMatrix
        get() {
            if (cachedT == null) {
                val m = mainDiagonalRef.size
                val ta = Array<DoubleArray?>(m) { DoubleArray(m) }
                for (i in 0 until m) {
                    ta[i]!![i] = mainDiagonalRef[i]
                    if (i > 0) {
                        ta[i]!![i - 1] = secondaryDiagonalRef[i - 1]
                    }
                    if (i < mainDiagonalRef.size - 1) {
                        ta[i]!![i + 1] = secondaryDiagonalRef[i]
                    }
                }
                cachedT = createRealMatrix(ta)
            }

            // return the cached matrix
            return cachedT!!
        }

    /**
     * Transform original matrix to tridiagonal form.
     *
     * Transformation is done using Householder transforms.
     */
    private fun transform() {
        val m = householderVectorsRef!!.size
        val z = DoubleArray(m)
        for (k in 0 until m - 1) {

            //zero-out a row and a column simultaneously
            val hK = householderVectorsRef[k]
            mainDiagonalRef[k] = hK!![k]
            var xNormSqr = 0.0
            for (j in k + 1 until m) {
                val c = hK[j]
                xNormSqr += c * c
            }
            val a = if (hK[k + 1] > 0) -sqrt(xNormSqr) else sqrt(xNormSqr)
            secondaryDiagonalRef[k] = a
            if (a != 0.0) {
                // apply Householder transform from left and right simultaneously
                hK[k + 1] -= a
                val beta = -1 / (a * hK[k + 1])

                // compute a = beta A v, where v is the Householder vector
                // this loop is written in such a way
                //   1) only the upper triangular part of the matrix is accessed
                //   2) access is cache-friendly for a matrix stored in rows
                ArrayUtil.fill(z, k + 1, m, 0.0)
                for (i in k + 1 until m) {
                    val hI = householderVectorsRef[i]
                    val hKI = hK[i]
                    var zI = hI!![i] * hKI
                    for (j in i + 1 until m) {
                        val hIJ = hI[j]
                        zI += hIJ * hK[j]
                        z[j] += hIJ * hKI
                    }
                    z[i] = beta * (z[i] + zI)
                }

                // compute gamma = beta vT z / 2
                var gamma = 0.0
                for (i in k + 1 until m) {
                    gamma += z[i] * hK[i]
                }
                gamma *= beta / 2

                // compute z = z - gamma v
                for (i in k + 1 until m) {
                    z[i] -= gamma * hK[i]
                }

                // update matrix: A = A - v zT - z vT
                // only the upper triangular part of the matrix is updated
                for (i in k + 1 until m) {
                    val hI = householderVectorsRef[i]
                    for (j in i until m) {
                        hI!![j] -= hK[i] * z[j] + z[i] * hK[j]
                    }
                }
            }
        }
        mainDiagonalRef[m - 1] = householderVectorsRef[m - 1]!![m - 1]
    }

    /**
     * Build the transformation to tridiagonal shape of a symmetrical matrix.
     *
     * The specified matrix is assumed to be symmetrical without any check.
     * Only the upper triangular part of the matrix is used.
     *
     * @param matrix Symmetrical matrix to transform.
     * @throws NonSquareMatrixException if the matrix is not square.
     */
    init {
        if (!matrix.isSquare) {
            throw NonSquareMatrixException(
                matrix.getRowDimension(),
                matrix.getColumnDimension()
            )
        }
        val m = matrix.getRowDimension()
        householderVectorsRef = matrix.getData()
        mainDiagonalRef = DoubleArray(m)
        secondaryDiagonalRef = DoubleArray(m - 1)
        cachedQ = null
        cachedQt = null
        cachedT = null

        // transform matrix
        transform()
    }
}