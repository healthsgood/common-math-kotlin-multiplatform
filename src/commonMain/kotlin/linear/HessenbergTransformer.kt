package linear

import linear.exception.NonSquareMatrixException
import util.FastMath.abs
import util.FastMath.sqrt
import util.Precision.equals

/**
 * Class transforming a general real matrix to Hessenberg form.
 *
 * A m  m matrix A can be written as the product of three matrices: A = P
 *  H  P<sup>T</sup> with P an orthogonal matrix and H a Hessenberg
 * matrix. Both P and H are m  m matrices.
 *
 * Transformation to Hessenberg form is often not a goal by itself, but it is an
 * intermediate step in more general decomposition algorithms like
 * [eigen decomposition][EigenDecomposition]. This class is therefore
 * intended for internal use by the library and is not public. As a consequence
 * of this explicitly limited scope, many methods directly returns references to
 * internal arrays, not copies.
 *
 * This class is based on the method orthes in class EigenvalueDecomposition
 * from the [JAMA](http://math.nist.gov/javanumerics/jama/) library.
 *
 * @see [MathWorld](http://mathworld.wolfram.com/HessenbergDecomposition.html)
 *
 * @see [Householder Transformations](http://en.wikipedia.org/wiki/Householder_transformation)
 *
 * @since 3.1
 */
internal class HessenbergTransformer(matrix: RealMatrix) {
    /**
     * Get the Householder vectors of the transform.
     *
     * Note that since this class is only intended for internal use, it returns
     * directly a reference to its internal arrays, not a copy.
     *
     * @return the main diagonal elements of the B matrix
     */
    /**
     * Householder vectors.
     */
    val householderVectorsRef: Array<DoubleArray?>?

    /**
     * Temporary storage vector.
     */
    private val ort: DoubleArray

    /**
     * Cached value of P.
     */
    private var cachedP: RealMatrix?

    /**
     * Cached value of Pt.
     */
    private var cachedPt: RealMatrix?

    /**
     * Cached value of H.
     */
    private var cachedH: RealMatrix?// Double division avoids possible underflow

    /**
     * Returns the matrix P of the transform.
     *
     * P is an orthogonal matrix, i.e. its inverse is also its transpose.
     *
     * @return the P matrix
     */
    val p: RealMatrix?
        get() {
            if (cachedP == null) {
                val n = householderVectorsRef!!.size
                val high = n - 1
                val pa = Array(n) { DoubleArray(n) }
                for (i in 0 until n) {
                    for (j in 0 until n) {
                        pa[i][j] = if (i == j) 1.0 else 0.0
                    }
                }
                for (m in high - 1 downTo 1) {
                    if (householderVectorsRef[m]!![m - 1] != 0.0) {
                        for (i in m + 1..high) {
                            ort[i] = householderVectorsRef[i]!![m - 1]
                        }
                        for (j in m..high) {
                            var g = 0.0
                            for (i in m..high) {
                                g += ort[i] * pa[i][j]
                            }

                            // Double division avoids possible underflow
                            g = g / ort[m] / householderVectorsRef[m]!![m - 1]
                            for (i in m..high) {
                                pa[i][j] += g * ort[i]
                            }
                        }
                    }
                }
                cachedP = MatrixUtils.createRealMatrix(pa as Array<DoubleArray?>?)
            }
            return cachedP
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
                cachedPt = p!!.transpose()
            }

            // return the cached matrix
            return cachedPt
        }// copy the entry of the lower sub-diagonal

    // copy upper triangular part of the matrix

    // return the cached matrix
    /**
     * Returns the Hessenberg matrix H of the transform.
     *
     * @return the H matrix
     */
    val h: RealMatrix?
        get() {
            if (cachedH == null) {
                val m = householderVectorsRef!!.size
                val h = Array(m) { DoubleArray(m) }
                for (i in 0 until m) {
                    if (i > 0) {
                        // copy the entry of the lower sub-diagonal
                        h[i][i - 1] = householderVectorsRef[i]!![i - 1]
                    }

                    // copy upper triangular part of the matrix
                    for (j in i until m) {
                        h[i][j] = householderVectorsRef[i]!![j]
                    }
                }
                cachedH = MatrixUtils.createRealMatrix(h as Array<DoubleArray?>?)
            }

            // return the cached matrix
            return cachedH
        }

    /**
     * Transform original matrix to Hessenberg form.
     *
     * Transformation is done using Householder transforms.
     */
    private fun transform() {
        val n = householderVectorsRef!!.size
        val high = n - 1
        for (m in 1..high - 1) {
            // Scale column.
            var scale = 0.0
            for (i in m..high) {
                scale += abs(householderVectorsRef[i]!![m - 1])
            }
            if (!equals(scale, 0.0)) {
                // Compute Householder transformation.
                var h = 0.0
                for (i in high downTo m) {
                    ort[i] = householderVectorsRef[i]!![m - 1] / scale
                    h += ort[i] * ort[i]
                }
                val g = if (ort[m] > 0) -sqrt(h) else sqrt(h)
                h -= ort[m] * g
                ort[m] -= g

                // Apply Householder similarity transformation
                // H = (I - u*u' / h) * H * (I - u*u' / h)
                for (j in m until n) {
                    var f = 0.0
                    for (i in high downTo m) {
                        f += ort[i] * householderVectorsRef[i]!![j]
                    }
                    f /= h
                    for (i in m..high) {
                        householderVectorsRef[i]!![j] -= f * ort[i]
                    }
                }
                for (i in 0..high) {
                    var f = 0.0
                    for (j in high downTo m) {
                        f += ort[j] * householderVectorsRef[i]!![j]
                    }
                    f /= h
                    for (j in m..high) {
                        householderVectorsRef[i]!![j] -= f * ort[j]
                    }
                }
                ort[m] = scale * ort[m]
                householderVectorsRef[m]!![m - 1] = scale * g
            }
        }
    }

    /**
     * Build the transformation to Hessenberg form of a general matrix.
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
        val m = matrix.getRowDimension()
        householderVectorsRef = matrix.getData()
        ort = DoubleArray(m)
        cachedP = null
        cachedPt = null
        cachedH = null

        // transform matrix
        transform()
    }
}