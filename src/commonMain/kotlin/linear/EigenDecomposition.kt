package linear

import complex.Complex
import exception.DimensionMismatchException
import exception.MathArithmeticException
import exception.MathUnsupportedOperationException
import exception.MaxCountExceededException
import exception.util.LocalizedFormats
import linear.exception.SingularMatrixException
import util.FastMath.abs
import util.FastMath.max
import util.FastMath.min
import util.FastMath.sqrt
import util.Precision
import util.Precision.compareTo
import util.Precision.equals
import util.Precision.equalsByEps

/**
 * Calculates the eigen decomposition of a real matrix.
 *
 * The eigen decomposition of matrix A is a set of two matrices:
 * V and D such that A = V  D  V<sup>T</sup>.
 * A, V and D are all m  m matrices.
 *
 * This class is similar in spirit to the `EigenvalueDecomposition`
 * class from the [JAMA](http://math.nist.gov/javanumerics/jama/)
 * library, with the following changes:
 *
 *  * a [getVt][.getVT] method has been added,
 *  * two [getRealEigenvalue][.getRealEigenvalue] and [   getImagEigenvalue][.getImagEigenvalue] methods to pick up a single eigenvalue have been added,
 *  * a [getEigenvector][.getEigenvector] method to pick up a single
 * eigenvector has been added,
 *  * a [getDeterminant][.getDeterminant] method has been added.
 *  * a [getSolver][.getSolver] method has been added.
 *
 *
 *
 * As of 3.1, this class supports general real matrices (both symmetric and non-symmetric):
 *
 *
 *
 * If A is symmetric, then A = V*D*V' where the eigenvalue matrix D is diagonal and the eigenvector
 * matrix V is orthogonal, i.e. A = V.multiply(D.multiply(V.transpose())) and
 * V.multiply(V.transpose()) equals the identity matrix.
 *
 *
 *
 * If A is not symmetric, then the eigenvalue matrix D is block diagonal with the real eigenvalues
 * in 1-by-1 blocks and any complex eigenvalues, lambda + i*mu, in 2-by-2 blocks:
 * <pre>
 * [lambda, mu    ]
 * [   -mu, lambda]
</pre> *
 * The columns of V represent the eigenvectors in the sense that A*V = V*D,
 * i.e. A.multiply(V) equals V.multiply(D).
 * The matrix V may be badly conditioned, or even singular, so the validity of the equation
 * A = V*D*inverse(V) depends upon the condition of V.
 *
 *
 *
 * This implementation is based on the paper by A. Drubrulle, R.S. Martin and
 * J.H. Wilkinson "The Implicit QL Algorithm" in Wilksinson and Reinsch (1971)
 * Handbook for automatic computation, vol. 2, Linear algebra, Springer-Verlag,
 * New-York
 *
 *
 * @see [MathWorld](http://mathworld.wolfram.com/EigenDecomposition.html)
 *
 * @see [Wikipedia](http://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix)
 *
 * @since 2.0 (changed to concrete class in 3.0)
 */
class EigenDecomposition {
    /**
     * Whether the matrix is symmetric.
     */
    private val isSymmetric: Boolean

    /**
     * Maximum number of iterations accepted in the implicit QL transformation
     */
    private val maxIter: Byte = 30

    /**
     * Main diagonal of the tridiagonal matrix.
     */
    private lateinit var main: DoubleArray

    /**
     * Secondary diagonal of the tridiagonal matrix.
     */
    private lateinit var secondary: DoubleArray

    /**
     * Transformer to tridiagonal (may be null if matrix is already
     * tridiagonal).
     */
    private var transformer: TriDiagonalTransformer? = null

    /**
     * Real part of the realEigenvalues.
     */
    private lateinit var realEigenvalues: DoubleArray

    /**
     * Imaginary part of the realEigenvalues.
     */
    private lateinit var imagEigenvalues: DoubleArray

    /**
     * Eigenvectors.
     */
    private lateinit var eigenvectors: Array<ArrayRealVector?>

    /**
     * Cached value of V.
     */
    private var cachedV: RealMatrix? = null

    /**
     * Cached value of D.
     */
    private var cachedD: RealMatrix? = null

    /**
     * Cached value of Vt.
     */
    private var cachedVt: RealMatrix? = null

    /**
     * Calculates the eigen decomposition of the given real matrix.
     *
     *
     * Supports decomposition of a general matrix since 3.1.
     *
     * @param matrix Matrix to decompose.
     * @throws MaxCountExceededException if the algorithm fails to converge.
     * @throws MathArithmeticException   if the decomposition of a general matrix
     * results in a matrix with zero norm
     * @since 3.1
     */
    constructor(matrix: RealMatrix) {
        val symTol = 10 * matrix.getRowDimension() * matrix.getColumnDimension() * Precision.EPSILON
        isSymmetric = MatrixUtils.isSymmetric(matrix, symTol)
        if (isSymmetric) {
            transformToTridiagonal(matrix)
            findEigenVectors(transformer!!.q!!.getData())
        } else {
            val t = transformToSchur(matrix)
            findEigenVectorsFromSchur(t)
        }
    }

    /**
     * Calculates the eigen decomposition of the given real matrix.
     *
     * @param matrix         Matrix to decompose.
     * @param splitTolerance Dummy parameter (present for backward
     * compatibility only).
     * @throws MathArithmeticException   if the decomposition of a general matrix
     * results in a matrix with zero norm
     * @throws MaxCountExceededException if the algorithm fails to converge.
     */
    @Deprecated("in 3.1 (to be removed in 4.0) due to unused parameter")
    constructor(
        matrix: RealMatrix,
        splitTolerance: Double
    ) : this(matrix)

    /**
     * Calculates the eigen decomposition of the symmetric tridiagonal
     * matrix.  The Householder matrix is assumed to be the identity matrix.
     *
     * @param main      Main diagonal of the symmetric tridiagonal form.
     * @param secondary Secondary of the tridiagonal form.
     * @throws MaxCountExceededException if the algorithm fails to converge.
     * @since 3.1
     */
    constructor(main: DoubleArray, secondary: DoubleArray) {
        isSymmetric = true
        this.main = main.copyOf()
        this.secondary = secondary.copyOf()
        transformer = null
        val size = main.size
        val z = Array<DoubleArray?>(size) { DoubleArray(size) }
        for (i in 0 until size) {
            z[i]!![i] = 1.0
        }
        findEigenVectors(z)
    }

    /**
     * Calculates the eigen decomposition of the symmetric tridiagonal
     * matrix.  The Householder matrix is assumed to be the identity matrix.
     *
     * @param main           Main diagonal of the symmetric tridiagonal form.
     * @param secondary      Secondary of the tridiagonal form.
     * @param splitTolerance Dummy parameter (present for backward
     * compatibility only).
     * @throws MaxCountExceededException if the algorithm fails to converge.
     */
    @Deprecated("in 3.1 (to be removed in 4.0) due to unused parameter")
    constructor(
        main: DoubleArray, secondary: DoubleArray,
        splitTolerance: Double
    ) : this(main, secondary) // return the cached matrix

    /**
     * Gets the matrix V of the decomposition.
     * V is an orthogonal matrix, i.e. its transpose is also its inverse.
     * The columns of V are the eigenvectors of the original matrix.
     * No assumption is made about the orientation of the system axes formed
     * by the columns of V (e.g. in a 3-dimension space, V can form a left-
     * or right-handed system).
     *
     * @return the V matrix.
     */
    val v: RealMatrix?
        get() {
            if (cachedV == null) {
                val m = eigenvectors.size
                cachedV = MatrixUtils.createRealMatrix(m, m)
                for (k in 0 until m) {
                    cachedV!!.setColumnVector(k, eigenvectors[k])
                }
            }
            // return the cached matrix
            return cachedV
        }// cache the matrix for subsequent calls

    /**
     * Gets the block diagonal matrix D of the decomposition.
     * D is a block diagonal matrix.
     * Real eigenvalues are on the diagonal while complex values are on
     * 2x2 blocks { {real +imaginary}, {-imaginary, real} }.
     *
     * @return the D matrix.
     * @see .getRealEigenvalues
     * @see .getImagEigenvalues
     */
    val d: RealMatrix?
        get() {
            if (cachedD == null) {
                // cache the matrix for subsequent calls
                cachedD = MatrixUtils.createRealDiagonalMatrix(realEigenvalues)
                for (i in imagEigenvalues.indices) {
                    if (compareTo(imagEigenvalues[i], 0.0, EPSILON) > 0) {
                        cachedD!!.setEntry(i, i + 1, imagEigenvalues[i])
                    } else if (compareTo(imagEigenvalues[i], 0.0, EPSILON) < 0) {
                        cachedD!!.setEntry(i, i - 1, imagEigenvalues[i])
                    }
                }
            }
            return cachedD
        }// return the cached matrix

    /**
     * Gets the transpose of the matrix V of the decomposition.
     * V is an orthogonal matrix, i.e. its transpose is also its inverse.
     * The columns of V are the eigenvectors of the original matrix.
     * No assumption is made about the orientation of the system axes formed
     * by the columns of V (e.g. in a 3-dimension space, V can form a left-
     * or right-handed system).
     *
     * @return the transpose of the V matrix.
     */
    val vT: RealMatrix?
        get() {
            if (cachedVt == null) {
                val m = eigenvectors.size
                cachedVt = MatrixUtils.createRealMatrix(m, m)
                for (k in 0 until m) {
                    cachedVt!!.setRowVector(k, eigenvectors[k])
                }
            }

            // return the cached matrix
            return cachedVt
        }

    /**
     * Returns whether the calculated eigen values are complex or real.
     *
     * The method performs a zero check for each element of the
     * [.getImagEigenvalues] array and returns `true` if any
     * element is not equal to zero.
     *
     * @return `true` if the eigen values are complex, `false` otherwise
     * @since 3.1
     */
    fun hasComplexEigenvalues(): Boolean {
        for (i in imagEigenvalues.indices) {
            if (!equalsByEps(imagEigenvalues[i], 0.0, EPSILON)) {
                return true
            }
        }
        return false
    }

    /**
     * Gets a copy of the real parts of the eigenvalues of the original matrix.
     *
     * @return a copy of the real parts of the eigenvalues of the original matrix.
     * @see .getD
     * @see .getRealEigenvalue
     * @see .getImagEigenvalues
     */
    fun getRealEigenvalues(): DoubleArray {
        return realEigenvalues.copyOf()
    }

    /**
     * Returns the real part of the i<sup>th</sup> eigenvalue of the original
     * matrix.
     *
     * @param i index of the eigenvalue (counting from 0)
     * @return real part of the i<sup>th</sup> eigenvalue of the original
     * matrix.
     * @see .getD
     * @see .getRealEigenvalues
     * @see .getImagEigenvalue
     */
    fun getRealEigenvalue(i: Int): Double {
        return realEigenvalues[i]
    }

    /**
     * Gets a copy of the imaginary parts of the eigenvalues of the original
     * matrix.
     *
     * @return a copy of the imaginary parts of the eigenvalues of the original
     * matrix.
     * @see .getD
     * @see .getImagEigenvalue
     * @see .getRealEigenvalues
     */
    fun getImagEigenvalues(): DoubleArray {
        return imagEigenvalues.copyOf()
    }

    /**
     * Gets the imaginary part of the i<sup>th</sup> eigenvalue of the original
     * matrix.
     *
     * @param i Index of the eigenvalue (counting from 0).
     * @return the imaginary part of the i<sup>th</sup> eigenvalue of the original
     * matrix.
     * @see .getD
     * @see .getImagEigenvalues
     * @see .getRealEigenvalue
     */
    fun getImagEigenvalue(i: Int): Double {
        return imagEigenvalues[i]
    }

    /**
     * Gets a copy of the i<sup>th</sup> eigenvector of the original matrix.
     *
     * @param i Index of the eigenvector (counting from 0).
     * @return a copy of the i<sup>th</sup> eigenvector of the original matrix.
     * @see .getD
     */
    fun getEigenvector(i: Int): RealVector {
        return eigenvectors[i]!!.copy()
    }

    /**
     * Computes the determinant of the matrix.
     *
     * @return the determinant of the matrix.
     */
    val determinant: Double
        get() {
            var determinant = 1.0
            for (lambda in realEigenvalues) {
                determinant *= lambda
            }
            return determinant
        }

    /**
     * Computes the square-root of the matrix.
     * This implementation assumes that the matrix is symmetric and positive
     * definite.
     *
     * @return the square-root of the matrix.
     * @throws MathUnsupportedOperationException if the matrix is not
     * symmetric or not positive definite.
     * @since 3.1
     */
    val squareRoot: RealMatrix?
        get() {
            if (!isSymmetric) {
                throw MathUnsupportedOperationException()
            }
            val sqrtEigenValues = DoubleArray(realEigenvalues.size)
            for (i in realEigenvalues.indices) {
                val eigen = realEigenvalues[i]
                if (eigen <= 0) {
                    throw MathUnsupportedOperationException()
                }
                sqrtEigenValues[i] = sqrt(eigen)
            }
            val sqrtEigen = MatrixUtils.createRealDiagonalMatrix(sqrtEigenValues)
            val v = v
            val vT = vT
            return v!!.multiply(sqrtEigen)!!.multiply(vT)
        }

    /**
     * Gets a solver for finding the A  X = B solution in exact
     * linear sense.
     *
     *
     * Since 3.1, eigen decomposition of a general matrix is supported,
     * but the [DecompositionSolver] only supports real eigenvalues.
     *
     * @return a solver
     * @throws MathUnsupportedOperationException if the decomposition resulted in
     * complex eigenvalues
     */
    val solver: DecompositionSolver
        get() {
            if (hasComplexEigenvalues()) {
                throw MathUnsupportedOperationException()
            }
            return Solver(realEigenvalues, imagEigenvalues, eigenvectors)
        }

    /**
     * Transforms the matrix to tridiagonal form.
     *
     * @param matrix Matrix to transform.
     */
    private fun transformToTridiagonal(matrix: RealMatrix) {
        // transform the matrix to tridiagonal
        transformer = TriDiagonalTransformer(matrix)
        main = transformer!!.mainDiagonalRef
        secondary = transformer!!.secondaryDiagonalRef
    }

    /**
     * Find eigenvalues and eigenvectors (Dubrulle et al., 1971)
     *
     * @param householderMatrix Householder matrix of the transformation
     * to tridiagonal form.
     */
    private fun findEigenVectors(householderMatrix: Array<DoubleArray?>?) {
        val z = householderMatrix!!.copyOf()
        val n = main.size
        realEigenvalues = DoubleArray(n)
        imagEigenvalues = DoubleArray(n)
        val e = DoubleArray(n)
        for (i in 0 until n - 1) {
            realEigenvalues[i] = main[i]
            e[i] = secondary[i]
        }
        realEigenvalues[n - 1] = main[n - 1]
        e[n - 1] = 0.0

        // Determine the largest main and secondary value in absolute term.
        var maxAbsoluteValue = 0.0
        for (i in 0 until n) {
            if (abs(realEigenvalues[i]) > maxAbsoluteValue) {
                maxAbsoluteValue = abs(realEigenvalues[i])
            }
            if (abs(e[i]) > maxAbsoluteValue) {
                maxAbsoluteValue = abs(e[i])
            }
        }
        // Make null any main and secondary value too small to be significant
        if (maxAbsoluteValue != 0.0) {
            for (i in 0 until n) {
                if (abs(realEigenvalues[i]) <= Precision.EPSILON * maxAbsoluteValue) {
                    realEigenvalues[i] = 0.0
                }
                if (abs(e[i]) <= Precision.EPSILON * maxAbsoluteValue) {
                    e[i] = 0.0
                }
            }
        }
        for (j in 0 until n) {
            var its = 0
            var m: Int
            do {
                m = j
                while (m < n - 1) {
                    val delta = abs(realEigenvalues[m]) +
                            abs(realEigenvalues[m + 1])
                    if (abs(e[m]) + delta == delta) {
                        break
                    }
                    m++
                }
                if (m != j) {
                    if (its == maxIter.toInt()) {
                        throw MaxCountExceededException(
                            LocalizedFormats.CONVERGENCE_FAILED,
                            maxIter
                        )
                    }
                    its++
                    var q = (realEigenvalues[j + 1] - realEigenvalues[j]) / (2 * e[j])
                    var t = sqrt(1 + q * q)
                    q = if (q < 0.0) {
                        realEigenvalues[m] - realEigenvalues[j] + e[j] / (q - t)
                    } else {
                        realEigenvalues[m] - realEigenvalues[j] + e[j] / (q + t)
                    }
                    var u = 0.0
                    var s = 1.0
                    var c = 1.0
                    var i: Int
                    i = m - 1
                    while (i >= j) {
                        var p = s * e[i]
                        val h = c * e[i]
                        if (abs(p) >= abs(q)) {
                            c = q / p
                            t = sqrt(c * c + 1.0)
                            e[i + 1] = p * t
                            s = 1.0 / t
                            c *= s
                        } else {
                            s = p / q
                            t = sqrt(s * s + 1.0)
                            e[i + 1] = q * t
                            c = 1.0 / t
                            s *= c
                        }
                        if (e[i + 1] == 0.0) {
                            realEigenvalues[i + 1] -= u
                            e[m] = 0.0
                            break
                        }
                        q = realEigenvalues[i + 1] - u
                        t = (realEigenvalues[i] - q) * s + 2.0 * c * h
                        u = s * t
                        realEigenvalues[i + 1] = q + u
                        q = c * t - h
                        for (ia in 0 until n) {
                            p = z[ia]!![i + 1]
                            z[ia]!![i + 1] = s * z[ia]!![i] + c * p
                            z[ia]!![i] = c * z[ia]!![i] - s * p
                        }
                        i--
                    }
                    if (t == 0.0 && i >= j) {
                        continue
                    }
                    realEigenvalues[j] -= u
                    e[j] = q
                    e[m] = 0.0
                }
            } while (m != j)
        }

        //Sort the eigen values (and vectors) in increase order
        for (i in 0 until n) {
            var k = i
            var p = realEigenvalues[i]
            for (j in i + 1 until n) {
                if (realEigenvalues[j] > p) {
                    k = j
                    p = realEigenvalues[j]
                }
            }
            if (k != i) {
                realEigenvalues[k] = realEigenvalues[i]
                realEigenvalues[i] = p
                for (j in 0 until n) {
                    p = z[j]!![i]
                    z[j]!![i] = z[j]!![k]
                    z[j]!![k] = p
                }
            }
        }

        // Determine the largest eigen value in absolute term.
        maxAbsoluteValue = 0.0
        for (i in 0 until n) {
            if (abs(realEigenvalues[i]) > maxAbsoluteValue) {
                maxAbsoluteValue = abs(realEigenvalues[i])
            }
        }
        // Make null any eigen value too small to be significant
        if (maxAbsoluteValue != 0.0) {
            for (i in 0 until n) {
                if (abs(realEigenvalues[i]) < Precision.EPSILON * maxAbsoluteValue) {
                    realEigenvalues[i] = 0.0
                }
            }
        }
        eigenvectors = arrayOfNulls(n)
        val tmp = DoubleArray(n)
        for (i in 0 until n) {
            for (j in 0 until n) {
                tmp[j] = z[j]!![i]
            }
            eigenvectors[i] = ArrayRealVector(tmp)
        }
    }

    /**
     * Transforms the matrix to Schur form and calculates the eigenvalues.
     *
     * @param matrix Matrix to transform.
     * @return the [Shur transform][SchurTransformer] for this matrix
     */
    private fun transformToSchur(matrix: RealMatrix): SchurTransformer {
        val schurTransform = SchurTransformer(matrix)
        val matT = schurTransform.t.getData()
        realEigenvalues = DoubleArray(matT!!.size)
        imagEigenvalues = DoubleArray(matT.size)
        var i = 0
        while (i < realEigenvalues.size) {
            if (i == realEigenvalues.size - 1 ||
                equalsByEps(matT[i + 1]!![i], 0.0, EPSILON)
            ) {
                realEigenvalues[i] = matT[i]!![i]
            } else {
                val x = matT[i + 1]!![i + 1]
                val p = 0.5 * (matT[i]!![i] - x)
                val z = sqrt(abs(p * p + matT[i + 1]!![i] * matT[i]!![i + 1]))
                realEigenvalues[i] = x + p
                imagEigenvalues[i] = z
                realEigenvalues[i + 1] = x + p
                imagEigenvalues[i + 1] = -z
                i++
            }
            i++
        }
        return schurTransform
    }

    /**
     * Performs a division of two complex numbers.
     *
     * @param xr real part of the first number
     * @param xi imaginary part of the first number
     * @param yr real part of the second number
     * @param yi imaginary part of the second number
     * @return result of the complex division
     */
    private fun cdiv(
        xr: Double, xi: Double,
        yr: Double, yi: Double
    ): Complex {
        return Complex(xr, xi).divide(Complex(yr, yi))
    }

    /**
     * Find eigenvectors from a matrix transformed to Schur form.
     *
     * @param schur the schur transformation of the matrix
     * @throws MathArithmeticException if the Schur form has a norm of zero
     */
    @Throws(MathArithmeticException::class)
    private fun findEigenVectorsFromSchur(schur: SchurTransformer) {
        val matrixT = schur.t.getData()
        val matrixP = schur.p.getData()
        val n = matrixT!!.size

        // compute matrix norm
        var norm = 0.0
        for (i in 0 until n) {
            for (j in max(i - 1, 0) until n) {
                norm += abs(matrixT[i]!![j])
            }
        }

        // we can not handle a matrix with zero norm
        if (equalsByEps(norm, 0.0, EPSILON)) {
            throw MathArithmeticException(LocalizedFormats.ZERO_NORM)
        }

        // Backsubstitute to find vectors of upper triangular form
        var r = 0.0
        var s = 0.0
        var z = 0.0
        for (idx in n - 1 downTo 0) {
            val p = realEigenvalues[idx]
            var q = imagEigenvalues[idx]
            if (equals(q, 0.0)) {
                // Real vector
                var l = idx
                matrixT[idx]!![idx] = 1.0
                for (i in idx - 1 downTo 0) {
                    val w = matrixT[i]!![i] - p
                    r = 0.0
                    for (j in l..idx) {
                        r += matrixT[i]!![j] * matrixT[j]!![idx]
                    }
                    if (compareTo(imagEigenvalues[i], 0.0, EPSILON) < 0) {
                        z = w
                        s = r
                    } else {
                        l = i
                        if (equals(imagEigenvalues[i], 0.0)) {
                            if (w != 0.0) {
                                matrixT[i]!![idx] = -r / w
                            } else {
                                matrixT[i]!![idx] = -r / (Precision.EPSILON * norm)
                            }
                        } else {
                            // Solve real equations
                            val x = matrixT[i]!![i + 1]
                            val y = matrixT[i + 1]!![i]
                            q = (realEigenvalues[i] - p) * (realEigenvalues[i] - p) +
                                    imagEigenvalues[i] * imagEigenvalues[i]
                            val t = (x * s - z * r) / q
                            matrixT[i]!![idx] = t
                            if (abs(x) > abs(z)) {
                                matrixT[i + 1]!![idx] = (-r - w * t) / x
                            } else {
                                matrixT[i + 1]!![idx] = (-s - y * t) / z
                            }
                        }

                        // Overflow control
                        val t = abs(matrixT[i]!![idx])
                        if (Precision.EPSILON * t * t > 1) {
                            for (j in i..idx) {
                                matrixT[j]!![idx] /= t
                            }
                        }
                    }
                }
            } else if (q < 0.0) {
                // Complex vector
                var l = idx - 1

                // Last vector component imaginary so matrix is triangular
                if (abs(matrixT[idx]!![idx - 1]) > abs(matrixT[idx - 1]!![idx])) {
                    matrixT[idx - 1]!![idx - 1] = q / matrixT[idx]!![idx - 1]
                    matrixT[idx - 1]!![idx] = -(matrixT[idx]!![idx] - p) / matrixT[idx]!![idx - 1]
                } else {
                    val result = cdiv(
                        0.0, -matrixT[idx - 1]!![idx],
                        matrixT[idx - 1]!![idx - 1] - p, q
                    )
                    matrixT[idx - 1]!![idx - 1] = result.real
                    matrixT[idx - 1]!![idx] = result.imaginary
                }
                matrixT[idx]!![idx - 1] = 0.0
                matrixT[idx]!![idx] = 1.0
                for (i in idx - 2 downTo 0) {
                    var ra = 0.0
                    var sa = 0.0
                    for (j in l..idx) {
                        ra += matrixT[i]!![j] * matrixT[j]!![idx - 1]
                        sa += matrixT[i]!![j] * matrixT[j]!![idx]
                    }
                    val w = matrixT[i]!![i] - p
                    if (compareTo(imagEigenvalues[i], 0.0, EPSILON) < 0) {
                        z = w
                        r = ra
                        s = sa
                    } else {
                        l = i
                        if (equals(imagEigenvalues[i], 0.0)) {
                            val c = cdiv(-ra, -sa, w, q)
                            matrixT[i]!![idx - 1] = c.real
                            matrixT[i]!![idx] = c.imaginary
                        } else {
                            // Solve complex equations
                            val x = matrixT[i]!![i + 1]
                            val y = matrixT[i + 1]!![i]
                            var vr = (realEigenvalues[i] - p) * (realEigenvalues[i] - p) +
                                    imagEigenvalues[i] * imagEigenvalues[i] - q * q
                            val vi = (realEigenvalues[i] - p) * 2.0 * q
                            if (equals(vr, 0.0) && equals(vi, 0.0)) {
                                vr = Precision.EPSILON * norm *
                                        (abs(w) + abs(q) + abs(x) +
                                                abs(y) + abs(z))
                            }
                            val c = cdiv(
                                x * r - z * ra + q * sa,
                                x * s - z * sa - q * ra, vr, vi
                            )
                            matrixT[i]!![idx - 1] = c.real
                            matrixT[i]!![idx] = c.imaginary
                            if (abs(x) > abs(z) + abs(q)) {
                                matrixT[i + 1]!![idx - 1] = (-ra - w * matrixT[i]!![idx - 1] +
                                        q * matrixT[i]!![idx]) / x
                                matrixT[i + 1]!![idx] = (-sa - w * matrixT[i]!![idx] - q * matrixT[i]!![idx - 1]) / x
                            } else {
                                val c2 = cdiv(
                                    -r - y * matrixT[i]!![idx - 1],
                                    -s - y * matrixT[i]!![idx], z, q
                                )
                                matrixT[i + 1]!![idx - 1] = c2.real
                                matrixT[i + 1]!![idx] = c2.imaginary
                            }
                        }

                        // Overflow control
                        val t = max(
                            abs(matrixT[i]!![idx - 1]),
                            abs(matrixT[i]!![idx])
                        )
                        if (Precision.EPSILON * t * t > 1) {
                            for (j in i..idx) {
                                matrixT[j]!![idx - 1] /= t
                                matrixT[j]!![idx] /= t
                            }
                        }
                    }
                }
            }
        }

        // Back transformation to get eigenvectors of original matrix
        for (j in n - 1 downTo 0) {
            for (i in 0..n - 1) {
                z = 0.0
                for (k in 0..min(j, n - 1)) {
                    z += matrixP!![i]!![k] * matrixT[k]!![j]
                }
                matrixP!![i]!![j] = z
            }
        }
        eigenvectors = arrayOfNulls(n)
        val tmp = DoubleArray(n)
        for (i in 0 until n) {
            for (j in 0 until n) {
                tmp[j] = matrixP!![j]!![i]
            }
            eigenvectors[i] = ArrayRealVector(tmp)
        }
    }

    /**
     * Specialized solver.
     */
    private class Solver
    /**
     * Builds a solver from decomposed matrix.
     *
     * @param realEigenvalues Real parts of the eigenvalues.
     * @param imagEigenvalues Imaginary parts of the eigenvalues.
     * @param eigenvectors    Eigenvectors.
     */(
        /**
         * Real part of the realEigenvalues.
         */
        private val realEigenvalues: DoubleArray,
        /**
         * Imaginary part of the realEigenvalues.
         */
        private val imagEigenvalues: DoubleArray,
        /**
         * Eigenvectors.
         */
        private val eigenvectors: Array<ArrayRealVector?>
    ) : DecompositionSolver {
        /**
         * Solves the linear equation A  X = B for symmetric matrices A.
         *
         *
         * This method only finds exact linear solutions, i.e. solutions for
         * which ||A  X - B|| is exactly 0.
         *
         *
         * @param b Right-hand side of the equation A  X = B.
         * @return a Vector X that minimizes the two norm of A  X - B.
         * @throws DimensionMismatchException if the matrices dimensions do not match.
         * @throws SingularMatrixException    if the decomposed matrix is singular.
         */
        override fun solve(b: RealVector?): RealVector? {
            if (!isNonSingular) {
                throw SingularMatrixException()
            }
            val m = realEigenvalues.size
            if (b!!.getDimension() != m) {
                throw DimensionMismatchException(b.getDimension(), m)
            }
            val bp = DoubleArray(m)
            for (i in 0 until m) {
                val v = eigenvectors[i]
                val vData = v!!.dataRef
                val s = v.dotProduct(b) / realEigenvalues[i]
                for (j in 0 until m) {
                    bp[j] += s * vData[j]
                }
            }
            return ArrayRealVector(bp, false)
        }

        /**
         * {@inheritDoc}
         */
        override fun solve(b: RealMatrix?): RealMatrix? {
            if (!isNonSingular) {
                throw SingularMatrixException()
            }
            val m = realEigenvalues.size
            if (b!!.getRowDimension() != m) {
                throw DimensionMismatchException(b.getRowDimension(), m)
            }
            val nColB = b.getColumnDimension()
            val bp = Array<DoubleArray?>(m) { DoubleArray(nColB) }
            val tmpCol = DoubleArray(m)
            for (k in 0 until nColB) {
                for (i in 0 until m) {
                    tmpCol[i] = b.getEntry(i, k)
                    bp[i]!![k] = 0.0
                }
                for (i in 0 until m) {
                    val v = eigenvectors[i]
                    val vData = v!!.dataRef
                    var s = 0.0
                    for (j in 0 until m) {
                        s += v.getEntry(j) * tmpCol[j]
                    }
                    s /= realEigenvalues[i]
                    for (j in 0 until m) {
                        bp[j]!![k] += s * vData[j]
                    }
                }
            }
            return Array2DRowRealMatrix(bp, false)
        }// Looking for eigenvalues that are 0, where we consider anything much much smaller
        // than the largest eigenvalue to be effectively 0.
// Looping over all values (in case they are not sorted in decreasing
        // order of their norm).
        // Corner case: zero matrix, all exactly 0 eigenvalues
        /**
         * Checks whether the decomposed matrix is non-singular.
         *
         * @return true if the decomposed matrix is non-singular.
         */
        override val isNonSingular: Boolean
            get() {
                var largestEigenvalueNorm = 0.0
                // Looping over all values (in case they are not sorted in decreasing
                // order of their norm).
                for (i in realEigenvalues.indices) {
                    largestEigenvalueNorm = max(largestEigenvalueNorm, eigenvalueNorm(i))
                }
                // Corner case: zero matrix, all exactly 0 eigenvalues
                if (largestEigenvalueNorm == 0.0) {
                    return false
                }
                for (i in realEigenvalues.indices) {
                    // Looking for eigenvalues that are 0, where we consider anything much much smaller
                    // than the largest eigenvalue to be effectively 0.
                    if (equalsByEps(eigenvalueNorm(i) / largestEigenvalueNorm, 0.0, EPSILON)) {
                        return false
                    }
                }
                return true
            }

        /**
         * @param i which eigenvalue to find the norm of
         * @return the norm of ith (complex) eigenvalue.
         */
        private fun eigenvalueNorm(i: Int): Double {
            val re = realEigenvalues[i]
            val im = imagEigenvalues[i]
            return sqrt(re * re + im * im)
        }

        /**
         * Get the inverse of the decomposed matrix.
         *
         * @return the inverse matrix.
         * @throws SingularMatrixException if the decomposed matrix is singular.
         */
        override fun getInverse(): RealMatrix? {
            if (!isNonSingular) {
                throw SingularMatrixException()
            }
            val m = realEigenvalues.size
            val invData = Array(m) { DoubleArray(m) }
            for (i in 0 until m) {
                val invI = invData[i]
                for (j in 0 until m) {
                    var invIJ = 0.0
                    for (k in 0 until m) {
                        val vK = eigenvectors[k]!!.dataRef
                        invIJ += vK[i] * vK[j] / realEigenvalues[k]
                    }
                    invI[j] = invIJ
                }
            }
            return MatrixUtils.createRealMatrix(invData as Array<DoubleArray?>?)
        }
    }

    companion object {
        /**
         * Internally used epsilon criteria.
         */
        private const val EPSILON = 1e-12
    }
}