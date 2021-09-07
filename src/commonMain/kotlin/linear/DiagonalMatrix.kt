package linear

import exception.DimensionMismatchException
import exception.NotStrictlyPositiveException
import exception.NumberIsTooLargeException
import exception.OutOfRangeException
import linear.exception.MatrixDimensionMismatchException
import linear.exception.SingularMatrixException
import util.FastMath.abs
import util.MathUtils.checkNotNull
import util.Precision.equalsByEps
import util.Precision.equalsByMaxUlps

/**
 * Implementation of a diagonal matrix.
 *
 * @author haokangkang
 * @since 3.1.1
 */
class DiagonalMatrix : AbstractRealMatrix {
    /**
     * Gets a reference to the underlying data array.
     *
     * @return 1-dimensional array of entries.
     */
    /**
     * Entries of the diagonal.
     */
    val dataRef: DoubleArray

    /**
     * Creates a matrix with the supplied dimension.
     *
     * @param dimension Number of rows and columns in the new matrix.
     * @throws NotStrictlyPositiveException if the dimension is
     * not positive.
     */
    constructor(dimension: Int) : super(dimension, dimension) {
        dataRef = DoubleArray(dimension)
    }
    /**
     * Creates a matrix using the input array as the underlying data.
     * <br></br>
     * If an array is created specially in order to be embedded in a
     * this instance and not used directly, the `copyArray` may be
     * set to `false`.
     * This will prevent the copying and improve performance as no new
     * array will be built and no data will be copied.
     *
     * @param d         Data for new matrix.
     * @param copyArray if `true`, the input array will be copied,
     * otherwise it will be referenced.
     * @throws NullArgumentException if d is null
     */
    /**
     * Creates a matrix using the input array as the underlying data.
     * <br></br>
     * The input array is copied, not referenced.
     *
     * @param d Data for the new matrix.
     */
    constructor(d: DoubleArray?, copyArray: Boolean = true) {
        checkNotNull(d)
        dataRef = if (copyArray) d!!.copyOf() else d!!
    }

    /**
     * {@inheritDoc}
     *
     * @throws DimensionMismatchException if the requested dimensions are not equal.
     */
    override fun createMatrix(
        rowDimension: Int,
        columnDimension: Int
    ): RealMatrix? {
        if (rowDimension != columnDimension) {
            throw DimensionMismatchException(rowDimension, columnDimension)
        }
        return DiagonalMatrix(rowDimension)
    }

    /**
     * {@inheritDoc}
     */
    override fun copy(): RealMatrix? {
        return DiagonalMatrix(dataRef)
    }

    /**
     * Compute the sum of `this` and `m`.
     *
     * @param m Matrix to be added.
     * @return `this + m`.
     * @throws MatrixDimensionMismatchException if `m` is not the same
     * size as `this`.
     */
    @Throws(MatrixDimensionMismatchException::class)
    fun add(m: DiagonalMatrix): DiagonalMatrix {
        // Safety check.
        MatrixUtils.checkAdditionCompatible(this, m)
        val dim = getRowDimension()
        val outData = DoubleArray(dim)
        for (i in 0 until dim) {
            outData[i] = dataRef[i] + m.dataRef[i]
        }
        return DiagonalMatrix(outData, false)
    }

    /**
     * Returns `this` minus `m`.
     *
     * @param m Matrix to be subtracted.
     * @return `this - m`
     * @throws MatrixDimensionMismatchException if `m` is not the same
     * size as `this`.
     */
    @Throws(MatrixDimensionMismatchException::class)
    fun subtract(m: DiagonalMatrix): DiagonalMatrix {
        MatrixUtils.checkSubtractionCompatible(this, m)
        val dim = getRowDimension()
        val outData = DoubleArray(dim)
        for (i in 0 until dim) {
            outData[i] = dataRef[i] - m.dataRef[i]
        }
        return DiagonalMatrix(outData, false)
    }

    /**
     * Returns the result of postmultiplying `this` by `m`.
     *
     * @param m matrix to postmultiply by
     * @return `this * m`
     * @throws DimensionMismatchException if
     * `columnDimension(this) != rowDimension(m)`
     */
    @Throws(DimensionMismatchException::class)
    fun multiply(m: DiagonalMatrix): DiagonalMatrix {
        MatrixUtils.checkMultiplicationCompatible(this, m)
        val dim = getRowDimension()
        val outData = DoubleArray(dim)
        for (i in 0 until dim) {
            outData[i] = dataRef[i] * m.dataRef[i]
        }
        return DiagonalMatrix(outData, false)
    }

    /**
     * Returns the result of postmultiplying `this` by `m`.
     *
     * @param m matrix to postmultiply by
     * @return `this * m`
     * @throws DimensionMismatchException if
     * `columnDimension(this) != rowDimension(m)`
     */
    @Throws(DimensionMismatchException::class)
    override fun multiply(m: RealMatrix?): RealMatrix? {
        return if (m is DiagonalMatrix) {
            multiply(m)
        } else {
            MatrixUtils.checkMultiplicationCompatible(this, m as AnyMatrix)
            val nRows = m.getRowDimension()
            val nCols = m.getColumnDimension()
            val product = Array<DoubleArray?>(nRows) { DoubleArray(nCols) }
            for (r in 0 until nRows) {
                for (c in 0 until nCols) {
                    product[r]!![c] = dataRef[r] * m.getEntry(r, c)
                }
            }
            Array2DRowRealMatrix(product, false)
        }
    }

    /**
     * {@inheritDoc}
     */
    override fun getData(): Array<DoubleArray?>? {
        val dim = getRowDimension()
        val out = Array<DoubleArray?>(dim) { DoubleArray(dim) }
        for (i in 0 until dim) {
            out[i]!![i] = dataRef[i]
        }
        return out
    }

    /**
     * {@inheritDoc}
     */
    override fun getEntry(row: Int, column: Int): Double {
        MatrixUtils.checkMatrixIndex(this, row, column)
        return if (row == column) dataRef[row] else 0.0
    }

    /**
     * {@inheritDoc}
     *
     * @throws NumberIsTooLargeException if `row != column` and value is non-zero.
     */
    override fun setEntry(row: Int, column: Int, value: Double) {
        if (row == column) {
            MatrixUtils.checkRowIndex(this, row)
            dataRef[row] = value
        } else {
            ensureZero(value)
        }
    }

    /**
     * {@inheritDoc}
     *
     * @throws NumberIsTooLargeException if `row != column` and increment is non-zero.
     */
    override fun addToEntry(
        row: Int,
        column: Int,
        increment: Double
    ) {
        if (row == column) {
            MatrixUtils.checkRowIndex(this, row)
            dataRef[row] += increment
        } else {
            ensureZero(increment)
        }
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class)
    override fun multiplyEntry(
        row: Int,
        column: Int,
        factor: Double
    ) {
        // we don't care about non-diagonal elements for multiplication
        if (row == column) {
            MatrixUtils.checkRowIndex(this, row)
            dataRef[row] *= factor
        }
    }

    /**
     * {@inheritDoc}
     */
    override fun getRowDimension(): Int {
        return dataRef.size
    }

    /**
     * {@inheritDoc}
     */
    override fun getColumnDimension(): Int {
        return dataRef.size
    }

    /**
     * {@inheritDoc}
     */
    @Throws(DimensionMismatchException::class)
    override fun operate(v: DoubleArray?): DoubleArray? {
        return multiply(DiagonalMatrix(v, false)).dataRef
    }

    /**
     * {@inheritDoc}
     */
    @Throws(DimensionMismatchException::class)
    override fun preMultiply(v: DoubleArray): DoubleArray? {
        return operate(v)
    }

    /**
     * {@inheritDoc}
     */
    @Throws(DimensionMismatchException::class)
    override fun preMultiply(v: RealVector): RealVector? {
        val vectorData: DoubleArray
        vectorData = if (v is ArrayRealVector) {
            v.dataRef
        } else {
            v.toArray()
        }
        return MatrixUtils.createRealVector(preMultiply(vectorData))
    }

    /**
     * Ensure a value is zero.
     *
     * @param value value to check
     * @throws NumberIsTooLargeException if value is not zero
     */
    @Throws(NumberIsTooLargeException::class)
    private fun ensureZero(value: Double) {
        if (!equalsByMaxUlps(0.0, value, 1)) {
            throw NumberIsTooLargeException(abs(value), 0, true)
        }
    }
    /**
     * Computes the inverse of this diagonal matrix.
     *
     * @param threshold Singularity threshold.
     * @return the inverse of `m`
     * @throws SingularMatrixException if the matrix is singular
     * @since 3.3
     */
    /**
     * Computes the inverse of this diagonal matrix.
     *
     *
     * Note: this method will use a singularity threshold of 0,
     * use [.inverse] if a different threshold is needed.
     *
     * @return the inverse of `m`
     * @throws SingularMatrixException if the matrix is singular
     * @since 3.3
     */

    @Throws(SingularMatrixException::class)
    fun inverse(threshold: Double = 0.0): DiagonalMatrix {
        if (isSingular(threshold)) {
            throw SingularMatrixException()
        }
        val result = DoubleArray(dataRef.size)
        for (i in dataRef.indices) {
            result[i] = 1.0 / dataRef[i]
        }
        return DiagonalMatrix(result, false)
    }

    /**
     * Returns whether this diagonal matrix is singular, i.e. any diagonal entry
     * is equal to `0` within the given threshold.
     *
     * @param threshold Singularity threshold.
     * @return `true` if the matrix is singular, `false` otherwise
     * @since 3.3
     */
    fun isSingular(threshold: Double): Boolean {
        for (i in dataRef.indices) {
            if (equalsByEps(dataRef[i], 0.0, threshold)) {
                return true
            }
        }
        return false
    }
}