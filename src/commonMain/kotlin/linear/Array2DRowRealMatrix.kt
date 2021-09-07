package linear

import exception.*
import exception.util.LocalizedFormats
import linear.exception.MatrixDimensionMismatchException
import linear.visitor.RealMatrixChangingVisitor
import linear.visitor.RealMatrixPreservingVisitor
import util.MathUtils.checkNotNull
import crossjvm.ArrayUtil

/**
 * Implementation of [RealMatrix] using a `double[][]` array to
 * store entries.
 *
 * @author haokangkang
 */
class Array2DRowRealMatrix : AbstractRealMatrix {
    /**
     * Get a reference to the underlying data array.
     *
     * @return 2-dimensional array of entries.
     */
    /**
     * Entries of the matrix.
     */
    var dataRef: Array<DoubleArray?>? = null
        private set

    /**
     * Create a new RealMatrix with the supplied row and column dimensions.
     *
     * @param rowDimension    Number of rows in the new matrix.
     * @param columnDimension Number of columns in the new matrix.
     * @throws NotStrictlyPositiveException if the row or column dimension is
     * not positive.
     */
    constructor(
        rowDimension: Int,
        columnDimension: Int
    ) : super(rowDimension, columnDimension) {
        dataRef = Array(rowDimension) { DoubleArray(columnDimension) }
    }

    /**
     * Create a new `RealMatrix` using the input array as the underlying
     * data array.
     *
     * The input array is copied, not referenced. This constructor has
     * the same effect as calling [.Array2DRowRealMatrix]
     * with the second argument set to `true`.
     *
     * @param d Data for the new matrix.
     * @throws DimensionMismatchException if `d` is not rectangular.
     * @throws NoDataException            if `d` row or column dimension is zero.
     * @throws NullArgumentException      if `d` is `null`.
     * @see .Array2DRowRealMatrix
     */
    constructor(d: Array<DoubleArray?>?) {
        copyIn(d)
    }

    /**
     * Create a new RealMatrix using the input array as the underlying
     * data array.
     * If an array is built specially in order to be embedded in a
     * RealMatrix and not used directly, the `copyArray` may be
     * set to `false`. This will prevent the copying and improve
     * performance as no new array will be built and no data will be copied.
     *
     * @param d         Data for new matrix.
     * @param copyArray if `true`, the input array will be copied,
     * otherwise it will be referenced.
     * @throws DimensionMismatchException if `d` is not rectangular.
     * @throws NoDataException            if `d` row or column dimension is zero.
     * @throws NullArgumentException      if `d` is `null`.
     * @see .Array2DRowRealMatrix
     */
    constructor(d: Array<DoubleArray?>?, copyArray: Boolean) {
        if (copyArray) {
            copyIn(d)
        } else {
            if (d == null) {
                throw NullArgumentException()
            }
            val nRows = d.size
            if (nRows == 0) {
                throw NoDataException(LocalizedFormats.AT_LEAST_ONE_ROW)
            }
            val nCols: Int = d[0]!!.size
            if (nCols == 0) {
                throw NoDataException(LocalizedFormats.AT_LEAST_ONE_COLUMN)
            }
            for (r in 1 until nRows) {
                if (d[r]!!.size != nCols) {
                    throw DimensionMismatchException(d[r]!!.size, nCols)
                }
            }
            dataRef = d
        }
    }

    /**
     * Create a new (column) RealMatrix using `v` as the
     * data for the unique column of the created matrix.
     * The input array is copied.
     *
     * @param v Column vector holding data for new matrix.
     */
    constructor(v: DoubleArray) {
        val nRows = v.size
        dataRef = Array(nRows) { DoubleArray(1) }
        for (row in 0 until nRows) {
            dataRef!![row]!![0] = v[row]
        }
    }

    /**
     * {@inheritDoc}
     */
    @Throws(NotStrictlyPositiveException::class)
    override fun createMatrix(
        rowDimension: Int,
        columnDimension: Int
    ): RealMatrix? {
        return Array2DRowRealMatrix(rowDimension, columnDimension)
    }

    /**
     * {@inheritDoc}
     */
    override fun copy(): RealMatrix? {
        return Array2DRowRealMatrix(copyOut(), false)
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
    fun add(m: Array2DRowRealMatrix): Array2DRowRealMatrix {
        // Safety check.
        MatrixUtils.checkAdditionCompatible(this, m)
        val rowCount = getRowDimension()
        val columnCount = getColumnDimension()
        val outData = Array<DoubleArray?>(rowCount) { DoubleArray(columnCount) }
        for (row in 0 until rowCount) {
            val dataRow = dataRef!![row]
            val mRow = m.dataRef!![row]
            val outDataRow = outData[row]
            for (col in 0 until columnCount) {
                outDataRow!![col] = dataRow!![col] + mRow!![col]
            }
        }
        return Array2DRowRealMatrix(outData, false)
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
    fun subtract(m: Array2DRowRealMatrix): Array2DRowRealMatrix {
        MatrixUtils.checkSubtractionCompatible(this, m)
        val rowCount = getRowDimension()
        val columnCount = getColumnDimension()
        val outData = Array<DoubleArray?>(rowCount) { DoubleArray(columnCount) }
        for (row in 0 until rowCount) {
            val dataRow = dataRef!![row]
            val mRow = m.dataRef!![row]
            val outDataRow = outData[row]
            for (col in 0 until columnCount) {
                outDataRow!![col] = dataRow!![col] - mRow!![col]
            }
        }
        return Array2DRowRealMatrix(outData, false)
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
    fun multiply(m: Array2DRowRealMatrix): Array2DRowRealMatrix {
        MatrixUtils.checkMultiplicationCompatible(this, m)
        val nRows = getRowDimension()
        val nCols = m.getColumnDimension()
        val nSum = getColumnDimension()
        val outData = Array<DoubleArray?>(nRows) { DoubleArray(nCols) }
        // Will hold a column of "m".
        val mCol = DoubleArray(nSum)
        val mData = m.dataRef

        // Multiply.
        for (col in 0 until nCols) {
            // Copy all elements of column "col" of "m" so that
            // will be in contiguous memory.
            for (mRow in 0 until nSum) {
                mCol[mRow] = mData!![mRow]!![col]
            }
            for (row in 0 until nRows) {
                val dataRow = dataRef!![row]
                var sum = 0.0
                for (i in 0 until nSum) {
                    sum += dataRow!![i] * mCol[i]
                }
                outData[row]!![col] = sum
            }
        }
        return Array2DRowRealMatrix(outData, false)
    }

    /**
     * {@inheritDoc}
     */
    override fun getData(): Array<DoubleArray?>? {
        return copyOut()
    }

    /**
     * {@inheritDoc}
     */
    @Throws(
        NoDataException::class,
        OutOfRangeException::class,
        DimensionMismatchException::class,
        NullArgumentException::class
    )
    override fun setSubMatrix(
        subMatrix: Array<DoubleArray?>?,
        row: Int,
        column: Int
    ) {
        if (dataRef == null) {
            if (row > 0) {
                throw MathIllegalStateException(LocalizedFormats.FIRST_ROWS_NOT_INITIALIZED_YET, row)
            }
            if (column > 0) {
                throw MathIllegalStateException(LocalizedFormats.FIRST_COLUMNS_NOT_INITIALIZED_YET, column)
            }
            checkNotNull(subMatrix)
            val nRows = subMatrix!!.size
            if (nRows == 0) {
                throw NoDataException(LocalizedFormats.AT_LEAST_ONE_ROW)
            }
            val nCols: Int = subMatrix[0]!!.size
            if (nCols == 0) {
                throw NoDataException(LocalizedFormats.AT_LEAST_ONE_COLUMN)
            }
            dataRef = Array(subMatrix.size) { DoubleArray(nCols) }
            for (i in dataRef!!.indices) {
                if (subMatrix[i]!!.size != nCols) {
                    throw DimensionMismatchException(subMatrix[i]!!.size, nCols)
                }
                ArrayUtil.arraycopy(subMatrix[i], 0, dataRef!![i + row], column, nCols)
            }
        } else {
            super.setSubMatrix(subMatrix, row, column)
        }
    }

    /**
     * {@inheritDoc}
     */
    override fun getEntry(row: Int, column: Int): Double {
        MatrixUtils.checkMatrixIndex(this, row, column)
        return dataRef!![row]!![column]
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class)
    override fun setEntry(row: Int, column: Int, value: Double) {
        MatrixUtils.checkMatrixIndex(this, row, column)
        dataRef!![row]!![column] = value
    }

    /**
     * {@inheritDoc}
     */
    override fun addToEntry(
        row: Int, column: Int,
        increment: Double
    ) {
        MatrixUtils.checkMatrixIndex(this, row, column)
        dataRef!![row]!![column] += increment
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class)
    override fun multiplyEntry(
        row: Int, column: Int,
        factor: Double
    ) {
        MatrixUtils.checkMatrixIndex(this, row, column)
        dataRef!![row]!![column] *= factor
    }

    /**
     * {@inheritDoc}
     */
    override fun getRowDimension(): Int {
        return if (dataRef == null) 0 else dataRef!!.size
    }

    /**
     * {@inheritDoc}
     */
    override fun getColumnDimension(): Int {
        return if (dataRef == null || dataRef!![0] == null) 0 else dataRef!![0]!!.size
    }

    /**
     * {@inheritDoc}
     */
    @Throws(DimensionMismatchException::class)
    override fun operate(v: DoubleArray?): DoubleArray? {
        val nRows = getRowDimension()
        val nCols = getColumnDimension()
        if (v!!.size != nCols) {
            throw DimensionMismatchException(v.size, nCols)
        }
        val out = DoubleArray(nRows)
        for (row in 0 until nRows) {
            val dataRow = dataRef!![row]
            var sum = 0.0
            for (i in 0 until nCols) {
                sum += dataRow!![i] * v[i]
            }
            out[row] = sum
        }
        return out
    }

    /**
     * {@inheritDoc}
     */
    @Throws(DimensionMismatchException::class)
    override fun preMultiply(v: DoubleArray): DoubleArray? {
        val nRows = getRowDimension()
        val nCols = getColumnDimension()
        if (v.size != nRows) {
            throw DimensionMismatchException(v.size, nRows)
        }
        val out = DoubleArray(nCols)
        for (col in 0 until nCols) {
            var sum = 0.0
            for (i in 0 until nRows) {
                sum += dataRef!![i]!![col] * v[i]
            }
            out[col] = sum
        }
        return out
    }

    /**
     * {@inheritDoc}
     */
    override fun walkInRowOrder(visitor: RealMatrixChangingVisitor): Double {
        val rows = getRowDimension()
        val columns = getColumnDimension()
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1)
        for (i in 0 until rows) {
            val rowI = dataRef!![i]
            for (j in 0 until columns) {
                rowI!![j] = visitor.visit(i, j, rowI[j])
            }
        }
        return visitor.end()
    }

    /**
     * {@inheritDoc}
     */
    override fun walkInRowOrder(visitor: RealMatrixPreservingVisitor?): Double {
        val rows = getRowDimension()
        val columns = getColumnDimension()
        visitor!!.start(rows, columns, 0, rows - 1, 0, columns - 1)
        for (i in 0 until rows) {
            val rowI = dataRef!![i]
            for (j in 0 until columns) {
                visitor.visit(i, j, rowI!![j])
            }
        }
        return visitor.end()
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class, NumberIsTooSmallException::class)
    override fun walkInRowOrder(
        visitor: RealMatrixChangingVisitor,
        startRow: Int,
        endRow: Int,
        startColumn: Int,
        endColumn: Int
    ): Double {
        MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn)
        visitor.start(
            getRowDimension(), getColumnDimension(),
            startRow, endRow, startColumn, endColumn
        )
        for (i in startRow..endRow) {
            val rowI = dataRef!![i]
            for (j in startColumn..endColumn) {
                rowI!![j] = visitor.visit(i, j, rowI[j])
            }
        }
        return visitor.end()
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class, NumberIsTooSmallException::class)
    override fun walkInRowOrder(
        visitor: RealMatrixPreservingVisitor?,
        startRow: Int,
        endRow: Int,
        startColumn: Int,
        endColumn: Int
    ): Double {
        MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn)
        visitor!!.start(
            getRowDimension(), getColumnDimension(),
            startRow, endRow, startColumn, endColumn
        )
        for (i in startRow..endRow) {
            val rowI = dataRef!![i]
            for (j in startColumn..endColumn) {
                visitor.visit(i, j, rowI!![j])
            }
        }
        return visitor.end()
    }

    /**
     * {@inheritDoc}
     */
    override fun walkInColumnOrder(visitor: RealMatrixChangingVisitor): Double {
        val rows = getRowDimension()
        val columns = getColumnDimension()
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1)
        for (j in 0 until columns) {
            for (i in 0 until rows) {
                val rowI = dataRef!![i]
                rowI!![j] = visitor.visit(i, j, rowI[j])
            }
        }
        return visitor.end()
    }

    /**
     * {@inheritDoc}
     */
    override fun walkInColumnOrder(visitor: RealMatrixPreservingVisitor): Double {
        val rows = getRowDimension()
        val columns = getColumnDimension()
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1)
        for (j in 0 until columns) {
            for (i in 0 until rows) {
                visitor.visit(i, j, dataRef!![i]!![j])
            }
        }
        return visitor.end()
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class, NumberIsTooSmallException::class)
    override fun walkInColumnOrder(
        visitor: RealMatrixChangingVisitor,
        startRow: Int,
        endRow: Int,
        startColumn: Int,
        endColumn: Int
    ): Double {
        MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn)
        visitor.start(
            getRowDimension(), getColumnDimension(),
            startRow, endRow, startColumn, endColumn
        )
        for (j in startColumn..endColumn) {
            for (i in startRow..endRow) {
                val rowI = dataRef!![i]
                rowI!![j] = visitor.visit(i, j, rowI[j])
            }
        }
        return visitor.end()
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class, NumberIsTooSmallException::class)
    override fun walkInColumnOrder(
        visitor: RealMatrixPreservingVisitor,
        startRow: Int,
        endRow: Int,
        startColumn: Int,
        endColumn: Int
    ): Double {
        MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn)
        visitor.start(
            getRowDimension(), getColumnDimension(),
            startRow, endRow, startColumn, endColumn
        )
        for (j in startColumn..endColumn) {
            for (i in startRow..endRow) {
                visitor.visit(i, j, dataRef!![i]!![j])
            }
        }
        return visitor.end()
    }

    /**
     * Get a fresh copy of the underlying data array.
     *
     * @return a copy of the underlying data array.
     */
    private fun copyOut(): Array<DoubleArray?> {
        val nRows = getRowDimension()
        val out = Array<DoubleArray?>(nRows) { DoubleArray(getColumnDimension()) }
        // can't copy 2-d array in one shot, otherwise get row references
        for (i in 0 until nRows) {
            ArrayUtil.arraycopy(dataRef!![i], 0, out[i], 0, dataRef!![i]!!.size)
        }
        return out
    }

    /**
     * Replace data with a fresh copy of the input array.
     *
     * @param in Data to copy.
     * @throws NoDataException            if the input array is empty.
     * @throws DimensionMismatchException if the input array is not rectangular.
     * @throws NullArgumentException      if the input array is `null`.
     */
    @Throws(DimensionMismatchException::class, NoDataException::class, NullArgumentException::class)
    private fun copyIn(`in`: Array<DoubleArray?>?) {
        setSubMatrix(`in`, 0, 0)
    }
}