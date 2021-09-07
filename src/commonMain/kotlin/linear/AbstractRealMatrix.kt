package linear

import exception.*
import exception.util.LocalizedFormats
import linear.exception.MatrixDimensionMismatchException
import linear.exception.NonSquareMatrixException
import linear.visitor.DefaultRealMatrixPreservingVisitor
import linear.visitor.RealMatrixChangingVisitor
import linear.visitor.RealMatrixPreservingVisitor
import util.FastMath.abs
import util.FastMath.max
import util.FastMath.sqrt
import util.MathUtils.checkNotNull
import util.MathUtils.hash

/**
 * Basic implementation of RealMatrix methods regardless of the underlying storage.
 *
 * All the methods implemented here use [.getEntry] to access
 * matrix elements. Derived class can provide faster implementations.
 *
 * @author haokangkang
 * @since 2.0
 */
abstract class AbstractRealMatrix : RealLinearOperator, RealMatrix {
    companion object {
        /**
         * Default format.
         */
        private val DEFAULT_FORMAT = RealMatrixFormat.instance

    }

    /**
     * Creates a matrix with no data
     */
    protected constructor()

    /**
     * Create a new RealMatrix with the supplied row and column dimensions.
     *
     * @param rowDimension    the number of rows in the new matrix
     * @param columnDimension the number of columns in the new matrix
     * @throws NotStrictlyPositiveException if row or column dimension is not positive
     */
    protected constructor(
        rowDimension: Int,
        columnDimension: Int
    ) {
        if (rowDimension < 1) {
            throw NotStrictlyPositiveException(rowDimension)
        }
        if (columnDimension < 1) {
            throw NotStrictlyPositiveException(columnDimension)
        }
    }

    /**
     * {@inheritDoc}
     */
    @Throws(MatrixDimensionMismatchException::class)
    override fun add(m: RealMatrix?): RealMatrix? {
        MatrixUtils.checkAdditionCompatible(this, m as RealMatrix)
        val rowCount = getRowDimension()
        val columnCount = getColumnDimension()
        val out = createMatrix(rowCount, columnCount)
        for (row in 0 until rowCount) {
            for (col in 0 until columnCount) {
                out!!.setEntry(row, col, getEntry(row, col) + m.getEntry(row, col))
            }
        }
        return out
    }

    /**
     * {@inheritDoc}
     */
    @Throws(MatrixDimensionMismatchException::class)
    override fun subtract(m: RealMatrix?): RealMatrix? {
        MatrixUtils.checkSubtractionCompatible(this, m as RealMatrix)
        val rowCount = getRowDimension()
        val columnCount = getColumnDimension()
        val out = createMatrix(rowCount, columnCount)
        for (row in 0 until rowCount) {
            for (col in 0 until columnCount) {
                out!!.setEntry(row, col, getEntry(row, col) - m.getEntry(row, col))
            }
        }
        return out
    }

    /**
     * {@inheritDoc}
     */
    override fun scalarAdd(d: Double): RealMatrix? {
        val rowCount = getRowDimension()
        val columnCount = getColumnDimension()
        val out = createMatrix(rowCount, columnCount)
        for (row in 0 until rowCount) {
            for (col in 0 until columnCount) {
                out!!.setEntry(row, col, getEntry(row, col) + d)
            }
        }
        return out
    }

    /**
     * {@inheritDoc}
     */
    override fun scalarMultiply(d: Double): RealMatrix? {
        val rowCount = getRowDimension()
        val columnCount = getColumnDimension()
        val out = createMatrix(rowCount, columnCount)
        for (row in 0 until rowCount) {
            for (col in 0 until columnCount) {
                out!!.setEntry(row, col, getEntry(row, col) * d)
            }
        }
        return out
    }

    /**
     * {@inheritDoc}
     */
    @Throws(DimensionMismatchException::class)
    override fun multiply(m: RealMatrix?): RealMatrix? {
        MatrixUtils.checkMultiplicationCompatible(this, m as RealMatrix)
        val nRows = getRowDimension()
        val nCols = m.getColumnDimension()
        val nSum = getColumnDimension()
        val out = createMatrix(nRows, nCols)
        for (row in 0 until nRows) {
            for (col in 0 until nCols) {
                var sum = 0.0
                for (i in 0 until nSum) {
                    sum += getEntry(row, i) * m.getEntry(i, col)
                }
                out!!.setEntry(row, col, sum)
            }
        }
        return out
    }

    /**
     * {@inheritDoc}
     */
    @Throws(DimensionMismatchException::class)
    override fun preMultiply(m: RealMatrix?): RealMatrix? {
        return m!!.multiply(this)
    }

    /**
     * {@inheritDoc}
     */
    override fun getData(): Array<DoubleArray?>? {
        val data = Array<DoubleArray?>(getRowDimension()) { DoubleArray(getColumnDimension()) }
        for (i in data.indices) {
            val dataI = data[i]
            for (j in dataI!!.indices) {
                dataI[j] = getEntry(i, j)
            }
        }
        return data
    }

    /**
     * {@inheritDoc}
     */
    override fun getNorm(): Double {
        return walkInColumnOrder(object : RealMatrixPreservingVisitor {
            /** Last row index.  */
            private var endRow = 0.0

            /** Sum of absolute values on one column.  */
            private var columnSum = 0.0

            /** Maximal sum across all columns.  */
            private var maxColSum = 0.0

            /** {@inheritDoc}  */
            override fun start(
                rows: Int, columns: Int,
                startRow: Int, endRow: Int,
                startColumn: Int, endColumn: Int
            ) {
                this.endRow = endRow.toDouble()
                columnSum = 0.0
                maxColSum = 0.0
            }

            /** {@inheritDoc}  */
            override fun visit(row: Int, column: Int, value: Double) {
                columnSum += abs(value)
                if (row.toDouble() == endRow) {
                    maxColSum = max(maxColSum, columnSum)
                    columnSum = 0.0
                }
            }

            /** {@inheritDoc}  */
            override fun end(): Double {
                return maxColSum
            }
        })
    }

    /**
     * {@inheritDoc}
     */
    override fun getFrobeniusNorm(): Double {
        return walkInOptimizedOrder(object : RealMatrixPreservingVisitor {
            /** Sum of squared entries.  */
            private var sum = 0.0

            /** {@inheritDoc}  */
            override fun start(
                rows: Int, columns: Int,
                startRow: Int, endRow: Int,
                startColumn: Int, endColumn: Int
            ) {
                sum = 0.0
            }

            /** {@inheritDoc}  */
            override fun visit(row: Int, column: Int, value: Double) {
                sum += value * value
            }

            /** {@inheritDoc}  */
            override fun end(): Double {
                return sqrt(sum)
            }
        })
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class, NumberIsTooSmallException::class)
    override fun getSubMatrix(
        startRow: Int, endRow: Int,
        startColumn: Int, endColumn: Int
    ): RealMatrix? {
        MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn)
        val subMatrix = createMatrix(endRow - startRow + 1, endColumn - startColumn + 1)
        for (i in startRow..endRow) {
            for (j in startColumn..endColumn) {
                subMatrix!!.setEntry(i - startRow, j - startColumn, getEntry(i, j))
            }
        }
        return subMatrix
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class, NumberIsTooSmallException::class, MatrixDimensionMismatchException::class)
    override fun copySubMatrix(
        startRow: Int,
        endRow: Int,
        startColumn: Int,
        endColumn: Int,
        destination: Array<DoubleArray?>?
    ) {
        MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn)
        val rowsCount = endRow + 1 - startRow
        val columnsCount = endColumn + 1 - startColumn
        if (destination!!.size < rowsCount || destination[0]!!.size < columnsCount) {
            throw MatrixDimensionMismatchException(
                destination.size, destination[0]!!.size,
                rowsCount, columnsCount
            )
        }
        for (i in 1 until rowsCount) {
            if (destination[i]!!.size < columnsCount) {
                throw MatrixDimensionMismatchException(
                    destination.size, destination[i]!!.size,
                    rowsCount, columnsCount
                )
            }
        }
        walkInOptimizedOrder(object : DefaultRealMatrixPreservingVisitor() {
            /** Initial row index.  */
            private var startRow = 0

            /** Initial column index.  */
            private var startColumn = 0

            /** {@inheritDoc}  */
            override fun start(
                rows: Int, columns: Int,
                startRow: Int, endRow: Int,
                startColumn: Int, endColumn: Int
            ) {
                this.startRow = startRow
                this.startColumn = startColumn
            }

            /** {@inheritDoc}  */
            override fun visit(row: Int, column: Int, value: Double) {
                destination[row - startRow]!![column - startColumn] = value
            }
        }, startRow, endRow, startColumn, endColumn)
    }

    /**
     * {@inheritDoc}
     */
    @Throws(
        OutOfRangeException::class,
        NullArgumentException::class,
        NoDataException::class,
        MatrixDimensionMismatchException::class
    )
    override fun copySubMatrix(
        selectedRows: IntArray?,
        selectedColumns: IntArray?,
        destination: Array<DoubleArray?>?
    ) {
        MatrixUtils.checkSubMatrixIndex(this, selectedRows, selectedColumns)
        val nCols = selectedColumns!!.size
        if (destination!!.size < selectedRows!!.size ||
            destination[0]!!.size < nCols
        ) {
            throw MatrixDimensionMismatchException(
                destination.size, destination[0]!!.size,
                selectedRows.size, selectedColumns.size
            )
        }
        for (i in selectedRows.indices) {
            val destinationI = destination[i]
            if (destinationI!!.size < nCols) {
                throw MatrixDimensionMismatchException(
                    destination.size, destinationI.size,
                    selectedRows.size, selectedColumns.size
                )
            }
            for (j in selectedColumns.indices) {
                destinationI[j] = getEntry(selectedRows[i], selectedColumns[j])
            }
        }
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
    override fun setSubMatrix(subMatrix: Array<DoubleArray?>?, row: Int, column: Int) {
        checkNotNull(subMatrix)
        val nRows = subMatrix!!.size
        if (nRows == 0) {
            throw NoDataException(LocalizedFormats.AT_LEAST_ONE_ROW)
        }
        val nCols: Int = subMatrix[0]!!.size
        if (nCols == 0) {
            throw NoDataException(LocalizedFormats.AT_LEAST_ONE_COLUMN)
        }
        for (r in 1 until nRows) {
            if (subMatrix[r]!!.size != nCols) {
                throw DimensionMismatchException(nCols, subMatrix[r]!!.size)
            }
        }
        MatrixUtils.checkRowIndex(this, row)
        MatrixUtils.checkColumnIndex(this, column)
        MatrixUtils.checkRowIndex(this, nRows + row - 1)
        MatrixUtils.checkColumnIndex(this, nCols + column - 1)
        for (i in 0 until nRows) {
            for (j in 0 until nCols) {
                setEntry(row + i, column + j, subMatrix[i]!![j])
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class)
    override fun getRowMatrix(row: Int): RealMatrix? {
        MatrixUtils.checkRowIndex(this, row)
        val nCols = getColumnDimension()
        val out = createMatrix(1, nCols)
        for (i in 0 until nCols) {
            out!!.setEntry(0, i, getEntry(row, i))
        }
        return out
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class, MatrixDimensionMismatchException::class)
    override fun setRowMatrix(row: Int, matrix: RealMatrix?) {
        MatrixUtils.checkRowIndex(this, row)
        val nCols = getColumnDimension()
        if (matrix!!.getRowDimension() != 1 ||
            matrix.getColumnDimension() != nCols
        ) {
            throw MatrixDimensionMismatchException(
                matrix.getRowDimension(),
                matrix.getColumnDimension(),
                1, nCols
            )
        }
        for (i in 0 until nCols) {
            setEntry(row, i, matrix.getEntry(0, i))
        }
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class)
    override fun getColumnMatrix(column: Int): RealMatrix? {
        MatrixUtils.checkColumnIndex(this, column)
        val nRows = getRowDimension()
        val out = createMatrix(nRows, 1)
        for (i in 0 until nRows) {
            out!!.setEntry(i, 0, getEntry(i, column))
        }
        return out
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class, MatrixDimensionMismatchException::class)
    override fun setColumnMatrix(column: Int, matrix: RealMatrix?) {
        MatrixUtils.checkColumnIndex(this, column)
        val nRows = getRowDimension()
        if (matrix!!.getRowDimension() != nRows ||
            matrix.getColumnDimension() != 1
        ) {
            throw MatrixDimensionMismatchException(
                matrix.getRowDimension(),
                matrix.getColumnDimension(),
                nRows, 1
            )
        }
        for (i in 0 until nRows) {
            setEntry(i, column, matrix.getEntry(i, 0))
        }
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class)
    override fun getRowVector(row: Int): RealVector? {
        return ArrayRealVector(getRow(row), false)
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class, MatrixDimensionMismatchException::class)
    override fun setRowVector(row: Int, vector: RealVector?) {
        MatrixUtils.checkRowIndex(this, row)
        val nCols = getColumnDimension()
        if (vector!!.getDimension() != nCols) {
            throw MatrixDimensionMismatchException(
                1, vector.getDimension(),
                1, nCols
            )
        }
        for (i in 0 until nCols) {
            setEntry(row, i, vector.getEntry(i))
        }
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class)
    override fun getColumnVector(column: Int): RealVector? {
        return ArrayRealVector(getColumn(column), false)
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class, MatrixDimensionMismatchException::class)
    override fun setColumnVector(column: Int, vector: RealVector?) {
        MatrixUtils.checkColumnIndex(this, column)
        val nRows = getRowDimension()
        if (vector!!.getDimension() != nRows) {
            throw MatrixDimensionMismatchException(
                vector.getDimension(), 1,
                nRows, 1
            )
        }
        for (i in 0 until nRows) {
            setEntry(i, column, vector.getEntry(i))
        }
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class)
    override fun getRow(row: Int): DoubleArray? {
        MatrixUtils.checkRowIndex(this, row)
        val nCols = getColumnDimension()
        val out = DoubleArray(nCols)
        for (i in 0 until nCols) {
            out[i] = getEntry(row, i)
        }
        return out
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class, MatrixDimensionMismatchException::class)
    override fun setRow(row: Int, array: DoubleArray?) {
        MatrixUtils.checkRowIndex(this, row)
        val nCols = getColumnDimension()
        if (array!!.size != nCols) {
            throw MatrixDimensionMismatchException(1, array.size, 1, nCols)
        }
        for (i in 0 until nCols) {
            setEntry(row, i, array[i])
        }
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class)
    override fun getColumn(column: Int): DoubleArray? {
        MatrixUtils.checkColumnIndex(this, column)
        val nRows = getRowDimension()
        val out = DoubleArray(nRows)
        for (i in 0 until nRows) {
            out[i] = getEntry(i, column)
        }
        return out
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class, MatrixDimensionMismatchException::class)
    override fun setColumn(column: Int, array: DoubleArray?) {
        MatrixUtils.checkColumnIndex(this, column)
        val nRows = getRowDimension()
        if (array!!.size != nRows) {
            throw MatrixDimensionMismatchException(array.size, 1, nRows, 1)
        }
        for (i in 0 until nRows) {
            setEntry(i, column, array[i])
        }
    }

    /**
     * {@inheritDoc}
     */
    @Throws(MathIllegalNumberException::class)
    override fun addToEntry(row: Int, column: Int, increment: Double) {
        MatrixUtils.checkMatrixIndex(this, row, column)
        setEntry(row, column, getEntry(row, column) + increment)
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class)
    override fun multiplyEntry(row: Int, column: Int, factor: Double) {
        MatrixUtils.checkMatrixIndex(this, row, column)
        setEntry(row, column, getEntry(row, column) * factor)
    }

    /**
     * {@inheritDoc}
     */
    override fun transpose(): RealMatrix? {
        val nRows = getRowDimension()
        val nCols = getColumnDimension()
        val out = createMatrix(nCols, nRows)
        walkInOptimizedOrder(object : DefaultRealMatrixPreservingVisitor() {
            /** {@inheritDoc}  */
            override fun visit(row: Int, column: Int, value: Double) {
                out!!.setEntry(column, row, value)
            }
        })
        return out
    }

    /**
     * {@inheritDoc}
     */
    override val isSquare: Boolean
        get() = getColumnDimension() == getRowDimension()

    /**
     * Returns the number of rows of this matrix.
     *
     * @return the number of rows.
     */
    abstract override fun getRowDimension(): Int

    /**
     * Returns the number of columns of this matrix.
     *
     * @return the number of columns.
     */
    abstract override fun getColumnDimension(): Int

    /**
     * {@inheritDoc}
     */
    @Throws(NonSquareMatrixException::class)
    override fun getTrace(): Double {
        val nRows = getRowDimension()
        val nCols = getColumnDimension()
        if (nRows != nCols) {
            throw NonSquareMatrixException(nRows, nCols)
        }
        var trace = 0.0
        for (i in 0 until nRows) {
            trace += getEntry(i, i)
        }
        return trace
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
            var sum = 0.0
            for (i in 0 until nCols) {
                sum += getEntry(row, i) * v[i]
            }
            out[row] = sum
        }
        return out
    }

    /**
     * {@inheritDoc}
     */
    @Throws(DimensionMismatchException::class)
    override fun operate(v: RealVector?): RealVector? {
        return try {
            ArrayRealVector(operate((v as ArrayRealVector).dataRef), false)
        } catch (cce: ClassCastException) {
            val nRows = getRowDimension()
            val nCols = getColumnDimension()
            if (v!!.getDimension() != nCols) {
                throw DimensionMismatchException(v.getDimension(), nCols)
            }
            val out = DoubleArray(nRows)
            var row = 0
            while (row < nRows) {
                var sum = 0.0
                var i = 0
                while (i < nCols) {
                    sum += getEntry(row, i) * v.getEntry(i)
                    ++i
                }
                out[row] = sum
                ++row
            }
            ArrayRealVector(out, false)
        }
    }

    /**
     * {@inheritDoc}
     */
    @Throws(DimensionMismatchException::class)
    open fun preMultiply(v: DoubleArray): DoubleArray? {
        val nRows = getRowDimension()
        val nCols = getColumnDimension()
        if (v.size != nRows) {
            throw DimensionMismatchException(v.size, nRows)
        }
        val out = DoubleArray(nCols)
        for (col in 0 until nCols) {
            var sum = 0.0
            for (i in 0 until nRows) {
                sum += getEntry(i, col) * v[i]
            }
            out[col] = sum
        }
        return out
    }

    /**
     * {@inheritDoc}
     */
    @Throws(DimensionMismatchException::class)
    open fun preMultiply(v: RealVector): RealVector? {
        return try {
            ArrayRealVector(preMultiply((v as ArrayRealVector).dataRef), false)
        } catch (cce: ClassCastException) {
            val nRows = getRowDimension()
            val nCols = getColumnDimension()
            if (v.getDimension() != nRows) {
                throw DimensionMismatchException(v.getDimension(), nRows)
            }
            val out = DoubleArray(nCols)
            var col = 0
            while (col < nCols) {
                var sum = 0.0
                var i = 0
                while (i < nRows) {
                    sum += getEntry(i, col) * v.getEntry(i)
                    ++i
                }
                out[col] = sum
                ++col
            }
            ArrayRealVector(out, false)
        }
    }

    /**
     * {@inheritDoc}
     */
    open fun walkInRowOrder(visitor: RealMatrixChangingVisitor): Double {
        val rows = getRowDimension()
        val columns = getColumnDimension()
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1)
        for (row in 0 until rows) {
            for (column in 0 until columns) {
                val oldValue = getEntry(row, column)
                val newValue = visitor.visit(row, column, oldValue)
                setEntry(row, column, newValue)
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
        for (row in 0 until rows) {
            for (column in 0 until columns) {
                visitor.visit(row, column, getEntry(row, column))
            }
        }
        return visitor.end()
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class, NumberIsTooSmallException::class)
    open fun walkInRowOrder(
        visitor: RealMatrixChangingVisitor,
        startRow: Int, endRow: Int,
        startColumn: Int, endColumn: Int
    ): Double {
        MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn)
        visitor.start(
            getRowDimension(), getColumnDimension(),
            startRow, endRow, startColumn, endColumn
        )
        for (row in startRow..endRow) {
            for (column in startColumn..endColumn) {
                val oldValue = getEntry(row, column)
                val newValue = visitor.visit(row, column, oldValue)
                setEntry(row, column, newValue)
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
        for (row in startRow..endRow) {
            for (column in startColumn..endColumn) {
                visitor.visit(row, column, getEntry(row, column))
            }
        }
        return visitor.end()
    }

    /**
     * {@inheritDoc}
     */
    open fun walkInColumnOrder(visitor: RealMatrixChangingVisitor): Double {
        val rows = getRowDimension()
        val columns = getColumnDimension()
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1)
        for (column in 0 until columns) {
            for (row in 0 until rows) {
                val oldValue = getEntry(row, column)
                val newValue = visitor.visit(row, column, oldValue)
                setEntry(row, column, newValue)
            }
        }
        return visitor.end()
    }

    /**
     * {@inheritDoc}
     */
    open fun walkInColumnOrder(visitor: RealMatrixPreservingVisitor): Double {
        val rows = getRowDimension()
        val columns = getColumnDimension()
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1)
        for (column in 0 until columns) {
            for (row in 0 until rows) {
                visitor.visit(row, column, getEntry(row, column))
            }
        }
        return visitor.end()
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class, NumberIsTooSmallException::class)
    open fun walkInColumnOrder(
        visitor: RealMatrixChangingVisitor,
        startRow: Int, endRow: Int,
        startColumn: Int, endColumn: Int
    ): Double {
        MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn)
        visitor.start(
            getRowDimension(), getColumnDimension(),
            startRow, endRow, startColumn, endColumn
        )
        for (column in startColumn..endColumn) {
            for (row in startRow..endRow) {
                val oldValue = getEntry(row, column)
                val newValue = visitor.visit(row, column, oldValue)
                setEntry(row, column, newValue)
            }
        }
        return visitor.end()
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class, NumberIsTooSmallException::class)
    open fun walkInColumnOrder(
        visitor: RealMatrixPreservingVisitor,
        startRow: Int, endRow: Int,
        startColumn: Int, endColumn: Int
    ): Double {
        MatrixUtils.checkSubMatrixIndex(this, startRow, endRow, startColumn, endColumn)
        visitor.start(
            getRowDimension(), getColumnDimension(),
            startRow, endRow, startColumn, endColumn
        )
        for (column in startColumn..endColumn) {
            for (row in startRow..endRow) {
                visitor.visit(row, column, getEntry(row, column))
            }
        }
        return visitor.end()
    }

    /**
     * {@inheritDoc}
     */
    open fun walkInOptimizedOrder(visitor: RealMatrixChangingVisitor): Double {
        return walkInRowOrder(visitor)
    }

    /**
     * {@inheritDoc}
     */
    override fun walkInOptimizedOrder(visitor: RealMatrixPreservingVisitor?): Double {
        return walkInRowOrder(visitor)
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class, NumberIsTooSmallException::class)
    open fun walkInOptimizedOrder(
        visitor: RealMatrixChangingVisitor,
        startRow: Int, endRow: Int,
        startColumn: Int,
        endColumn: Int
    ): Double {
        return walkInRowOrder(visitor, startRow, endRow, startColumn, endColumn)
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class, NumberIsTooSmallException::class)
    override fun walkInOptimizedOrder(
        visitor: RealMatrixPreservingVisitor?,
        startRow: Int,
        endRow: Int,
        startColumn: Int,
        endColumn: Int
    ): Double {
        return walkInRowOrder(visitor, startRow, endRow, startColumn, endColumn)
    }

    /**
     * Get a string representation for this matrix.
     *
     * @return a string representation for this matrix
     */
    override fun toString(): String {
        val res = StringBuilder()
        val className = this::class.simpleName
        val shortClassName = className?.substring(className.lastIndexOf('.') + 1)
        res.append(shortClassName)
        res.append(DEFAULT_FORMAT.format(this))
        return res.toString()
    }

    /**
     * Returns true iff `object` is a
     * `RealMatrix` instance with the same dimensions as this
     * and all corresponding matrix entries are equal.
     *
     * @param object the object to test equality against.
     * @return true if object equals this
     */
    override fun equals(`object`: Any?): Boolean {
        if (`object` === this) {
            return true
        }
        if (`object` is RealMatrix == false) {
            return false
        }
        val m = `object`
        val nRows = getRowDimension()
        val nCols = getColumnDimension()
        if (m.getColumnDimension() != nCols || m.getRowDimension() != nRows) {
            return false
        }
        for (row in 0 until nRows) {
            for (col in 0 until nCols) {
                if (getEntry(row, col) != m.getEntry(row, col)) {
                    return false
                }
            }
        }
        return true
    }

    /**
     * Computes a hashcode for the matrix.
     *
     * @return hashcode for matrix
     */
    override fun hashCode(): Int {
        var ret = 7
        val nRows = getRowDimension()
        val nCols = getColumnDimension()
        ret = ret * 31 + nRows
        ret = ret * 31 + nCols
        for (row in 0 until nRows) {
            for (col in 0 until nCols) {
                ret = ret * 31 + (11 * (row + 1) + 17 * (col + 1)) *
                        hash(getEntry(row, col))
            }
        }
        return ret
    }
    /*
     * Empty implementations of these methods are provided in order to allow for
     * the use of the @Override tag with Java 1.5.
     */
    /**
     * {@inheritDoc}
     */
    @Throws(NotStrictlyPositiveException::class)
    abstract override fun createMatrix(rowDimension: Int, columnDimension: Int): RealMatrix?

    /**
     * {@inheritDoc}
     */
    abstract override fun copy(): RealMatrix?

}