package linear

import exception.*
import exception.util.LocalizedFormats
import linear.exception.MatrixDimensionMismatchException
import linear.exception.NonSquareMatrixException
import linear.exception.NonSymmetricMatrixException
import linear.exception.SingularMatrixException
import util.FastMath.abs
import util.FastMath.max
import util.MathUtils.checkNotNull
import util.Precision

/**
 * A collection of static methods that operate on or return matrices.
 * @author haokangkang
 */
object MatrixUtils {
    /**
     * Returns a [RealMatrix] with specified dimensions.
     *
     * The type of matrix returned depends on the dimension. Below
     * 2<sup>12</sup> elements (i.e. 4096 elements or 6464 for a
     * square matrix) which can be stored in a 32kB array, a [ ] instance is built. Above this threshold a [ ] instance is built.
     *
     * The matrix elements are all set to 0.0.
     *
     * @param rows    number of rows of the matrix
     * @param columns number of columns of the matrix
     * @return RealMatrix with specified dimensions
     * @see .createRealMatrix
     */

    fun createRealMatrix(rows: Int, columns: Int): RealMatrix {
        return if (rows * columns <= 4096) Array2DRowRealMatrix(rows, columns) else BlockRealMatrix(rows, columns)
    }

    /**
     * Returns a [RealMatrix] whose entries are the the values in the
     * the input array.
     *
     * The type of matrix returned depends on the dimension. Below
     * 2<sup>12</sup> elements (i.e. 4096 elements or 6464 for a
     * square matrix) which can be stored in a 32kB array, a [ ] instance is built. Above this threshold a [ ] instance is built.
     *
     * The input array is copied, not referenced.
     *
     * @param data input array
     * @return RealMatrix containing the values of the array
     * @throws exception.DimensionMismatchException if `data` is not rectangular (not all rows have the same length).
     * @throws NoDataException                                               if a row or column is empty.
     * @throws NullArgumentException                                         if either `data` or `data[0]`
     * is `null`.
     * @throws DimensionMismatchException                                    if `data` is not rectangular.
     * @see .createRealMatrix
     */

    @Throws(NullArgumentException::class, DimensionMismatchException::class, NoDataException::class)
    fun createRealMatrix(data: Array<DoubleArray?>?): RealMatrix {
        if (data == null ||
            data[0] == null
        ) {
            throw NullArgumentException()
        }
        return if (data.size * data[0]!!.size <= 4096) Array2DRowRealMatrix(data) else BlockRealMatrix(data as Array<DoubleArray>)
    }

    /**
     * Returns `dimension x dimension` identity matrix.
     *
     * @param dimension dimension of identity matrix to generate
     * @return identity matrix
     * @throws IllegalArgumentException if dimension is not positive
     * @since 1.1
     */

    fun createRealIdentityMatrix(dimension: Int): RealMatrix {
        val m = createRealMatrix(dimension, dimension)
        for (i in 0 until dimension) {
            m.setEntry(i, i, 1.0)
        }
        return m
    }

    /**
     * Returns a diagonal matrix with specified elements.
     *
     * @param diagonal diagonal elements of the matrix (the array elements
     * will be copied)
     * @return diagonal matrix
     * @since 2.0
     */
    fun createRealDiagonalMatrix(diagonal: DoubleArray): RealMatrix {
        val m = createRealMatrix(diagonal.size, diagonal.size)
        for (i in diagonal.indices) {
            m.setEntry(i, i, diagonal[i])
        }
        return m
    }

    /**
     * Creates a [RealVector] using the data from the input array.
     *
     * @param data the input data
     * @return a data.length RealVector
     * @throws NoDataException       if `data` is empty.
     * @throws NullArgumentException if `data` is `null`.
     */
    @Throws(NoDataException::class, NullArgumentException::class)
    fun createRealVector(data: DoubleArray?): RealVector {
        if (data == null) {
            throw NullArgumentException()
        }
        return ArrayRealVector(data, true)
    }

    /**
     * Create a row [RealMatrix] using the data from the input
     * array.
     *
     * @param rowData the input row data
     * @return a 1 x rowData.length RealMatrix
     * @throws NoDataException       if `rowData` is empty.
     * @throws NullArgumentException if `rowData` is `null`.
     */
    @Throws(NoDataException::class, NullArgumentException::class)
    fun createRowRealMatrix(rowData: DoubleArray?): RealMatrix {
        if (rowData == null) {
            throw NullArgumentException()
        }
        val nCols = rowData.size
        val m = createRealMatrix(1, nCols)
        for (i in 0 until nCols) {
            m.setEntry(0, i, rowData[i])
        }
        return m
    }

    /**
     * Creates a column [RealMatrix] using the data from the input
     * array.
     *
     * @param columnData the input column data
     * @return a columnData x 1 RealMatrix
     * @throws NoDataException       if `columnData` is empty.
     * @throws NullArgumentException if `columnData` is `null`.
     */
    @Throws(NoDataException::class, NullArgumentException::class)
    fun createColumnRealMatrix(columnData: DoubleArray?): RealMatrix {
        if (columnData == null) {
            throw NullArgumentException()
        }
        val nRows = columnData.size
        val m = createRealMatrix(nRows, 1)
        for (i in 0 until nRows) {
            m.setEntry(i, 0, columnData[i])
        }
        return m
    }

    /**
     * Checks whether a matrix is symmetric, within a given relative tolerance.
     *
     * @param matrix            Matrix to check.
     * @param relativeTolerance Tolerance of the symmetry check.
     * @param raiseException    If `true`, an exception will be raised if
     * the matrix is not symmetric.
     * @return `true` if `matrix` is symmetric.
     * @throws NonSquareMatrixException    if the matrix is not square.
     * @throws NonSymmetricMatrixException if the matrix is not symmetric.
     */
    private fun isSymmetricInternal(
        matrix: RealMatrix,
        relativeTolerance: Double,
        raiseException: Boolean
    ): Boolean {
        val rows = matrix.getRowDimension()
        if (rows != matrix.getColumnDimension()) {
            return if (raiseException) {
                throw NonSquareMatrixException(rows, matrix.getColumnDimension())
            } else {
                false
            }
        }
        for (i in 0 until rows) {
            for (j in i + 1 until rows) {
                val mij = matrix.getEntry(i, j)
                val mji = matrix.getEntry(j, i)
                if (abs(mij - mji) >
                    max(abs(mij), abs(mji)) * relativeTolerance
                ) {
                    return if (raiseException) {
                        throw NonSymmetricMatrixException(i, j, relativeTolerance)
                    } else {
                        false
                    }
                }
            }
        }
        return true
    }

    /**
     * Checks whether a matrix is symmetric.
     *
     * @param matrix Matrix to check.
     * @param eps    Relative tolerance.
     * @throws NonSquareMatrixException    if the matrix is not square.
     * @throws NonSymmetricMatrixException if the matrix is not symmetric.
     * @since 3.1
     */
    fun checkSymmetric(
        matrix: RealMatrix,
        eps: Double
    ) {
        isSymmetricInternal(matrix, eps, true)
    }

    /**
     * Checks whether a matrix is symmetric.
     *
     * @param matrix Matrix to check.
     * @param eps    Relative tolerance.
     * @return `true` if `matrix` is symmetric.
     * @since 3.1
     */
    fun isSymmetric(
        matrix: RealMatrix,
        eps: Double
    ): Boolean {
        return isSymmetricInternal(matrix, eps, false)
    }

    /**
     * Check if matrix indices are valid.
     *
     * @param m      Matrix.
     * @param row    Row index to check.
     * @param column Column index to check.
     * @throws OutOfRangeException if `row` or `column` is not
     * a valid index.
     */
    @Throws(OutOfRangeException::class)
    fun checkMatrixIndex(
        m: AnyMatrix,
        row: Int, column: Int
    ) {
        checkRowIndex(m, row)
        checkColumnIndex(m, column)
    }

    /**
     * Check if a row index is valid.
     *
     * @param m   Matrix.
     * @param row Row index to check.
     * @throws OutOfRangeException if `row` is not a valid index.
     */
    @Throws(OutOfRangeException::class)
    fun checkRowIndex(m: AnyMatrix, row: Int) {
        if (row < 0 ||
            row >= m.getRowDimension()
        ) {
            throw OutOfRangeException(
                LocalizedFormats.ROW_INDEX,
                row, 0, m.getRowDimension() - 1
            )
        }
    }

    /**
     * Check if a column index is valid.
     *
     * @param m      Matrix.
     * @param column Column index to check.
     * @throws OutOfRangeException if `column` is not a valid index.
     */
    @Throws(OutOfRangeException::class)
    fun checkColumnIndex(m: AnyMatrix, column: Int) {
        if (column < 0 || column >= m.getColumnDimension()) {
            throw OutOfRangeException(
                LocalizedFormats.COLUMN_INDEX,
                column, 0, m.getColumnDimension() - 1
            )
        }
    }

    /**
     * Check if submatrix ranges indices are valid.
     * Rows and columns are indicated counting from 0 to `n - 1`.
     *
     * @param m           Matrix.
     * @param startRow    Initial row index.
     * @param endRow      Final row index.
     * @param startColumn Initial column index.
     * @param endColumn   Final column index.
     * @throws OutOfRangeException       if the indices are invalid.
     * @throws NumberIsTooSmallException if `endRow < startRow` or
     * `endColumn < startColumn`.
     */
    @Throws(NumberIsTooSmallException::class, OutOfRangeException::class)
    fun checkSubMatrixIndex(
        m: AnyMatrix,
        startRow: Int, endRow: Int,
        startColumn: Int, endColumn: Int
    ) {
        checkRowIndex(m, startRow)
        checkRowIndex(m, endRow)
        if (endRow < startRow) {
            throw NumberIsTooSmallException(
                LocalizedFormats.INITIAL_ROW_AFTER_FINAL_ROW,
                endRow, startRow, false
            )
        }
        checkColumnIndex(m, startColumn)
        checkColumnIndex(m, endColumn)
        if (endColumn < startColumn) {
            throw NumberIsTooSmallException(
                LocalizedFormats.INITIAL_COLUMN_AFTER_FINAL_COLUMN,
                endColumn, startColumn, false
            )
        }
    }

    /**
     * Check if submatrix ranges indices are valid.
     * Rows and columns are indicated counting from 0 to n-1.
     *
     * @param m               Matrix.
     * @param selectedRows    Array of row indices.
     * @param selectedColumns Array of column indices.
     * @throws NullArgumentException if `selectedRows` or
     * `selectedColumns` are `null`.
     * @throws NoDataException       if the row or column selections are empty (zero
     * length).
     * @throws OutOfRangeException   if row or column selections are not valid.
     */
    @Throws(NoDataException::class, NullArgumentException::class, OutOfRangeException::class)
    fun checkSubMatrixIndex(
        m: AnyMatrix,
        selectedRows: IntArray?,
        selectedColumns: IntArray?
    ) {
        if (selectedRows == null) {
            throw NullArgumentException()
        }
        if (selectedColumns == null) {
            throw NullArgumentException()
        }
        if (selectedRows.size == 0) {
            throw NoDataException(LocalizedFormats.EMPTY_SELECTED_ROW_INDEX_ARRAY)
        }
        if (selectedColumns.size == 0) {
            throw NoDataException(LocalizedFormats.EMPTY_SELECTED_COLUMN_INDEX_ARRAY)
        }
        for (row in selectedRows) {
            checkRowIndex(m, row)
        }
        for (column in selectedColumns) {
            checkColumnIndex(m, column)
        }
    }

    /**
     * Check if matrices are addition compatible.
     *
     * @param left  Left hand side matrix.
     * @param right Right hand side matrix.
     * @throws MatrixDimensionMismatchException if the matrices are not addition
     * compatible.
     */
    @Throws(MatrixDimensionMismatchException::class)
    fun checkAdditionCompatible(left: AnyMatrix, right: AnyMatrix) {
        if (left.getRowDimension() != right.getRowDimension() ||
            left.getColumnDimension() != right.getColumnDimension()
        ) {
            throw MatrixDimensionMismatchException(
                left.getRowDimension(), left.getColumnDimension(),
                right.getRowDimension(), right.getColumnDimension()
            )
        }
    }

    /**
     * Check if matrices are subtraction compatible
     *
     * @param left  Left hand side matrix.
     * @param right Right hand side matrix.
     * @throws MatrixDimensionMismatchException if the matrices are not addition
     * compatible.
     */
    @Throws(MatrixDimensionMismatchException::class)
    fun checkSubtractionCompatible(left: AnyMatrix, right: AnyMatrix) {
        if (left.getRowDimension() != right.getRowDimension() ||
            left.getColumnDimension() != right.getColumnDimension()
        ) {
            throw MatrixDimensionMismatchException(
                left.getRowDimension(), left.getColumnDimension(),
                right.getRowDimension(), right.getColumnDimension()
            )
        }
    }

    /**
     * Check if matrices are multiplication compatible
     *
     * @param left  Left hand side matrix.
     * @param right Right hand side matrix.
     * @throws DimensionMismatchException if matrices are not multiplication
     * compatible.
     */
    @Throws(DimensionMismatchException::class)
    fun checkMultiplicationCompatible(left: AnyMatrix, right: AnyMatrix) {
        if (left.getColumnDimension() != right.getRowDimension()) {
            throw DimensionMismatchException(
                left.getColumnDimension(),
                right.getRowDimension()
            )
        }
    }

    /**
     * Solve  a  system of composed of a Lower Triangular Matrix
     * [RealMatrix].
     *
     *
     * This method is called to solve systems of equations which are
     * of the lower triangular form. The matrix [RealMatrix]
     * is assumed, though not checked, to be in lower triangular form.
     * The vector [RealVector] is overwritten with the solution.
     * The matrix is checked that it is square and its dimensions match
     * the length of the vector.
     *
     *
     * @param rm RealMatrix which is lower triangular
     * @param b  RealVector this is overwritten
     * @throws DimensionMismatchException if the matrix and vector are not
     * conformable
     * @throws NonSquareMatrixException   if the matrix `rm` is not square
     * @throws MathArithmeticException    if the absolute value of one of the diagonal
     * coefficient of `rm` is lower than [Precision.SAFE_MIN]
     */
    @Throws(DimensionMismatchException::class, MathArithmeticException::class, NonSquareMatrixException::class)
    fun solveLowerTriangularSystem(rm: RealMatrix?, b: RealVector?) {
        if (rm == null || b == null || rm.getRowDimension() != b.getDimension()) {
            throw DimensionMismatchException(
                rm?.getRowDimension() ?: 0,
                b?.getDimension() ?: 0
            )
        }
        if (rm.getColumnDimension() != rm.getRowDimension()) {
            throw NonSquareMatrixException(
                rm.getRowDimension(),
                rm.getColumnDimension()
            )
        }
        val rows = rm.getRowDimension()
        for (i in 0 until rows) {
            val diag = rm.getEntry(i, i)
            if (abs(diag) < Precision.SAFE_MIN) {
                throw MathArithmeticException(LocalizedFormats.ZERO_DENOMINATOR)
            }
            val bi = b.getEntry(i) / diag
            b.setEntry(i, bi)
            for (j in i + 1 until rows) {
                b.setEntry(j, b.getEntry(j) - bi * rm.getEntry(j, i))
            }
        }
    }

    /**
     * Solver a  system composed  of an Upper Triangular Matrix
     * [RealMatrix].
     *
     *
     * This method is called to solve systems of equations which are
     * of the lower triangular form. The matrix [RealMatrix]
     * is assumed, though not checked, to be in upper triangular form.
     * The vector [RealVector] is overwritten with the solution.
     * The matrix is checked that it is square and its dimensions match
     * the length of the vector.
     *
     *
     * @param rm RealMatrix which is upper triangular
     * @param b  RealVector this is overwritten
     * @throws DimensionMismatchException if the matrix and vector are not
     * conformable
     * @throws NonSquareMatrixException   if the matrix `rm` is not
     * square
     * @throws MathArithmeticException    if the absolute value of one of the diagonal
     * coefficient of `rm` is lower than [Precision.SAFE_MIN]
     */
    @Throws(DimensionMismatchException::class, MathArithmeticException::class, NonSquareMatrixException::class)
    fun solveUpperTriangularSystem(rm: RealMatrix?, b: RealVector?) {
        if (rm == null || b == null || rm.getRowDimension() != b.getDimension()) {
            throw DimensionMismatchException(
                rm?.getRowDimension() ?: 0,
                b?.getDimension() ?: 0
            )
        }
        if (rm.getColumnDimension() != rm.getRowDimension()) {
            throw NonSquareMatrixException(
                rm.getRowDimension(),
                rm.getColumnDimension()
            )
        }
        val rows = rm.getRowDimension()
        for (i in rows - 1 downTo -1 + 1) {
            val diag = rm.getEntry(i, i)
            if (abs(diag) < Precision.SAFE_MIN) {
                throw MathArithmeticException(LocalizedFormats.ZERO_DENOMINATOR)
            }
            val bi = b.getEntry(i) / diag
            b.setEntry(i, bi)
            for (j in i - 1 downTo -1 + 1) {
                b.setEntry(j, b.getEntry(j) - bi * rm.getEntry(j, i))
            }
        }
    }
    /**
     * Computes the inverse of the given matrix.
     *
     *
     * By default, the inverse of the matrix is computed using the QR-decomposition,
     * unless a more efficient method can be determined for the input matrix.
     *
     * @param matrix    Matrix whose inverse shall be computed
     * @param threshold Singularity threshold
     * @return the inverse of `m`
     * @throws NullArgumentException    if `matrix` is `null`
     * @throws SingularMatrixException  if matrix is singular
     * @throws NonSquareMatrixException if matrix is not square
     * @since 3.3
     */
    /**
     * Computes the inverse of the given matrix.
     *
     *
     * By default, the inverse of the matrix is computed using the QR-decomposition,
     * unless a more efficient method can be determined for the input matrix.
     *
     *
     * Note: this method will use a singularity threshold of 0,
     * use [.inverse] if a different threshold is needed.
     *
     * @param matrix Matrix whose inverse shall be computed
     * @return the inverse of `matrix`
     * @throws NullArgumentException    if `matrix` is `null`
     * @throws SingularMatrixException  if m is singular
     * @throws NonSquareMatrixException if matrix is not square
     * @since 3.3
     */

    @Throws(NullArgumentException::class, SingularMatrixException::class, NonSquareMatrixException::class)
    fun inverse(matrix: RealMatrix, threshold: Double = 0.0): RealMatrix? {
        checkNotNull(matrix)
        if (!matrix.isSquare) {
            throw NonSquareMatrixException(
                matrix.getRowDimension(),
                matrix.getColumnDimension()
            )
        }
        return if (matrix is DiagonalMatrix) {
            matrix.inverse(threshold)
        } else {
            val decomposition = QRDecomposition(matrix, threshold)
            decomposition.solver!!.getInverse()
        }
    }
}