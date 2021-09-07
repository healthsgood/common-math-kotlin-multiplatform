package linear

import exception.*
import linear.exception.MatrixDimensionMismatchException
import linear.exception.NonSquareMatrixException
import linear.visitor.RealMatrixPreservingVisitor

/**
 * Interface defining a real-valued matrix with basic algebraic operations.
 *
 *
 * Matrix element indexing is 0-based -- e.g., `getEntry(0, 0)`
 * returns the element in the first row, first column of the matrix.
 *
 * @author haokangkang
 */
interface RealMatrix : AnyMatrix {
    /**
     * Create a new RealMatrix of the same type as the instance with the
     * supplied
     * row and column dimensions.
     *
     * @param rowDimension    the number of rows in the new matrix
     * @param columnDimension the number of columns in the new matrix
     * @return a new matrix of the same type as the instance
     * @throws NotStrictlyPositiveException if row or column dimension is not
     * positive.
     * @since 2.0
     */
    @Throws(NotStrictlyPositiveException::class)
    fun createMatrix(rowDimension: Int, columnDimension: Int): RealMatrix?

    /**
     * Returns a (deep) copy of this.
     *
     * @return matrix copy
     */
    fun copy(): RealMatrix?

    /**
     * Returns the sum of `this` and `m`.
     *
     * @param m matrix to be added
     * @return `this + m`
     * @throws MatrixDimensionMismatchException if `m` is not the same
     * size as `this`.
     */
    @Throws(MatrixDimensionMismatchException::class)
    fun add(m: RealMatrix?): RealMatrix?

    /**
     * Returns `this` minus `m`.
     *
     * @param m matrix to be subtracted
     * @return `this - m`
     * @throws MatrixDimensionMismatchException if `m` is not the same
     * size as `this`.
     */
    @Throws(MatrixDimensionMismatchException::class)
    fun subtract(m: RealMatrix?): RealMatrix?

    /**
     * Returns the result of adding `d` to each entry of `this`.
     *
     * @param d value to be added to each entry
     * @return `d + this`
     */
    fun scalarAdd(d: Double): RealMatrix?

    /**
     * Returns the result of multiplying each entry of `this` by
     * `d`.
     *
     * @param d value to multiply all entries by
     * @return `d * this`
     */
    fun scalarMultiply(d: Double): RealMatrix?

    /**
     * Returns the result of postmultiplying `this` by `m`.
     *
     * @param m matrix to postmultiply by
     * @return `this * m`
     * @throws DimensionMismatchException if
     * `columnDimension(this) != rowDimension(m)`
     */
    @Throws(DimensionMismatchException::class)
    fun multiply(m: RealMatrix?): RealMatrix?

    /**
     * Returns the result of premultiplying `this` by `m`.
     *
     * @param m matrix to premultiply by
     * @return `m * this`
     * @throws DimensionMismatchException if
     * `rowDimension(this) != columnDimension(m)`
     */
    @Throws(DimensionMismatchException::class)
    fun preMultiply(m: RealMatrix?): RealMatrix?

    /**
     * Returns matrix entries as a two-dimensional array.
     *
     * @return 2-dimensional array of entries
     */
    fun getData(): Array<DoubleArray?>?

    /**
     * Returns the [
 * maximum absolute row sum norm](http://mathworld.wolfram.com/MaximumAbsoluteRowSumNorm.html) of the matrix.
     *
     * @return norm
     */
    fun getNorm(): Double

    /**
     * Returns the [
 * Frobenius norm](http://mathworld.wolfram.com/FrobeniusNorm.html) of the matrix.
     *
     * @return norm
     */
    fun getFrobeniusNorm(): Double

    /**
     * Gets a submatrix. Rows and columns are indicated
     * counting from 0 to n-1.
     *
     * @param startRow    Initial row index
     * @param endRow      Final row index (inclusive)
     * @param startColumn Initial column index
     * @param endColumn   Final column index (inclusive)
     * @return The subMatrix containing the data of the
     * specified rows and columns.
     * @throws OutOfRangeException       if the indices are not valid.
     * @throws NumberIsTooSmallException if `endRow < startRow` or
     * `endColumn < startColumn`.
     */
    @Throws(OutOfRangeException::class, NumberIsTooSmallException::class)
    fun getSubMatrix(
        startRow: Int, endRow: Int, startColumn: Int,
        endColumn: Int
    ): RealMatrix?

    /**
     * Copy a submatrix. Rows and columns are indicated counting from 0 to n-1.
     *
     * @param startRow    Initial row index
     * @param endRow      Final row index (inclusive)
     * @param startColumn Initial column index
     * @param endColumn   Final column index (inclusive)
     * @param destination The arrays where the submatrix data should be copied
     * (if larger than rows/columns counts, only the upper-left part will be
     * used)
     * @throws OutOfRangeException              if the indices are not valid.
     * @throws NumberIsTooSmallException        if `endRow < startRow` or
     * `endColumn < startColumn`.
     * @throws MatrixDimensionMismatchException if the destination array is too
     * small.
     */
    @Throws(OutOfRangeException::class, NumberIsTooSmallException::class, MatrixDimensionMismatchException::class)
    fun copySubMatrix(
        startRow: Int, endRow: Int, startColumn: Int,
        endColumn: Int, destination: Array<DoubleArray?>?
    )

    /**
     * Copy a submatrix. Rows and columns are indicated counting from 0 to n-1.
     *
     * @param selectedRows    Array of row indices.
     * @param selectedColumns Array of column indices.
     * @param destination     The arrays where the submatrix data should be copied
     * (if larger than rows/columns counts, only the upper-left part will be
     * used)
     * @throws NullArgumentException            if the row or column selections are
     * `null`
     * @throws NoDataException                  if the row or column selections are empty (zero
     * length).
     * @throws OutOfRangeException              if the indices are not valid.
     * @throws MatrixDimensionMismatchException if the destination array is too
     * small.
     */
    @Throws(
        OutOfRangeException::class,
        NullArgumentException::class,
        NoDataException::class,
        MatrixDimensionMismatchException::class
    )
    fun copySubMatrix(
        selectedRows: IntArray?, selectedColumns: IntArray?,
        destination: Array<DoubleArray?>?
    )

    /**
     * Replace the submatrix starting at `row, column` using data in the
     * input `subMatrix` array. Indexes are 0-based.
     *
     *
     * Example:<br></br>
     * Starting with <pre>
     * 1  2  3  4
     * 5  6  7  8
     * 9  0  1  2
    </pre> *
     * and `subMatrix = {{3, 4} {5,6}}`, invoking
     * `setSubMatrix(subMatrix,1,1))` will result in <pre>
     * 1  2  3  4
     * 5  3  4  8
     * 9  5  6  2
    </pre> *
     *
     * @param subMatrix array containing the submatrix replacement data
     * @param row       row coordinate of the top, left element to be replaced
     * @param column    column coordinate of the top, left element to be replaced
     * @throws NoDataException            if `subMatrix` is empty.
     * @throws OutOfRangeException        if `subMatrix` does not fit into
     * this matrix from element in `(row, column)`.
     * @throws DimensionMismatchException if `subMatrix` is not rectangular
     * (not all rows have the same length) or empty.
     * @throws NullArgumentException      if `subMatrix` is `null`.
     * @since 2.0
     */
    @Throws(
        NoDataException::class,
        OutOfRangeException::class,
        DimensionMismatchException::class,
        NullArgumentException::class
    )
    fun setSubMatrix(subMatrix: Array<DoubleArray?>?, row: Int, column: Int)

    /**
     * Get the entries at the given row index as a row matrix.  Row indices start
     * at 0.
     *
     * @param row Row to be fetched.
     * @return row Matrix.
     * @throws OutOfRangeException if the specified row index is invalid.
     */
    @Throws(OutOfRangeException::class)
    fun getRowMatrix(row: Int): RealMatrix?

    /**
     * Sets the specified `row` of `this` matrix to the entries of
     * the specified row `matrix`. Row indices start at 0.
     *
     * @param row    Row to be set.
     * @param matrix Row matrix to be copied (must have one row and the same
     * number of columns as the instance).
     * @throws OutOfRangeException              if the specified row index is invalid.
     * @throws MatrixDimensionMismatchException if the row dimension of the
     * `matrix` is not `1`, or the column dimensions of `this`
     * and `matrix` do not match.
     */
    @Throws(OutOfRangeException::class, MatrixDimensionMismatchException::class)
    fun setRowMatrix(row: Int, matrix: RealMatrix?)

    /**
     * Get the entries at the given column index as a column matrix. Column
     * indices start at 0.
     *
     * @param column Column to be fetched.
     * @return column Matrix.
     * @throws OutOfRangeException if the specified column index is invalid.
     */
    @Throws(OutOfRangeException::class)
    fun getColumnMatrix(column: Int): RealMatrix?

    /**
     * Sets the specified `column` of `this` matrix to the entries
     * of the specified column `matrix`. Column indices start at 0.
     *
     * @param column Column to be set.
     * @param matrix Column matrix to be copied (must have one column and the
     * same number of rows as the instance).
     * @throws OutOfRangeException              if the specified column index is invalid.
     * @throws MatrixDimensionMismatchException if the column dimension of the
     * `matrix` is not `1`, or the row dimensions of `this`
     * and `matrix` do not match.
     */
    @Throws(OutOfRangeException::class, MatrixDimensionMismatchException::class)
    fun setColumnMatrix(column: Int, matrix: RealMatrix?)

    /**
     * Returns the entries in row number `row` as a vector. Row indices
     * start at 0.
     *
     * @param row Row to be fetched.
     * @return a row vector.
     * @throws OutOfRangeException if the specified row index is invalid.
     */
    @Throws(OutOfRangeException::class)
    fun getRowVector(row: Int): RealVector?

    /**
     * Sets the specified `row` of `this` matrix to the entries of
     * the specified `vector`. Row indices start at 0.
     *
     * @param row    Row to be set.
     * @param vector row vector to be copied (must have the same number of
     * column as the instance).
     * @throws OutOfRangeException              if the specified row index is invalid.
     * @throws MatrixDimensionMismatchException if the `vector` dimension
     * does not match the column dimension of `this` matrix.
     */
    @Throws(OutOfRangeException::class, MatrixDimensionMismatchException::class)
    fun setRowVector(row: Int, vector: RealVector?)

    /**
     * Get the entries at the given column index as a vector. Column indices
     * start at 0.
     *
     * @param column Column to be fetched.
     * @return a column vector.
     * @throws OutOfRangeException if the specified column index is invalid
     */
    @Throws(OutOfRangeException::class)
    fun getColumnVector(column: Int): RealVector?

    /**
     * Sets the specified `column` of `this` matrix to the entries
     * of the specified `vector`. Column indices start at 0.
     *
     * @param column Column to be set.
     * @param vector column vector to be copied (must have the same number of
     * rows as the instance).
     * @throws OutOfRangeException              if the specified column index is invalid.
     * @throws MatrixDimensionMismatchException if the `vector` dimension
     * does not match the row dimension of `this` matrix.
     */
    @Throws(OutOfRangeException::class, MatrixDimensionMismatchException::class)
    fun setColumnVector(column: Int, vector: RealVector?)

    /**
     * Get the entries at the given row index. Row indices start at 0.
     *
     * @param row Row to be fetched.
     * @return the array of entries in the row.
     * @throws OutOfRangeException if the specified row index is not valid.
     */
    @Throws(OutOfRangeException::class)
    fun getRow(row: Int): DoubleArray?

    /**
     * Sets the specified `row` of `this` matrix to the entries
     * of the specified `array`. Row indices start at 0.
     *
     * @param row   Row to be set.
     * @param array Row matrix to be copied (must have the same number of
     * columns as the instance)
     * @throws OutOfRangeException              if the specified row index is invalid.
     * @throws MatrixDimensionMismatchException if the `array` length does
     * not match the column dimension of `this` matrix.
     */
    @Throws(OutOfRangeException::class, MatrixDimensionMismatchException::class)
    fun setRow(row: Int, array: DoubleArray?)

    /**
     * Get the entries at the given column index as an array. Column indices
     * start at 0.
     *
     * @param column Column to be fetched.
     * @return the array of entries in the column.
     * @throws OutOfRangeException if the specified column index is not valid.
     */
    @Throws(OutOfRangeException::class)
    fun getColumn(column: Int): DoubleArray?

    /**
     * Sets the specified `column` of `this` matrix to the entries
     * of the specified `array`. Column indices start at 0.
     *
     * @param column Column to be set.
     * @param array  Column array to be copied (must have the same number of
     * rows as the instance).
     * @throws OutOfRangeException              if the specified column index is invalid.
     * @throws MatrixDimensionMismatchException if the `array` length does
     * not match the row dimension of `this` matrix.
     */
    @Throws(OutOfRangeException::class, MatrixDimensionMismatchException::class)
    fun setColumn(column: Int, array: DoubleArray?)

    /**
     * Get the entry in the specified row and column. Row and column indices
     * start at 0.
     *
     * @param row    Row index of entry to be fetched.
     * @param column Column index of entry to be fetched.
     * @return the matrix entry at `(row, column)`.
     * @throws OutOfRangeException if the row or column index is not valid.
     */
    @Throws(MathIllegalNumberException::class, OutOfRangeException::class)
    fun getEntry(row: Int, column: Int): Double

    /**
     * Set the entry in the specified row and column. Row and column indices
     * start at 0.
     *
     * @param row    Row index of entry to be set.
     * @param column Column index of entry to be set.
     * @param value  the new value of the entry.
     * @throws OutOfRangeException if the row or column index is not valid
     * @since 2.0
     */
    @Throws(OutOfRangeException::class)
    fun setEntry(row: Int, column: Int, value: Double)

    /**
     * Adds (in place) the specified value to the specified entry of
     * `this` matrix. Row and column indices start at 0.
     *
     * @param row       Row index of the entry to be modified.
     * @param column    Column index of the entry to be modified.
     * @param increment value to add to the matrix entry.
     * @throws OutOfRangeException if the row or column index is not valid.
     * @since 2.0
     */
    @Throws(MathIllegalNumberException::class)
    fun addToEntry(row: Int, column: Int, increment: Double)

    /**
     * Multiplies (in place) the specified entry of `this` matrix by the
     * specified value. Row and column indices start at 0.
     *
     * @param row    Row index of the entry to be modified.
     * @param column Column index of the entry to be modified.
     * @param factor Multiplication factor for the matrix entry.
     * @throws OutOfRangeException if the row or column index is not valid.
     * @since 2.0
     */
    @Throws(OutOfRangeException::class)
    fun multiplyEntry(row: Int, column: Int, factor: Double)

    /**
     * Returns the transpose of this matrix.
     *
     * @return transpose matrix
     */
    fun transpose(): RealMatrix?

    /**
     * Returns the [
 * trace](http://mathworld.wolfram.com/MatrixTrace.html) of the matrix (the sum of the elements on the main diagonal).
     *
     * @return the trace.
     * @throws NonSquareMatrixException if the matrix is not square.
     */
    @Throws(NonSquareMatrixException::class)
    fun getTrace(): Double

    /**
     * Returns the result of multiplying this by the vector `v`.
     *
     * @param v the vector to operate on
     * @return `this * v`
     * @throws DimensionMismatchException if the length of `v` does not
     * match the column dimension of `this`.
     */
    @Throws(DimensionMismatchException::class)
    fun operate(v: DoubleArray?): DoubleArray?

    /**
     * Returns the result of multiplying this by the vector `v`.
     *
     * @param v the vector to operate on
     * @return `this * v`
     * @throws DimensionMismatchException if the dimension of `v` does not
     * match the column dimension of `this`.
     */
    @Throws(DimensionMismatchException::class)
    fun operate(v: RealVector?): RealVector?

    /**
     * Visit (but don't change) all matrix entries in row order.
     *
     * Row order starts at upper left and iterating through all elements
     * of a row from left to right before going to the leftmost element
     * of the next row.
     *
     * @param visitor visitor used to process all matrix entries
     * @return the value returned by [RealMatrixPreservingVisitor.end] at the end
     * of the walk
     * @see . walkInRowOrder
     * @see . walkInRowOrder
     * @see .walkInRowOrder
     * @see . walkInColumnOrder
     * @see . walkInColumnOrder
     * @see . walkInColumnOrder
     * @see . walkInColumnOrder
     * @see . walkInOptimizedOrder
     * @see .walkInOptimizedOrder
     * @see . walkInOptimizedOrder
     * @see .walkInOptimizedOrder
     */
    fun walkInRowOrder(visitor: RealMatrixPreservingVisitor?): Double

    /**
     * Visit (but don't change) some matrix entries in row order.
     *
     * Row order starts at upper left and iterating through all elements
     * of a row from left to right before going to the leftmost element
     * of the next row.
     *
     * @param visitor     visitor used to process all matrix entries
     * @param startRow    Initial row index
     * @param endRow      Final row index (inclusive)
     * @param startColumn Initial column index
     * @param endColumn   Final column index
     * @return the value returned by [RealMatrixPreservingVisitor.end] at the end
     * of the walk
     * @throws OutOfRangeException       if the indices are not valid.
     * @throws NumberIsTooSmallException if `endRow < startRow` or
     * `endColumn < startColumn`.
     * @see . walkInRowOrder
     * @see .walkInRowOrder
     * @see . walkInRowOrder
     * @see . walkInColumnOrder
     * @see . walkInColumnOrder
     * @see . walkInColumnOrder
     * @see . walkInColumnOrder
     * @see . walkInOptimizedOrder
     * @see .walkInOptimizedOrder
     * @see . walkInOptimizedOrder
     * @see .walkInOptimizedOrder
     */
    @Throws(OutOfRangeException::class, NumberIsTooSmallException::class)
    fun walkInRowOrder(
        visitor: RealMatrixPreservingVisitor?, startRow: Int,
        endRow: Int, startColumn: Int, endColumn: Int
    ): Double

    /**
     * Visit (but don't change) all matrix entries using the fastest possible order.
     *
     * The fastest walking order depends on the exact matrix class. It may be
     * different from traditional row or column orders.
     *
     * @param visitor visitor used to process all matrix entries
     * @return the value returned by [RealMatrixPreservingVisitor.end] at the end
     * of the walk
     * @see . walkInRowOrder
     * @see .walkInRowOrder
     * @see . walkInRowOrder
     * @see .walkInRowOrder
     * @see . walkInColumnOrder
     * @see . walkInColumnOrder
     * @see . walkInColumnOrder
     * @see . walkInColumnOrder
     * @see . walkInOptimizedOrder
     * @see . walkInOptimizedOrder
     * @see .walkInOptimizedOrder
     */
    fun walkInOptimizedOrder(visitor: RealMatrixPreservingVisitor?): Double

    /**
     * Visit (but don't change) some matrix entries using the fastest possible order.
     *
     * The fastest walking order depends on the exact matrix class. It may be
     * different from traditional row or column orders.
     *
     * @param visitor     visitor used to process all matrix entries
     * @param startRow    Initial row index
     * @param endRow      Final row index (inclusive)
     * @param startColumn Initial column index
     * @param endColumn   Final column index (inclusive)
     * @return the value returned by [RealMatrixPreservingVisitor.end] at the end
     * of the walk
     * @throws OutOfRangeException       if the indices are not valid.
     * @throws NumberIsTooSmallException if `endRow < startRow` or
     * `endColumn < startColumn`.
     * @see . walkInRowOrder
     * @see .walkInRowOrder
     * @see . walkInRowOrder
     * @see .walkInRowOrder
     * @see . walkInColumnOrder
     * @see . walkInColumnOrder
     * @see . walkInColumnOrder
     * @see . walkInColumnOrder
     * @see . walkInOptimizedOrder
     * @see .walkInOptimizedOrder
     * @see . walkInOptimizedOrder
     */
    @Throws(OutOfRangeException::class, NumberIsTooSmallException::class)
    fun walkInOptimizedOrder(
        visitor: RealMatrixPreservingVisitor?,
        startRow: Int, endRow: Int, startColumn: Int, endColumn: Int
    ): Double
}