package linear

/**
 * Interface defining very basic matrix operations.
 *
 * @author haokangkang
 * @since 2.0
 */
interface AnyMatrix {
    /**
     * Is this a square matrix?
     *
     * @return true if the matrix is square (rowDimension = columnDimension)
     */
    val isSquare: Boolean

    /**
     * Returns the number of rows in the matrix.
     *
     * @return rowDimension
     */
    fun getRowDimension(): Int

    /**
     * Returns the number of columns in the matrix.
     *
     * @return columnDimension
     */
    fun getColumnDimension(): Int
}