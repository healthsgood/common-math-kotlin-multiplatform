package linear.visitor

/**
 * Interface defining a visitor for matrix entries.
 *
 * @author haokangkang
 * @see DefaultRealMatrixChangingVisitor
 *
 * @since 2.0
 */
interface RealMatrixChangingVisitor {
    /**
     * Start visiting a matrix.
     *
     * This method is called once before any entry of the matrix is visited.
     *
     * @param rows        number of rows of the matrix
     * @param columns     number of columns of the matrix
     * @param startRow    Initial row index
     * @param endRow      Final row index (inclusive)
     * @param startColumn Initial column index
     * @param endColumn   Final column index (inclusive)
     */
    fun start(
        rows: Int, columns: Int,
        startRow: Int, endRow: Int, startColumn: Int, endColumn: Int
    )

    /**
     * Visit one matrix entry.
     *
     * @param row    row index of the entry
     * @param column column index of the entry
     * @param value  current value of the entry
     * @return the new value to be set for the entry
     */
    fun visit(row: Int, column: Int, value: Double): Double

    /**
     * End visiting a matrix.
     *
     * This method is called once after all entries of the matrix have been visited.
     *
     * @return the value that the `walkInXxxOrder` must return
     */
    fun end(): Double
}