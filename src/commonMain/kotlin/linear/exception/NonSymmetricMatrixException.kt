package linear.exception

import exception.MathIllegalArgumentException
import exception.util.LocalizedFormats

/**
 * Exception to be thrown when a symmetric matrix is expected.
 *
 * @author haokangkang
 * @since 3.0
 */
class NonSymmetricMatrixException
/**
 * Construct an exception.
 *
 * @param row       Row index.
 * @param column    Column index.
 * @param threshold Relative symmetry threshold.
 */(
    /**
     * Row.
     */
    val row: Int,
    /**
     * Column.
     */
    val column: Int,
    /**
     * Threshold.
     */
    val threshold: Double
) : MathIllegalArgumentException(LocalizedFormats.NON_SYMMETRIC_MATRIX, row, column, threshold)