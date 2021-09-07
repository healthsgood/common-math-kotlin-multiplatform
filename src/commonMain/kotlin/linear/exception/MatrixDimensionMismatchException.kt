package linear.exception

import exception.MultiDimensionMismatchException
import exception.util.LocalizedFormats

/**
 * Exception to be thrown when either the number of rows or the number of
 * columns of a matrix do not match the expected values.
 *
 * @author haokangkang
 * @since 3.0
 */
class MatrixDimensionMismatchException
/**
 * Construct an exception from the mismatched dimensions.
 *
 * @param wrongRowDim    Wrong row dimension.
 * @param wrongColDim    Wrong column dimension.
 * @param expectedRowDim Expected row dimension.
 * @param expectedColDim Expected column dimension.
 */
    (
    wrongRowDim: Int,
    wrongColDim: Int,
    expectedRowDim: Int,
    expectedColDim: Int
) : MultiDimensionMismatchException(
    LocalizedFormats.DIMENSIONS_MISMATCH_2x2,
    arrayOf(wrongRowDim, wrongColDim),
    arrayOf(expectedRowDim, expectedColDim)
) {
    /**
     * @return the expected row dimension.
     */
    val wrongRowDimension: Int
        get() = getWrongDimension(0)

    /**
     * @return the expected row dimension.
     */
    val expectedRowDimension: Int
        get() = getExpectedDimension(0)

    /**
     * @return the wrong column dimension.
     */
    val wrongColumnDimension: Int
        get() = getWrongDimension(1)

    /**
     * @return the expected column dimension.
     */
    val expectedColumnDimension: Int
        get() = getExpectedDimension(1)
}