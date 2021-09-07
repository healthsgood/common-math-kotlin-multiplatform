package linear.exception

import exception.DimensionMismatchException
import exception.util.LocalizedFormats

/**
 * Exception to be thrown when a square matrix is expected.
 *
 * @author haokangkang
 * @since 3.0
 */
class NonSquareMatrixException
/**
 * Construct an exception from the mismatched dimensions.
 *
 * @param wrong    Row dimension.
 * @param expected Column dimension.
 */
    (
    wrong: Int,
    expected: Int
) : DimensionMismatchException(LocalizedFormats.NON_SQUARE_MATRIX, wrong, expected) {
    companion object {
        /**
         * Serializable version Id.
         */
        private const val serialVersionUID = -660069396594485772L
    }
}