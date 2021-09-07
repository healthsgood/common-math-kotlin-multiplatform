package exception

import exception.util.Localizable
import exception.util.LocalizedFormats


/**
 * Exception to be thrown when two dimensions differ.
 *
 * @since 2.2
 */
open class DimensionMismatchException
/**
 * Construct an exception from the mismatched dimensions.
 *
 * @param specific Specific context information pattern.
 * @param wrong    Wrong dimension.
 * @param expected Expected dimension.
 */(
    specific: Localizable,
    wrong: Int,
    /**
     * Correct dimension.
     */
    val dimension: Int
) : MathIllegalNumberException(
    specific, wrong, dimension
) {
    /**
     * @return the expected dimension.
     */

    /**
     * Construct an exception from the mismatched dimensions.
     *
     * @param wrong    Wrong dimension.
     * @param expected Expected dimension.
     */
    constructor(
        wrong: Int,
        expected: Int
    ) : this(LocalizedFormats.DIMENSIONS_MISMATCH_SIMPLE, wrong, expected)

}