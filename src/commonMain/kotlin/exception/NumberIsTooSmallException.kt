package exception

import exception.util.Localizable
import exception.util.LocalizedFormats

/**
 * Exception to be thrown when a number is too small.
 *
 * @since 2.2
 */
open class NumberIsTooSmallException
/**
 * Construct the exception with a specific context.
 *
 * @param specific       Specific context pattern.
 * @param wrong          Value that is smaller than the minimum.
 * @param min            Minimum.
 * @param boundIsAllowed Whether `min` is included in the allowed range.
 */(
    specific: Localizable,
    wrong: Number?,
    /**
     * Higher bound.
     */
    val min: Number?,
    /**
     * Whether the maximum is included in the allowed range.
     */
    val boundIsAllowed: Boolean
) : MathIllegalNumberException(specific, wrong, min) {
    /**
     * @return the minimum.
     */
    /**
     * @return `true` if the minimum is included in the allowed range.
     */

    /**
     * Construct the exception.
     *
     * @param wrong          Value that is smaller than the minimum.
     * @param min            Minimum.
     * @param boundIsAllowed Whether `min` is included in the allowed range.
     */
    constructor(
        wrong: Number?,
        min: Number?,
        boundIsAllowed: Boolean
    ) : this(
        if (boundIsAllowed) LocalizedFormats.NUMBER_TOO_SMALL else LocalizedFormats.NUMBER_TOO_SMALL_BOUND_EXCLUDED,
        wrong, min, boundIsAllowed
    )

}