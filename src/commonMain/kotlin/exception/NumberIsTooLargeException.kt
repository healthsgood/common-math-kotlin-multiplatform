package exception

import exception.util.Localizable
import exception.util.LocalizedFormats

/**
 * Exception to be thrown when a number is too large.
 *
 * @since 2.2
 */
class NumberIsTooLargeException
/**
 * Construct the exception with a specific context.
 *
 * @param specific       Specific context pattern.
 * @param wrong          Value that is larger than the maximum.
 * @param max            Maximum.
 * @param boundIsAllowed if true the maximum is included in the allowed range.
 */(
    specific: Localizable,
    wrong: Number?,
    /**
     * Higher bound.
     */
    val max: Number?,
    /**
     * Whether the maximum is included in the allowed range.
     */
    val boundIsAllowed: Boolean
) : MathIllegalNumberException(specific, wrong, max) {
    /**
     * @return the maximum.
     */
    /**
     * @return `true` if the maximum is included in the allowed range.
     */

    /**
     * Construct the exception.
     *
     * @param wrong          Value that is larger than the maximum.
     * @param max            Maximum.
     * @param boundIsAllowed if true the maximum is included in the allowed range.
     */
    constructor(
        wrong: Number?,
        max: Number?,
        boundIsAllowed: Boolean
    ) : this(
        if (boundIsAllowed) LocalizedFormats.NUMBER_TOO_LARGE else LocalizedFormats.NUMBER_TOO_LARGE_BOUND_EXCLUDED,
        wrong, max, boundIsAllowed
    )

}