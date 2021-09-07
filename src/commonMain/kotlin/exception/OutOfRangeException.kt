package exception

import exception.util.Localizable
import exception.util.LocalizedFormats

/**
 * Exception to be thrown when some argument is out of range.
 *
 * @since 2.2
 */
class OutOfRangeException
/**
 * Construct an exception from the mismatched dimensions with a
 * specific context information.
 *
 * @param specific Context information.
 * @param wrong    Requested value.
 * @param lo       Lower bound.
 * @param hi       Higher bound.
 */(
    specific: Localizable,
    wrong: Number?,
    /**
     * Lower bound.
     */
    val lo: Number?,
    /**
     * Higher bound.
     */
    val hi: Number?
) : MathIllegalNumberException(specific, wrong, lo, hi) {
    /**
     * @return the lower bound.
     */
    /**
     * @return the higher bound.
     */

    /**
     * Construct an exception from the mismatched dimensions.
     *
     * @param wrong Requested value.
     * @param lo    Lower bound.
     * @param hi    Higher bound.
     */
    constructor(
        wrong: Number?,
        lo: Number?,
        hi: Number?
    ) : this(LocalizedFormats.OUT_OF_RANGE_SIMPLE, wrong, lo, hi)

}