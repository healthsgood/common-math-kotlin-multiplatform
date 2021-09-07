package exception

import exception.util.Localizable
import exception.util.LocalizedFormats

/**
 * All conditions checks that fail due to a `null` argument must throw
 * this exception.
 * This class is meant to signal a precondition violation ("null is an illegal
 * argument") and so does not extend the standard `NullPointerException`.
 * Propagation of `NullPointerException` from within Commons-Math is
 * construed to be a bug.
 *
 * @since 2.2
 */
class NullArgumentException
/**
 * @param pattern   Message pattern providing the specific context of
 * the error.
 * @param arguments Values for replacing the placeholders in `pattern`.
 */
    (
    pattern: Localizable,
    vararg arguments: Any?
) : MathIllegalArgumentException(pattern, *arguments) {
    /**
     * Default constructor.
     */
    constructor() : this(LocalizedFormats.NULL_NOT_ALLOWED)
}