package exception

import exception.util.Localizable

/**
 * Exception to be thrown when the argument is not greater than 0.
 *
 * @since 2.2
 */
class NotStrictlyPositiveException : NumberIsTooSmallException {
    /**
     * Construct the exception.
     *
     * @param value Argument.
     */
    constructor(value: Number?) : super(value, INTEGER_ZERO, false)

    /**
     * Construct the exception with a specific context.
     *
     * @param specific Specific context where the error occurred.
     * @param value    Argument.
     */
    constructor(
        specific: Localizable,
        value: Number?
    ) : super(specific, value, INTEGER_ZERO, false)
}