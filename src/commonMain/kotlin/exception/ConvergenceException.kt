package exception

import exception.util.Localizable
import exception.util.LocalizedFormats

/**
 * Error thrown when a numerical computation can not be performed because the
 * numerical result failed to converge to a finite value.
 *
 * @since 2.2
 */
class ConvergenceException(
    pattern: Localizable,
    vararg args: Any?
) : MathIllegalStateException() {
    /**
     * Construct the exception.
     */
    constructor() : this(LocalizedFormats.CONVERGENCE_FAILED)


    /**
     * Construct the exception with a specific context and arguments.
     *
     * @param pattern Message pattern providing the specific context of
     * the error.
     * @param args    Arguments.
     */
    init {
        context.addMessage(pattern, *args)
    }
}