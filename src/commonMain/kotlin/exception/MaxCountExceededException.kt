package exception

import exception.util.Localizable
import exception.util.LocalizedFormats

/**
 * Exception to be thrown when some counter maximum value is exceeded.
 *
 * @since 3.0
 */
open class MaxCountExceededException(
    specific: Localizable,
    max: Number?,
    vararg args: Any?
) : MathIllegalStateException() {
    /**
     * @return the maximum number of evaluations.
     */
    /**
     * Maximum number of evaluations.
     */
    val max: Number?

    /**
     * Construct the exception.
     *
     * @param max Maximum.
     */
    constructor(max: Number?) : this(LocalizedFormats.MAX_COUNT_EXCEEDED, max)


    /**
     * Construct the exception with a specific context.
     *
     * @param specific Specific context pattern.
     * @param max      Maximum.
     * @param args     Additional arguments.
     */
    init {
        context.addMessage(specific, max, args)
        this.max = max
    }
}