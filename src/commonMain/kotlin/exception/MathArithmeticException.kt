package exception

import exception.util.ExceptionContext
import exception.util.ExceptionContextProvider
import exception.util.Localizable
import exception.util.LocalizedFormats

/**
 * Base class for arithmetic exceptions.
 * It is used for all the exceptions that have the semantics of the standard
 * [ArithmeticException], but must also provide a localized
 * message.
 *
 * @since 3.0
 */
class MathArithmeticException : ArithmeticException, ExceptionContextProvider {
    /**
     * Context.
     */
    override val context: ExceptionContext

    /**
     * Default constructor.
     */
    constructor() {
        context = ExceptionContext()
        context.addMessage(LocalizedFormats.ARITHMETIC_EXCEPTION)
    }

    /**
     * Constructor with a specific message.
     *
     * @param pattern Message pattern providing the specific context of
     * the error.
     * @param args    Arguments.
     */
    constructor(
        pattern: Localizable,
        vararg args: Any?
    ) {
        context = ExceptionContext()
        context.addMessage(pattern, *args)
    }

    /**
     * {@inheritDoc}
     */
    override val message: String
        get() = context.message

    /**
     * {@inheritDoc}
     */
    fun getLocalizedMessage(): String {
        return context.localizedMessage
    }
}