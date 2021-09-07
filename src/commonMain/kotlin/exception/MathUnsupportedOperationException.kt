package exception

import exception.util.ExceptionContext
import exception.util.ExceptionContextProvider
import exception.util.Localizable
import exception.util.LocalizedFormats

/**
 * Base class for all unsupported features.
 * It is used for all the exceptions that have the semantics of the standard
 * [UnsupportedOperationException], but must also provide a localized
 * message.
 *
 * @since 2.2
 */
class MathUnsupportedOperationException(
    pattern: Localizable,
    vararg args: Any?
) : UnsupportedOperationException(), ExceptionContextProvider {
    /**
     * {@inheritDoc}
     */
    /**
     * Context.
     */
    override val context: ExceptionContext = ExceptionContext()

    /**
     * Default constructor.
     */
    constructor() : this(LocalizedFormats.UNSUPPORTED_OPERATION)

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


    /**
     * @param pattern Message pattern providing the specific context of
     * the error.
     * @param args    Arguments.
     */
    init {
        context.addMessage(pattern, *args)
    }
}