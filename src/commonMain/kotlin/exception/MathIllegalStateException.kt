package exception

import exception.util.ExceptionContext
import exception.util.ExceptionContextProvider
import exception.util.Localizable
import exception.util.LocalizedFormats

/**
 * Base class for all exceptions that signal that the process
 * throwing the exception is in a state that does not comply with
 * the set of states that it is designed to be in.
 *
 * @since 2.2
 */
open class MathIllegalStateException : IllegalStateException, ExceptionContextProvider {
    /**
     * {@inheritDoc}
     */
    /**
     * Context.
     */
    override val context: ExceptionContext

    /**
     * Simple constructor.
     *
     * @param pattern Message pattern explaining the cause of the error.
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
     * Simple constructor.
     *
     * @param cause   Root cause.
     * @param pattern Message pattern explaining the cause of the error.
     * @param args    Arguments.
     */
    constructor(
        cause: Throwable?,
        pattern: Localizable,
        vararg args: Any?
    ) : super(cause) {
        context = ExceptionContext()
        context.addMessage(pattern, *args)
    }

    /**
     * Default constructor.
     */
    constructor() : this(LocalizedFormats.ILLEGAL_STATE)

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