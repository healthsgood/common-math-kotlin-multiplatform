package exception

import exception.util.ExceptionContext
import exception.util.ExceptionContextProvider
import exception.util.Localizable

/**
 * Base class for all preconditions violation exceptions.
 * In most cases, this class should not be instantiated directly: it should
 * serve as a base class to create all the exceptions that have the semantics
 * of the standard [IllegalArgumentException].
 *
 * @since 2.2
 */
open class MathIllegalArgumentException(
    pattern: Localizable,
    vararg args: Any?
) : IllegalArgumentException(), ExceptionContextProvider {
    /**
     * {@inheritDoc}
     */
    /**
     * Context.
     */
    override val context: ExceptionContext

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
     * @param pattern Message pattern explaining the cause of the error.
     * @param args    Arguments.
     */
    init {
        context = ExceptionContext()
        context.addMessage(pattern, *args)
    }
}