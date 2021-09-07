package exception.util

/**
 * Class that contains the actual implementation of the functionality mandated
 * by the [ExceptionContext] interface.
 * All Commons Math exceptions delegate the interface's methods to this class.
 *
 * @since 3.0
 */
class ExceptionContext {
    /**
     * Get a reference to the exception to which the context relates.
     *
     * @return a reference to the exception to which the context relates
     */
    /**
     * Various informations that enrich the informative message.
     */
    private var msgPatterns: MutableList<Localizable> = ArrayList()

    /**
     * Various informations that enrich the informative message.
     * The arguments will replace the corresponding place-holders in
     * [.msgPatterns].
     */
    private var msgArguments: MutableList<Array<Any?>?> = ArrayList()

    /**
     * Arbitrary context information.
     */
    private var context: MutableMap<String, Any> = HashMap()

    /**
     * Adds a message.
     *
     * @param pattern   Message pattern.
     * @param arguments Values for replacing the placeholders in the message
     * pattern.
     */
    fun addMessage(
        pattern: Localizable,
        vararg arguments: Any?
    ) {
        msgPatterns.add(pattern)
        msgArguments.add(ArgUtils.flatten(arguments))
    }

    /**
     * Sets the context (key, value) pair.
     * Keys are assumed to be unique within an instance. If the same key is
     * assigned a new value, the previous one will be lost.
     *
     * @param key   Context key (not null).
     * @param value Context value.
     */
    fun setValue(key: String, value: Any) {
        context[key] = value
    }

    /**
     * Gets the value associated to the given context key.
     *
     * @param key Context key.
     * @return the context value or `null` if the key does not exist.
     */
    fun getValue(key: String): Any? {
        return context[key]
    }

    /**
     * Gets all the keys stored in the exception
     *
     * @return the set of keys.
     */
    val keys: Set<String>
        get() = context.keys

    /**
     * Gets the default message.
     *
     * @return the message.
     */
    val message: String
        get() = buildMessage(": ")

    /**
     * Gets the message in the default locale.
     *
     * @return the localized message.
     */
    val localizedMessage: String
        get() = buildMessage(": ")


    /**
     * Gets the message in a specified locale.
     *
     * @param locale    Locale in which the message should be translated.
     * @param separator Separator inserted between the message parts.
     * @return the localized message.
     */
    fun getMessage(
        separator: String
    ): String {
        return buildMessage(separator)
    }

    /**
     * Builds a message string.
     *
     * @param locale    Locale in which the message should be translated.
     * @param separator Message separator.
     * @return a localized message string.
     */
    private fun buildMessage(
        separator: String
    ): String {
        val sb = StringBuilder()
        var count = 0
        val len = msgPatterns.size
        for (i in 0 until len) {
            val pat = msgPatterns[i]
            val args = msgArguments[i]
            sb.append(pat.getLocalizedString())
            args?.forEach {
                sb.append(", $it")
            }
            if (++count < len) {
                // Add a separator if there are other messages.
                sb.append(separator)
            }
        }
        return sb.toString()
    }
}