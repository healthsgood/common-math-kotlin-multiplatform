package linear

/**
 * Formats a vector in components list format "{v0; v1; ...; vk-1}".
 *
 * The prefix and suffix "{" and "}" and the separator "; " can be replaced by
 * any user-defined strings. The number format for components can be configured.
 *
 * White space is ignored at parse time, even if it is in the prefix, suffix
 * or separator specifications. So even if the default separator does include a space
 * character that is used at format time, both input string "{1;1;1}" and
 * " { 1 ; 1 ; 1 } " will be parsed without error and the same vector will be
 * returned. In the second case, however, the parse position after parsing will be
 * just after the closing curly brace, i.e. just before the trailing space.
 *
 * @author haokangkang
 * @since 2.0
 */
class RealVectorFormat constructor(
    /**
     * Prefix.
     */
    private val prefix: String = DEFAULT_PREFIX,
    /**
     * Suffix.
     */
    private val suffix: String = DEFAULT_SUFFIX,
    /**
     * Separator.
     */
    private val separator: String = DEFAULT_SEPARATOR
) {


    /**
     * This method calls [.format].
     *
     * @param v RealVector object to format.
     * @return a formatted vector.
     */
    fun format(v: RealVector): String {
        return format(v, StringBuilder()).toString()
    }

    /**
     * Formats a [RealVector] object to produce a string.
     *
     * @param vector     the object to format.
     * @param toAppendTo where the text is to be appended
     * @return the value passed in as toAppendTo.
     */
    private fun format(
        vector: RealVector, toAppendTo: StringBuilder,
    ): StringBuilder {

        // format prefix
        toAppendTo.append(prefix)

        // format components
        for (i in 0 until vector.getDimension()) {
            if (i > 0) {
                toAppendTo.append(separator)
            }
            toAppendTo.append(vector.getEntry(i))
        }

        // format suffix
        toAppendTo.append(suffix)
        return toAppendTo
    }


    companion object {
        /**
         * The default prefix: "{".
         */
        private const val DEFAULT_PREFIX = "{"

        /**
         * The default suffix: "}".
         */
        private const val DEFAULT_SUFFIX = "}"

        /**
         * The default separator: ", ".
         */
        private const val DEFAULT_SEPARATOR = "; "

        /**
         * Returns the default real vector format
         *
         * @return the real vector format
         */

        val instance: RealVectorFormat
            get() = RealVectorFormat()
    }
}