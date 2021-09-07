package linear

/**
 * Formats a `nxm` matrix in components list format
 * "{{a<sub>0</sub><sub>0</sub>,a<sub>0</sub><sub>1</sub>, ...,
 * a<sub>0</sub><sub>m-1</sub>},{a<sub>1</sub><sub>0</sub>,
 * a<sub>1</sub><sub>1</sub>, ..., a<sub>1</sub><sub>m-1</sub>},{...},{
 * a<sub>n-1</sub><sub>0</sub>, a<sub>n-1</sub><sub>1</sub>, ...,
 * a<sub>n-1</sub><sub>m-1</sub>}}".
 *
 * The prefix and suffix "{" and "}", the row prefix and suffix "{" and "}",
 * the row separator "," and the column separator "," can be replaced by any
 * user-defined strings. The number format for components can be configured.
 *
 *
 * White space is ignored at parse time, even if it is in the prefix, suffix
 * or separator specifications. So even if the default separator does include a space
 * character that is used at format time, both input string "{{1,1,1}}" and
 * " { { 1 , 1 , 1 } } " will be parsed without error and the same matrix will be
 * returned. In the second case, however, the parse position after parsing will be
 * just after the closing curly brace, i.e. just before the trailing space.
 *
 *
 * **Note:** the grouping functionality of the used "NumberFormat" is
 * disabled to prevent problems when parsing (e.g. 1,345.34 would be a valid number
 * but conflicts with the default column separator).
 *
 * @author haokangkang
 * @since 3.1
 */
class RealMatrixFormat constructor(
    /**
     * Prefix.
     */
    private val prefix: String = DEFAULT_PREFIX,
    /**
     * Suffix.
     */
    private val suffix: String = DEFAULT_SUFFIX,
    /**
     * Row prefix.
     */
    private val rowPrefix: String = DEFAULT_ROW_PREFIX,
    /**
     * Row suffix.
     */
    private val rowSuffix: String = DEFAULT_ROW_SUFFIX,
    /**
     * Row separator.
     */
    private val rowSeparator: String = DEFAULT_ROW_SEPARATOR,
    /**
     * Column separator.
     */
    private val columnSeparator: String = DEFAULT_COLUMN_SEPARATOR
) {

    /**
     * This method calls [.format].
     *
     * @param m RealMatrix object to format.
     * @return a formatted matrix.
     */
    fun format(m: RealMatrix): String {
        return format(m, StringBuilder()).toString()
    }

    /**
     * Formats a [RealMatrix] object to produce a string.
     *
     * @param matrix     the object to format.
     * @param toAppendTo where the text is to be appended
     * @return the value passed in as toAppendTo.
     */
    private fun format(
        matrix: RealMatrix, toAppendTo: StringBuilder
    ): StringBuilder {

        // format prefix
        toAppendTo.append(prefix)

        // format rows
        val rows = matrix.getRowDimension()
        for (i in 0 until rows) {
            toAppendTo.append(rowPrefix)
            for (j in 0 until matrix.getColumnDimension()) {
                if (j > 0) {
                    toAppendTo.append(columnSeparator)
                }
                toAppendTo.append(matrix.getEntry(i, j))
            }
            toAppendTo.append(rowSuffix)
            if (i < rows - 1) {
                toAppendTo.append(rowSeparator)
            }
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
         * The default row prefix: "{".
         */
        private const val DEFAULT_ROW_PREFIX = "{"

        /**
         * The default row suffix: "}".
         */
        private const val DEFAULT_ROW_SUFFIX = "}"

        /**
         * The default row separator: ",".
         */
        private const val DEFAULT_ROW_SEPARATOR = ","

        /**
         * The default column separator: ",".
         */
        private const val DEFAULT_COLUMN_SEPARATOR = ","

        /**
         * Returns the default real vector format for the given locale.
         *
         * @return the real vector format specific to the given locale.
         */
        val instance: RealMatrixFormat
            get() = RealMatrixFormat()
    }

}