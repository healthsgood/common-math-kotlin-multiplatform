package exception

import exception.util.Localizable
import exception.util.LocalizedFormats


/**
 * Exception to be thrown when two sets of dimensions differ.
 *
 * @since 3.0
 */
open class MultiDimensionMismatchException(
    specific: Localizable,
    wrong: Array<Int>,
    expected: Array<Int>
) : MathIllegalArgumentException(specific, wrong, expected) {
    /**
     * Wrong dimensions.
     */
    private val wrong: Array<Int>

    /**
     * Correct dimensions.
     */
    private val expected: Array<Int>

    /**
     * Construct an exception from the mismatched dimensions.
     *
     * @param wrong    Wrong dimensions.
     * @param expected Expected dimensions.
     */
    constructor(
        wrong: Array<Int>,
        expected: Array<Int>
    ) : this(LocalizedFormats.DIMENSIONS_MISMATCH, wrong, expected)

    /**
     * @return an array containing the wrong dimensions.
     */
    val wrongDimensions: Array<Int>
        get() = wrong.copyOf()

    /**
     * @return an array containing the expected dimensions.
     */
    val expectedDimensions: Array<Int>
        get() = expected.copyOf()

    /**
     * @param index Dimension index.
     * @return the wrong dimension stored at `index`.
     */
    fun getWrongDimension(index: Int): Int {
        return wrong[index]
    }

    /**
     * @param index Dimension index.
     * @return the expected dimension stored at `index`.
     */
    fun getExpectedDimension(index: Int): Int {
        return expected[index]
    }


    init {
        this.wrong = wrong.copyOf()
        this.expected = expected.copyOf()
    }
}