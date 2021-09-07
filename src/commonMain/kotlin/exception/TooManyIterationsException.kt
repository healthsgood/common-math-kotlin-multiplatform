package exception

import exception.util.LocalizedFormats

/**
 * Exception to be thrown when the maximal number of iterations is exceeded.
 *
 * @since 3.1
 */
class TooManyIterationsException(max: Number?) : MaxCountExceededException(max) {

    /**
     * Construct the exception.
     *
     * @param max Maximum number of evaluations.
     */
    init {
        context.addMessage(LocalizedFormats.ITERATIONS)
    }
}