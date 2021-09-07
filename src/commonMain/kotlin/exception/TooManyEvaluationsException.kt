package exception

import exception.util.LocalizedFormats

/**
 * Exception to be thrown when the maximal number of evaluations is exceeded.
 *
 * @since 3.0
 */
class TooManyEvaluationsException(max: Number?) : MaxCountExceededException(max) {

    /**
     * Construct the exception.
     *
     * @param max Maximum number of evaluations.
     */
    init {
        context.addMessage(LocalizedFormats.EVALUATIONS)
    }
}