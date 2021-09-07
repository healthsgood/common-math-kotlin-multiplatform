package util

import exception.MaxCountExceededException
import exception.NullArgumentException
import util.Incrementor.MaxCountExceededCallback

/**
 * Utility that increments a counter until a maximum is reached, at
 * which point, the instance will by default throw a
 * [MaxCountExceededException].
 * However, the user is able to override this behaviour by defining a
 * custom [callback][MaxCountExceededCallback], in order to e.g.
 * select which exception must be thrown.
 *
 * @author haokangkang
 * @since 3.0
 */
class Incrementor constructor(
    max: Int = 0, cb: MaxCountExceededCallback? =
        object : MaxCountExceededCallback {
            /** {@inheritDoc}  */
            @Throws(MaxCountExceededException::class)
            override fun trigger(max: Int) {
                throw MaxCountExceededException(max)
            }
        }
) {
    /**
     * Function called at counter exhaustion.
     */
    private val maxCountCallback: MaxCountExceededCallback
    /**
     * Gets the upper limit of the counter.
     *
     * @return the counter upper limit.
     */
    /**
     * Sets the upper limit for the counter.
     * This does not automatically reset the current count to zero (see
     * [.resetCount]).
     *
     * @param max Upper limit of the counter.
     */
    /**
     * Upper limit for the counter.
     */
    var maximalCount: Int
    /**
     * Gets the current count.
     *
     * @return the current count.
     */
    /**
     * Current count.
     */
    var count = 0
        private set

    /**
     * Checks whether a single increment is allowed.
     *
     * @return `false` if the next call to [ incrementCount][.incrementCount] will trigger a `MaxCountExceededException`,
     * `true` otherwise.
     */
    fun canIncrement(): Boolean {
        return count < maximalCount
    }

    /**
     * Performs multiple increments.
     * See the other [incrementCount][.incrementCount] method).
     *
     * @param value Number of increments.
     * @throws MaxCountExceededException at counter exhaustion.
     */
    @Throws(MaxCountExceededException::class)
    fun incrementCount(value: Int) {
        for (i in 0 until value) {
            incrementCount()
        }
    }

    /**
     * Adds one to the current iteration count.
     * At counter exhaustion, this method will call the
     * [trigger][MaxCountExceededCallback.trigger] method of the
     * callback object passed to the
     * [constructor][.Incrementor].
     * If not explictly set, a default callback is used that will throw
     * a `MaxCountExceededException`.
     *
     * @throws MaxCountExceededException at counter exhaustion, unless a
     * custom [callback][MaxCountExceededCallback] has been set at
     * construction.
     */
    @Throws(MaxCountExceededException::class)
    fun incrementCount() {
        if (++count > maximalCount) {
            maxCountCallback.trigger(maximalCount)
        }
    }

    /**
     * Resets the counter to 0.
     */
    fun resetCount() {
        count = 0
    }

    /**
     * Defines a method to be called at counter exhaustion.
     * The [trigger][.trigger] method should usually throw an exception.
     */
    interface MaxCountExceededCallback {
        /**
         * Function called when the maximal count has been reached.
         *
         * @param maximalCount Maximal count.
         * @throws MaxCountExceededException at counter exhaustion
         */
        @Throws(MaxCountExceededException::class)
        fun trigger(maximalCount: Int)
    }
    /**
     * Defines a maximal count and a callback method to be triggered at
     * counter exhaustion.
     *
     * @param max Maximal count.
     * @param cb  Function to be called when the maximal count has been reached.
     * @throws NullArgumentException if `cb` is `null`
     */
    /**
     * Defines a maximal count.
     *
     * @param max Maximal count.
     */
    /**
     * Default constructor.
     * For the new instance to be useful, the maximal count must be set
     * by calling [setMaximalCount][.setMaximalCount].
     */
    init {
        if (cb == null) {
            throw NullArgumentException()
        }
        maximalCount = max
        maxCountCallback = cb
    }
}