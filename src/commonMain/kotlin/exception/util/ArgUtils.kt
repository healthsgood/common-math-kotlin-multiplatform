package exception.util

/**
 * Utility class for transforming the list of arguments passed to
 * constructors of exceptions.
 */
object ArgUtils {
    /**
     * Transform a multidimensional array into a one-dimensional list.
     *
     * @param array Array (possibly multidimensional).
     * @return a list of all the `Object` instances contained in
     * `array`.
     */
    fun flatten(array: Array<out Any?>): Array<Any?> {
        val list: MutableList<Any?> = ArrayList()
        for (o in array) {
            if (o is Array<out Any?>) {
                for (oR in flatten(o)) {
                    list.add(oR)
                }
            } else {
                list.add(o)
            }
        }
        return list.toTypedArray()
    }
}