package crossjvm

/**
 * 需要用到的数组工具类函数
 * @author HaoKangKang healthhealthgood@gmail.com 2021-08-23 16:55
 */
object ArrayUtil {
    fun fill(a: DoubleArray, `val`: Double) {
        var i = 0
        val len = a.size
        while (i < len) {
            a[i] = `val`
            i++
        }
    }

    fun fill(a: DoubleArray, fromIndex: Int, toIndex: Int, `val`: Double) {
        rangeCheck(a.size, fromIndex, toIndex)
        for (i in fromIndex until toIndex) a[i] = `val`
    }

    fun arraycopy(source: DoubleArray?, startIndex: Int, destination: DoubleArray?, destinationOffset: Int, length: Int) {
        source!!.copyInto(destination!!, destinationOffset, startIndex, startIndex + length)
    }

    private fun rangeCheck(arrayLength: Int, fromIndex: Int, toIndex: Int) {
        require(fromIndex <= toIndex) { "fromIndex($fromIndex) > toIndex($toIndex)" }
        if (fromIndex < 0) {
            throw IndexOutOfBoundsException(fromIndex.toString())
        }
        if (toIndex > arrayLength) {
            throw IndexOutOfBoundsException(toIndex.toString())
        }
    }
}