package util

/**
 * @author hao kang kang
 */
object FastMath {


    fun sqrt(a: Double): Double {
        return kotlin.math.sqrt(a)
    }


    fun abs(x: Long): Long {
        return kotlin.math.abs(x)
    }


    fun abs(x: Double): Double {
        return kotlin.math.abs(x)
    }


    fun min(a: Int, b: Int): Int {
        return a.coerceAtMost(b)
    }


    fun min(a: Double, b: Double): Double {
        return kotlin.math.min(a, b)
    }


    fun max(a: Int, b: Int): Int {
        return a.coerceAtLeast(b)
    }


    fun max(a: Double, b: Double): Double {
        return kotlin.math.max(a, b)
    }
}