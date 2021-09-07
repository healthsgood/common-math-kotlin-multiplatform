package linear

import analysis.UnivariateFunction
import crossjvm.ArrayUtil
import exception.*
import exception.util.LocalizedFormats
import linear.RealVectorFormat.Companion.instance
import linear.visitor.RealVectorChangingVisitor
import linear.visitor.RealVectorPreservingVisitor
import util.FastMath.abs
import util.FastMath.max
import util.FastMath.sqrt
import util.MathUtils.hash

/**
 * This class implements the [RealVector] interface with a double array.
 *
 * @author haokangkang
 * @since 2.0
 */
class ArrayRealVector : RealVector {
    /**
     * Get a reference to the underlying data array.
     * This method does not make a fresh copy of the underlying data.
     *
     * @return the array of entries.
     */
    /**
     * Entries of the vector.
     */
    val dataRef: DoubleArray

    /**
     * Build a 0-length vector.
     * Zero-length vectors may be used to initialized construction of vectors
     * by data gathering. We start with zero-length and use either the [ ][.ArrayRealVector] constructor
     * or one of the `append` method ([.append],
     * [.append]) to gather data into this vector.
     */
    constructor() {
        dataRef = DoubleArray(0)
    }

    /**
     * Construct a vector of zeroes.
     *
     * @param size Size of the vector.
     */
    constructor(size: Int) {
        dataRef = DoubleArray(size)
    }

    /**
     * Construct a vector with preset values.
     *
     * @param size   Size of the vector
     * @param preset All entries will be set with this value.
     */
    constructor(size: Int, preset: Double) {
        dataRef = DoubleArray(size)
        ArrayUtil.fill(dataRef, preset)
    }

    /**
     * Construct a vector from an array, copying the input array.
     *
     * @param d Array.
     */
    constructor(d: DoubleArray) {
        dataRef = d.copyOf()
    }

    /**
     * Create a new ArrayRealVector using the input array as the underlying
     * data array.
     * If an array is built specially in order to be embedded in a
     * ArrayRealVector and not used directly, the `copyArray` may be
     * set to `false`. This will prevent the copying and improve
     * performance as no new array will be built and no data will be copied.
     *
     * @param d         Data for the new vector.
     * @param copyArray if `true`, the input array will be copied,
     * otherwise it will be referenced.
     * @throws NullArgumentException if `d` is `null`.
     * @see .ArrayRealVector
     */
    constructor(d: DoubleArray?, copyArray: Boolean) {
        if (d == null) {
            throw NullArgumentException()
        }
        dataRef = if (copyArray) d.copyOf() else d
    }

    /**
     * Construct a vector from part of a array.
     *
     * @param d    Array.
     * @param pos  Position of first entry.
     * @param size Number of entries to copy.
     * @throws NullArgumentException     if `d` is `null`.
     * @throws NumberIsTooLargeException if the size of `d` is less
     * than `pos + size`.
     */
    constructor(d: DoubleArray?, pos: Int, size: Int) {
        if (d == null) {
            throw NullArgumentException()
        }
        if (d.size < pos + size) {
            throw NumberIsTooLargeException(pos + size, d.size, true)
        }
        dataRef = DoubleArray(size)
        ArrayUtil.arraycopy(d, pos, dataRef, 0, size)
    }

    /**
     * Construct a vector from an array.
     *
     * @param d Array of `Double`s.
     */
    constructor(d: Array<Double>) {
        dataRef = DoubleArray(d.size)
        for (i in d.indices) {
            dataRef[i] = d[i]
        }
    }

    /**
     * Construct a vector from part of an array.
     *
     * @param d    Array.
     * @param pos  Position of first entry.
     * @param size Number of entries to copy.
     * @throws NullArgumentException     if `d` is `null`.
     * @throws NumberIsTooLargeException if the size of `d` is less
     * than `pos + size`.
     */
    constructor(d: Array<Double>?, pos: Int, size: Int) {
        if (d == null) {
            throw NullArgumentException()
        }
        if (d.size < pos + size) {
            throw NumberIsTooLargeException(pos + size, d.size, true)
        }
        dataRef = DoubleArray(size)
        for (i in pos until pos + size) {
            dataRef[i - pos] = d[i]
        }
    }

    /**
     * Construct a vector from another vector, using a deep copy.
     *
     * @param v vector to copy.
     * @throws NullArgumentException if `v` is `null`.
     */
    constructor(v: RealVector?) {
        if (v == null) {
            throw NullArgumentException()
        }
        dataRef = DoubleArray(v.getDimension())
        for (i in dataRef.indices) {
            dataRef[i] = v.getEntry(i)
        }
    }

    /**
     * Construct a vector from another vector.
     *
     * @param v    Vector to copy.
     * @param deep If `true` perform a deep copy, otherwise perform a
     * shallow copy.
     */

    constructor(v: ArrayRealVector, deep: Boolean = true) {
        dataRef = if (deep) v.dataRef.copyOf() else v.dataRef
    }

    /**
     * Construct a vector by appending one vector to another vector.
     *
     * @param v1 First vector (will be put in front of the new vector).
     * @param v2 Second vector (will be put at back of the new vector).
     */
    constructor(v1: ArrayRealVector, v2: ArrayRealVector) {
        dataRef = DoubleArray(v1.dataRef.size + v2.dataRef.size)
        ArrayUtil.arraycopy(v1.dataRef, 0, dataRef, 0, v1.dataRef.size)
        ArrayUtil.arraycopy(v2.dataRef, 0, dataRef, v1.dataRef.size, v2.dataRef.size)
    }

    /**
     * Construct a vector by appending one vector to another vector.
     *
     * @param v1 First vector (will be put in front of the new vector).
     * @param v2 Second vector (will be put at back of the new vector).
     */
    constructor(v1: ArrayRealVector, v2: RealVector) {
        val l1 = v1.dataRef.size
        val l2 = v2.getDimension()
        dataRef = DoubleArray(l1 + l2)
        ArrayUtil.arraycopy(v1.dataRef, 0, dataRef, 0, l1)
        for (i in 0 until l2) {
            dataRef[l1 + i] = v2.getEntry(i)
        }
    }

    /**
     * Construct a vector by appending one vector to another vector.
     *
     * @param v1 First vector (will be put in front of the new vector).
     * @param v2 Second vector (will be put at back of the new vector).
     */
    constructor(v1: RealVector, v2: ArrayRealVector) {
        val l1 = v1.getDimension()
        val l2 = v2.dataRef.size
        dataRef = DoubleArray(l1 + l2)
        for (i in 0 until l1) {
            dataRef[i] = v1.getEntry(i)
        }
        ArrayUtil.arraycopy(v2.dataRef, 0, dataRef, l1, l2)
    }

    /**
     * Construct a vector by appending one vector to another vector.
     *
     * @param v1 First vector (will be put in front of the new vector).
     * @param v2 Second vector (will be put at back of the new vector).
     */
    constructor(v1: ArrayRealVector, v2: DoubleArray) {
        val l1 = v1.getDimension()
        val l2 = v2.size
        dataRef = DoubleArray(l1 + l2)
        ArrayUtil.arraycopy(v1.dataRef, 0, dataRef, 0, l1)
        ArrayUtil.arraycopy(v2, 0, dataRef, l1, l2)
    }

    /**
     * Construct a vector by appending one vector to another vector.
     *
     * @param v1 First vector (will be put in front of the new vector).
     * @param v2 Second vector (will be put at back of the new vector).
     */
    constructor(v1: DoubleArray, v2: ArrayRealVector) {
        val l1 = v1.size
        val l2 = v2.getDimension()
        dataRef = DoubleArray(l1 + l2)
        ArrayUtil.arraycopy(v1, 0, dataRef, 0, l1)
        ArrayUtil.arraycopy(v2.dataRef, 0, dataRef, l1, l2)
    }

    /**
     * Construct a vector by appending one vector to another vector.
     *
     * @param v1 first vector (will be put in front of the new vector)
     * @param v2 second vector (will be put at back of the new vector)
     */
    constructor(v1: DoubleArray, v2: DoubleArray) {
        val l1 = v1.size
        val l2 = v2.size
        dataRef = DoubleArray(l1 + l2)
        ArrayUtil.arraycopy(v1, 0, dataRef, 0, l1)
        ArrayUtil.arraycopy(v2, 0, dataRef, l1, l2)
    }

    /**
     * {@inheritDoc}
     */
    override fun copy(): ArrayRealVector {
        return ArrayRealVector(this, true)
    }

    /**
     * {@inheritDoc}
     */
    @Throws(DimensionMismatchException::class)
    override fun add(v: RealVector): ArrayRealVector {
        return if (v is ArrayRealVector) {
            val vData = v.dataRef
            val dim = vData.size
            checkVectorDimensions(dim)
            val result = ArrayRealVector(dim)
            val resultData = result.dataRef
            for (i in 0 until dim) {
                resultData[i] = dataRef[i] + vData[i]
            }
            result
        } else {
            checkVectorDimensions(v)
            val out = dataRef.copyOf()
            val it = v.iterator()
            while (it.hasNext()) {
                val e = it.next()
                out[e.index] += e.value
            }
            ArrayRealVector(out, false)
        }
    }

    /**
     * {@inheritDoc}
     */
    @Throws(DimensionMismatchException::class)
    override fun subtract(v: RealVector): ArrayRealVector {
        return if (v is ArrayRealVector) {
            val vData = v.dataRef
            val dim = vData.size
            checkVectorDimensions(dim)
            val result = ArrayRealVector(dim)
            val resultData = result.dataRef
            for (i in 0 until dim) {
                resultData[i] = dataRef[i] - vData[i]
            }
            result
        } else {
            checkVectorDimensions(v)
            val out = dataRef.copyOf()
            val it = v.iterator()
            while (it.hasNext()) {
                val e = it.next()
                out[e.index] -= e.value
            }
            ArrayRealVector(out, false)
        }
    }

    /**
     * {@inheritDoc}
     */
    override fun map(function: UnivariateFunction): ArrayRealVector {
        return copy().mapToSelf(function)
    }

    /**
     * {@inheritDoc}
     */
    override fun mapToSelf(function: UnivariateFunction): ArrayRealVector {
        for (i in dataRef.indices) {
            dataRef[i] = function.value(dataRef[i])
        }
        return this
    }

    /**
     * {@inheritDoc}
     */
    override fun mapAddToSelf(d: Double): RealVector {
        for (i in dataRef.indices) {
            dataRef[i] += d
        }
        return this
    }

    /**
     * {@inheritDoc}
     */
    override fun mapSubtractToSelf(d: Double): RealVector {
        for (i in dataRef.indices) {
            dataRef[i] -= d
        }
        return this
    }

    /**
     * {@inheritDoc}
     */
    override fun mapMultiplyToSelf(d: Double): RealVector {
        for (i in dataRef.indices) {
            dataRef[i] *= d
        }
        return this
    }

    /**
     * {@inheritDoc}
     */
    override fun mapDivideToSelf(d: Double): RealVector {
        for (i in dataRef.indices) {
            dataRef[i] /= d
        }
        return this
    }

    /**
     * {@inheritDoc}
     */
    @Throws(DimensionMismatchException::class)
    override fun ebeMultiply(v: RealVector?): ArrayRealVector {
        return if (v is ArrayRealVector) {
            val vData = v.dataRef
            val dim = vData.size
            checkVectorDimensions(dim)
            val result = ArrayRealVector(dim)
            val resultData = result.dataRef
            for (i in 0 until dim) {
                resultData[i] = dataRef[i] * vData[i]
            }
            result
        } else {
            checkVectorDimensions(v!!)
            val out = dataRef.copyOf()
            for (i in dataRef.indices) {
                out[i] *= v.getEntry(i)
            }
            ArrayRealVector(out, false)
        }
    }

    /**
     * {@inheritDoc}
     */
    @Throws(DimensionMismatchException::class)
    override fun ebeDivide(v: RealVector?): ArrayRealVector {
        return if (v is ArrayRealVector) {
            val vData = v.dataRef
            val dim = vData.size
            checkVectorDimensions(dim)
            val result = ArrayRealVector(dim)
            val resultData = result.dataRef
            for (i in 0 until dim) {
                resultData[i] = dataRef[i] / vData[i]
            }
            result
        } else {
            checkVectorDimensions(v!!)
            val out = dataRef.copyOf()
            for (i in dataRef.indices) {
                out[i] /= v.getEntry(i)
            }
            ArrayRealVector(out, false)
        }
    }

    /**
     * {@inheritDoc}
     */
    @Throws(DimensionMismatchException::class)
    override fun dotProduct(v: RealVector): Double {
        if (v is ArrayRealVector) {
            val vData = v.dataRef
            checkVectorDimensions(vData.size)
            var dot = 0.0
            for (i in dataRef.indices) {
                dot += dataRef[i] * vData[i]
            }
            return dot
        }
        return super.dotProduct(v)
    }

    /**
     * {@inheritDoc}
     */
    override fun getNorm(): Double {
        var sum = 0.0
        for (a in dataRef) {
            sum += a * a
        }
        return sqrt(sum)
    }

    /**
     * {@inheritDoc}
     */
    override fun getL1Norm(): Double {
        var sum = 0.0
        for (a in dataRef) {
            sum += abs(a)
        }
        return sum
    }

    /**
     * {@inheritDoc}
     */
    override fun getLInfNorm(): Double {
        var max = 0.0
        for (a in dataRef) {
            max = max(max, abs(a))
        }
        return max
    }

    /**
     * {@inheritDoc}
     */
    @Throws(DimensionMismatchException::class)
    override fun getDistance(v: RealVector): Double {
        return if (v is ArrayRealVector) {
            val vData = v.dataRef
            checkVectorDimensions(vData.size)
            var sum = 0.0
            for (i in dataRef.indices) {
                val delta = dataRef[i] - vData[i]
                sum += delta * delta
            }
            sqrt(sum)
        } else {
            checkVectorDimensions(v)
            var sum = 0.0
            for (i in dataRef.indices) {
                val delta = dataRef[i] - v.getEntry(i)
                sum += delta * delta
            }
            sqrt(sum)
        }
    }

    /**
     * {@inheritDoc}
     */
    @Throws(DimensionMismatchException::class)
    override fun getL1Distance(v: RealVector): Double {
        return if (v is ArrayRealVector) {
            val vData = v.dataRef
            checkVectorDimensions(vData.size)
            var sum = 0.0
            for (i in dataRef.indices) {
                val delta = dataRef[i] - vData[i]
                sum += abs(delta)
            }
            sum
        } else {
            checkVectorDimensions(v)
            var sum = 0.0
            for (i in dataRef.indices) {
                val delta = dataRef[i] - v.getEntry(i)
                sum += abs(delta)
            }
            sum
        }
    }

    /**
     * {@inheritDoc}
     */
    @Throws(DimensionMismatchException::class)
    override fun getLInfDistance(v: RealVector): Double {
        return if (v is ArrayRealVector) {
            val vData = v.dataRef
            checkVectorDimensions(vData.size)
            var max = 0.0
            for (i in dataRef.indices) {
                val delta = dataRef[i] - vData[i]
                max = max(max, abs(delta))
            }
            max
        } else {
            checkVectorDimensions(v)
            var max = 0.0
            for (i in dataRef.indices) {
                val delta = dataRef[i] - v.getEntry(i)
                max = max(max, abs(delta))
            }
            max
        }
    }

    /**
     * {@inheritDoc}
     */
    fun outerProduct(v: RealVector): RealMatrix {
        return if (v is ArrayRealVector) {
            val vData = v.dataRef
            val m = dataRef.size
            val n = vData.size
            val out = MatrixUtils.createRealMatrix(m, n)
            for (i in 0 until m) {
                for (j in 0 until n) {
                    out.setEntry(i, j, dataRef[i] * vData[j])
                }
            }
            out
        } else {
            val m = dataRef.size
            val n = v.getDimension()
            val out = MatrixUtils.createRealMatrix(m, n)
            for (i in 0 until m) {
                for (j in 0 until n) {
                    out.setEntry(i, j, dataRef[i] * v.getEntry(j))
                }
            }
            out
        }
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class)
    override fun getEntry(index: Int): Double {
        return try {
            dataRef[index]
        } catch (e: IndexOutOfBoundsException) {
            throw OutOfRangeException(
                LocalizedFormats.INDEX, index, 0,
                getDimension() - 1
            )
        }
    }

    /**
     * {@inheritDoc}
     */
    override fun getDimension(): Int {
        return dataRef.size
    }

    /**
     * {@inheritDoc}
     */
    override fun append(v: RealVector?): RealVector {
        return try {
            ArrayRealVector(this, v as ArrayRealVector)
        } catch (cce: ClassCastException) {
            ArrayRealVector(this, v!!)
        }
    }

    /**
     * Construct a vector by appending a vector to this vector.
     *
     * @param v Vector to append to this one.
     * @return a new vector.
     */
    fun append(v: ArrayRealVector): ArrayRealVector {
        return ArrayRealVector(this, v)
    }

    /**
     * {@inheritDoc}
     */
    override fun append(`in`: Double): RealVector {
        val out = DoubleArray(dataRef.size + 1)
        ArrayUtil.arraycopy(dataRef, 0, out, 0, dataRef.size)
        out[dataRef.size] = `in`
        return ArrayRealVector(out, false)
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class, NotPositiveException::class)
    override fun getSubVector(index: Int, n: Int): RealVector {
        if (n < 0) {
            throw NotPositiveException(LocalizedFormats.NUMBER_OF_ELEMENTS_SHOULD_BE_POSITIVE, n)
        }
        val out = ArrayRealVector(n)
        try {
            ArrayUtil.arraycopy(dataRef, index, out.dataRef, 0, n)
        } catch (e: IndexOutOfBoundsException) {
            checkIndex(index)
            checkIndex(index + n - 1)
        }
        return out
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class)
    override fun setEntry(index: Int, value: Double) {
        try {
            dataRef[index] = value
        } catch (e: IndexOutOfBoundsException) {
            checkIndex(index)
        }
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class)
    override fun addToEntry(index: Int, increment: Double) {
        try {
            dataRef[index] += increment
        } catch (e: IndexOutOfBoundsException) {
            throw OutOfRangeException(
                LocalizedFormats.INDEX,
                index, 0, dataRef.size - 1
            )
        }
    }

    /**
     * {@inheritDoc}
     */
    @Throws(OutOfRangeException::class)
    override fun setSubVector(index: Int, v: RealVector?) {
        if (v is ArrayRealVector) {
            setSubVector(index, v.dataRef)
        } else {
            try {
                for (i in index until index + v!!.getDimension()) {
                    dataRef[i] = v.getEntry(i - index)
                }
            } catch (e: IndexOutOfBoundsException) {
                checkIndex(index)
                checkIndex(index + v!!.getDimension() - 1)
            }
        }
    }

    /**
     * Set a set of consecutive elements.
     *
     * @param index Index of first element to be set.
     * @param v     Vector containing the values to set.
     * @throws OutOfRangeException if the index is inconsistent with the vector
     * size.
     */
    @Throws(OutOfRangeException::class)
    fun setSubVector(index: Int, v: DoubleArray) {
        try {
            ArrayUtil.arraycopy(v, 0, dataRef, index, v.size)
        } catch (e: IndexOutOfBoundsException) {
            checkIndex(index)
            checkIndex(index + v.size - 1)
        }
    }

    /**
     * {@inheritDoc}
     */
    override fun set(value: Double) {
        ArrayUtil.fill(dataRef, value)
    }

    /**
     * {@inheritDoc}
     */
    override fun toArray(): DoubleArray {
        return dataRef.copyOf()
    }

    /**
     * {@inheritDoc}
     */
    override fun toString(): String {
        return DEFAULT_FORMAT.format(this)
    }

    /**
     * Check if instance and specified vectors have the same dimension.
     *
     * @param v Vector to compare instance with.
     * @throws DimensionMismatchException if the vectors do not
     * have the same dimension.
     */
    @Throws(DimensionMismatchException::class)
    override fun checkVectorDimensions(v: RealVector) {
        checkVectorDimensions(v.getDimension())
    }

    /**
     * Check if instance dimension is equal to some expected value.
     *
     * @param n Expected dimension.
     * @throws DimensionMismatchException if the dimension is
     * inconsistent with vector size.
     */
    @Throws(DimensionMismatchException::class)
    override fun checkVectorDimensions(n: Int) {
        if (dataRef.size != n) {
            throw DimensionMismatchException(dataRef.size, n)
        }
    }

    /**
     * Check if any coordinate of this vector is `NaN`.
     *
     * @return `true` if any coordinate of this vector is `NaN`,
     * `false` otherwise.
     */
    override fun isNaN(): Boolean {
        for (v in dataRef) {
            if (v.isNaN()) {
                return true
            }
        }
        return false
    }

    /**
     * Check whether any coordinate of this vector is infinite and none
     * are `NaN`.
     *
     * @return `true` if any coordinate of this vector is infinite and
     * none are `NaN`, `false` otherwise.
     */
    override fun isInfinite(): Boolean {
        if (isNaN()) {
            return false
        }
        for (v in dataRef) {
            if (v.isInfinite()) {
                return true
            }
        }
        return false
    }

    /**
     * {@inheritDoc}
     */
    override fun equals(other: Any?): Boolean {
        if (this === other) {
            return true
        }
        if (other !is RealVector) {
            return false
        }
        if (dataRef.size != other.getDimension()) {
            return false
        }
        if (other.isNaN()) {
            return this.isNaN()
        }
        for (i in dataRef.indices) {
            if (dataRef[i] != other.getEntry(i)) {
                return false
            }
        }
        return true
    }

    /**
     * {@inheritDoc} All `NaN` values have the same hash code.
     */
    override fun hashCode(): Int {
        return if (isNaN()) {
            9
        } else hash(dataRef)
    }

    /**
     * {@inheritDoc}
     */
    @Throws(DimensionMismatchException::class)
    override fun combine(a: Double, b: Double, y: RealVector): ArrayRealVector {
        return copy().combineToSelf(a, b, y)
    }

    /**
     * {@inheritDoc}
     */
    @Throws(DimensionMismatchException::class)
    override fun combineToSelf(a: Double, b: Double, y: RealVector): ArrayRealVector {
        if (y is ArrayRealVector) {
            val yData = y.dataRef
            checkVectorDimensions(yData.size)
            for (i in dataRef.indices) {
                dataRef[i] = a * dataRef[i] + b * yData[i]
            }
        } else {
            checkVectorDimensions(y)
            for (i in dataRef.indices) {
                dataRef[i] = a * dataRef[i] + b * y.getEntry(i)
            }
        }
        return this
    }

    /**
     * {@inheritDoc}
     */
    override fun walkInDefaultOrder(visitor: RealVectorPreservingVisitor): Double {
        visitor.start(dataRef.size, 0, dataRef.size - 1)
        for (i in dataRef.indices) {
            visitor.visit(i, dataRef[i])
        }
        return visitor.end()
    }

    /**
     * {@inheritDoc}
     */
    @Throws(NumberIsTooSmallException::class, OutOfRangeException::class)
    override fun walkInDefaultOrder(
        visitor: RealVectorPreservingVisitor,
        start: Int, end: Int
    ): Double {
        checkIndices(start, end)
        visitor.start(dataRef.size, start, end)
        for (i in start..end) {
            visitor.visit(i, dataRef[i])
        }
        return visitor.end()
    }

    /**
     * {@inheritDoc}
     *
     *
     * In this implementation, the optimized order is the default order.
     */
    override fun walkInOptimizedOrder(visitor: RealVectorPreservingVisitor): Double {
        return walkInDefaultOrder(visitor)
    }

    /**
     * {@inheritDoc}
     *
     *
     * In this implementation, the optimized order is the default order.
     */
    @Throws(NumberIsTooSmallException::class, OutOfRangeException::class)
    override fun walkInOptimizedOrder(
        visitor: RealVectorPreservingVisitor,
        start: Int,
        end: Int
    ): Double {
        return walkInDefaultOrder(visitor, start, end)
    }

    /**
     * {@inheritDoc}
     */
    override fun walkInDefaultOrder(visitor: RealVectorChangingVisitor): Double {
        visitor.start(dataRef.size, 0, dataRef.size - 1)
        for (i in dataRef.indices) {
            dataRef[i] = visitor.visit(i, dataRef[i])
        }
        return visitor.end()
    }

    /**
     * {@inheritDoc}
     */
    @Throws(NumberIsTooSmallException::class, OutOfRangeException::class)
    override fun walkInDefaultOrder(
        visitor: RealVectorChangingVisitor,
        start: Int, end: Int
    ): Double {
        checkIndices(start, end)
        visitor.start(dataRef.size, start, end)
        for (i in start..end) {
            dataRef[i] = visitor.visit(i, dataRef[i])
        }
        return visitor.end()
    }

    /**
     * {@inheritDoc}
     *
     *
     * In this implementation, the optimized order is the default order.
     */
    override fun walkInOptimizedOrder(visitor: RealVectorChangingVisitor): Double {
        return walkInDefaultOrder(visitor)
    }

    /**
     * {@inheritDoc}
     *
     *
     * In this implementation, the optimized order is the default order.
     */
    @Throws(NumberIsTooSmallException::class, OutOfRangeException::class)
    override fun walkInOptimizedOrder(
        visitor: RealVectorChangingVisitor,
        start: Int,
        end: Int
    ): Double {
        return walkInDefaultOrder(visitor, start, end)
    }

    companion object {
        /**
         * Default format.
         */
        private val DEFAULT_FORMAT = instance
    }
}