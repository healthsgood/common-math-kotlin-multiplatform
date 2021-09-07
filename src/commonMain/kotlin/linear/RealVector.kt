package linear

import analysis.FunctionUtils.fix2ndArgument
import analysis.UnivariateFunction
import analysis.function.Add
import analysis.function.Divide
import analysis.function.Multiply
import exception.*
import exception.util.LocalizedFormats
import linear.visitor.RealVectorChangingVisitor
import linear.visitor.RealVectorPreservingVisitor
import util.FastMath.abs
import util.FastMath.max
import util.FastMath.sqrt

/**
 * Class defining a real-valued vector with basic algebraic operations.
 *
 *
 * vector element indexing is 0-based -- e.g., `getEntry(0)`
 * returns the first element of the vector.
 *
 *
 *
 * The `code map` and `mapToSelf` methods operate
 * on vectors element-wise, i.e. they perform the same operation (adding a scalar,
 * applying a function ...) on each element in turn. The `map`
 * versions create a new vector to hold the result and do not change the instance.
 * The `mapToSelf` version uses the instance itself to store the
 * results, so the instance is changed by this method. In all cases, the result
 * vector is returned by the methods, allowing the *fluent API*
 * style, like this:
 *
 * <pre>
 * RealVector result = v.mapAddToSelf(3.4).mapToSelf(new Tan()).mapToSelf(new Power(2.3));
</pre> *
 *
 * @author haokangkang
 * @since 2.1
 */
abstract class RealVector {
    /**
     * Returns the size of the vector.
     *
     * @return the size of this vector.
     */
    abstract fun getDimension(): Int

    /**
     * Return the entry at the specified index.
     *
     * @param index Index location of entry to be fetched.
     * @return the vector entry at `index`.
     * @throws OutOfRangeException if the index is not valid.
     * @see .setEntry
     */
    @Throws(OutOfRangeException::class)
    abstract fun getEntry(index: Int): Double

    /**
     * Set a single element.
     *
     * @param index element index.
     * @param value new value for the element.
     * @throws OutOfRangeException if the index is not valid.
     * @see .getEntry
     */
    @Throws(OutOfRangeException::class)
    abstract fun setEntry(index: Int, value: Double)

    /**
     * Change an entry at the specified index.
     *
     * @param index     Index location of entry to be set.
     * @param increment Value to add to the vector entry.
     * @throws OutOfRangeException if the index is not valid.
     * @since 3.0
     */
    @Throws(OutOfRangeException::class)
    open fun addToEntry(index: Int, increment: Double) {
        setEntry(index, getEntry(index) + increment)
    }

    /**
     * Construct a new vector by appending a vector to this vector.
     *
     * @param v vector to append to this one.
     * @return a new vector.
     */
    abstract fun append(v: RealVector?): RealVector

    /**
     * Construct a new vector by appending a double to this vector.
     *
     * @param d double to append.
     * @return a new vector.
     */
    abstract fun append(d: Double): RealVector

    /**
     * Get a subvector from consecutive elements.
     *
     * @param index index of first element.
     * @param n     number of elements to be retrieved.
     * @return a vector containing n elements.
     * @throws OutOfRangeException  if the index is not valid.
     * @throws NotPositiveException if the number of elements is not positive.
     */
    @Throws(NotPositiveException::class, OutOfRangeException::class)
    abstract fun getSubVector(index: Int, n: Int): RealVector

    /**
     * Set a sequence of consecutive elements.
     *
     * @param index index of first element to be set.
     * @param v     vector containing the values to set.
     * @throws OutOfRangeException if the index is not valid.
     */
    @Throws(OutOfRangeException::class)
    abstract fun setSubVector(index: Int, v: RealVector?)

    /**
     * Check whether any coordinate of this vector is `NaN`.
     *
     * @return `true` if any coordinate of this vector is `NaN`,
     * `false` otherwise.
     */
    abstract fun isNaN(): Boolean

    /**
     * Check whether any coordinate of this vector is infinite and none are `NaN`.
     *
     * @return `true` if any coordinate of this vector is infinite and
     * none are `NaN`, `false` otherwise.
     */
    abstract fun isInfinite(): Boolean

    /**
     * Check if instance and specified vectors have the same dimension.
     *
     * @param v Vector to compare instance with.
     * @throws DimensionMismatchException if the vectors do not
     * have the same dimension.
     */
    @Throws(DimensionMismatchException::class)
    protected open fun checkVectorDimensions(v: RealVector) {
        checkVectorDimensions(v.getDimension())
    }

    /**
     * Check if instance dimension is equal to some expected value.
     *
     * @param n Expected dimension.
     * @throws DimensionMismatchException if the dimension is
     * inconsistent with the vector size.
     */
    @Throws(DimensionMismatchException::class)
    protected open fun checkVectorDimensions(n: Int) {
        val d = getDimension()
        if (d != n) {
            throw DimensionMismatchException(d, n)
        }
    }

    /**
     * Check if an index is valid.
     *
     * @param index Index to check.
     * @throws OutOfRangeException if `index` is not valid.
     */
    @Throws(OutOfRangeException::class)
    protected fun checkIndex(index: Int) {
        if (index < 0 ||
            index >= getDimension()
        ) {
            throw OutOfRangeException(
                LocalizedFormats.INDEX,
                index, 0, getDimension() - 1
            )
        }
    }

    /**
     * Checks that the indices of a subvector are valid.
     *
     * @param start the index of the first entry of the subvector
     * @param end   the index of the last entry of the subvector (inclusive)
     * @throws OutOfRangeException       if `start` of `end` are not valid
     * @throws NumberIsTooSmallException if `end < start`
     * @since 3.1
     */
    @Throws(NumberIsTooSmallException::class, OutOfRangeException::class)
    protected fun checkIndices(start: Int, end: Int) {
        val dim = getDimension()
        if (start < 0 || start >= dim) {
            throw OutOfRangeException(
                LocalizedFormats.INDEX, start, 0,
                dim - 1
            )
        }
        if (end < 0 || end >= dim) {
            throw OutOfRangeException(
                LocalizedFormats.INDEX, end, 0,
                dim - 1
            )
        }
        if (end < start) {
            // TODO Use more specific error message
            throw NumberIsTooSmallException(
                LocalizedFormats.INITIAL_ROW_AFTER_FINAL_ROW,
                end, start, false
            )
        }
    }

    /**
     * Compute the sum of this vector and `v`.
     * Returns a new vector. Does not change instance data.
     *
     * @param v Vector to be added.
     * @return `this` + `v`.
     * @throws DimensionMismatchException if `v` is not the same size as
     * `this` vector.
     */
    @Throws(DimensionMismatchException::class)
    open fun add(v: RealVector): RealVector {
        checkVectorDimensions(v)
        val result = v.copy()
        val it = iterator()
        while (it.hasNext()) {
            val e = it.next()
            val index = e.index
            result.setEntry(index, e.value + result.getEntry(index))
        }
        return result
    }

    /**
     * Subtract `v` from this vector.
     * Returns a new vector. Does not change instance data.
     *
     * @param v Vector to be subtracted.
     * @return `this` - `v`.
     * @throws DimensionMismatchException if `v` is not the same size as
     * `this` vector.
     */
    @Throws(DimensionMismatchException::class)
    open fun subtract(v: RealVector): RealVector {
        checkVectorDimensions(v)
        val result = v.mapMultiply(-1.0)
        val it = iterator()
        while (it.hasNext()) {
            val e = it.next()
            val index = e.index
            result.setEntry(index, e.value + result.getEntry(index))
        }
        return result
    }

    /**
     * Add a value to each entry.
     * Returns a new vector. Does not change instance data.
     *
     * @param d Value to be added to each entry.
     * @return `this` + `d`.
     */
    open fun mapAdd(d: Double): RealVector? {
        return copy().mapAddToSelf(d)
    }

    /**
     * Add a value to each entry.
     * The instance is changed in-place.
     *
     * @param d Value to be added to each entry.
     * @return `this`.
     */
    open fun mapAddToSelf(d: Double): RealVector? {
        return if (d != 0.0) {
            mapToSelf(fix2ndArgument(Add(), d))
        } else this
    }

    /**
     * Returns a (deep) copy of this vector.
     *
     * @return a vector copy.
     */
    abstract fun copy(): RealVector

    /**
     * Compute the dot product of this vector with `v`.
     *
     * @param v Vector with which dot product should be computed
     * @return the scalar dot product between this instance and `v`.
     * @throws DimensionMismatchException if `v` is not the same size as
     * `this` vector.
     */
    @Throws(DimensionMismatchException::class)
    open fun dotProduct(v: RealVector): Double {
        checkVectorDimensions(v)
        var d = 0.0
        val n = getDimension()
        for (i in 0 until n) {
            d += getEntry(i) * v.getEntry(i)
        }
        return d
    }

    /**
     * Computes the cosine of the angle between this vector and the
     * argument.
     *
     * @param v Vector.
     * @return the cosine of the angle between this vector and `v`.
     * @throws MathArithmeticException    if `this` or `v` is the null
     * vector
     * @throws DimensionMismatchException if the dimensions of `this` and
     * `v` do not match
     */
    @Throws(DimensionMismatchException::class, MathArithmeticException::class)
    open fun cosine(v: RealVector): Double {
        val norm = getNorm()
        val vNorm = v.getNorm()
        if (norm == 0.0 || vNorm == 0.0) {
            throw MathArithmeticException(LocalizedFormats.ZERO_NORM)
        }
        return dotProduct(v) / (norm * vNorm)
    }

    /**
     * Element-by-element division.
     *
     * @param v Vector by which instance elements must be divided.
     * @return a vector containing this[i] / v[i] for all i.
     * @throws DimensionMismatchException if `v` is not the same size as
     * `this` vector.
     */
    @Throws(DimensionMismatchException::class)
    abstract fun ebeDivide(v: RealVector?): RealVector

    /**
     * Element-by-element multiplication.
     *
     * @param v Vector by which instance elements must be multiplied
     * @return a vector containing this[i] * v[i] for all i.
     * @throws DimensionMismatchException if `v` is not the same size as
     * `this` vector.
     */
    @Throws(DimensionMismatchException::class)
    abstract fun ebeMultiply(v: RealVector?): RealVector

    /**
     * Distance between two vectors.
     *
     * This method computes the distance consistent with the
     * L<sub>2</sub> norm, i.e. the square root of the sum of
     * element differences, or Euclidean distance.
     *
     * @param v Vector to which distance is requested.
     * @return the distance between two vectors.
     * @throws DimensionMismatchException if `v` is not the same size as
     * `this` vector.
     * @see .getL1Distance
     * @see .getLInfDistance
     * @see .getNorm
     */
    @Throws(DimensionMismatchException::class)
    open fun getDistance(v: RealVector): Double {
        checkVectorDimensions(v)
        var d = 0.0
        val it = iterator()
        while (it.hasNext()) {
            val e = it.next()
            val diff = e.value - v.getEntry(e.index)
            d += diff * diff
        }
        return sqrt(d)
    }

    /**
     * Returns the L<sub>2</sub> norm of the vector.
     *
     * The L<sub>2</sub> norm is the root of the sum of
     * the squared elements.
     *
     * @return the norm.
     * @see .getL1Norm
     * @see .getLInfNorm
     * @see .getDistance
     */
    open fun getNorm(): Double {
        var sum = 0.0
        val it = iterator()
        while (it.hasNext()) {
            val e = it.next()
            val value = e.value
            sum += value * value
        }
        return sqrt(sum)
    }

    /**
     * Returns the L<sub>1</sub> norm of the vector.
     *
     * The L<sub>1</sub> norm is the sum of the absolute
     * values of the elements.
     *
     * @return the norm.
     * @see .getNorm
     * @see .getLInfNorm
     * @see .getL1Distance
     */
    open fun getL1Norm(): Double {
        var norm = 0.0
        val it = iterator()
        while (it.hasNext()) {
            val e = it.next()
            norm += abs(e.value)
        }
        return norm
    }

    /**
     * Returns the L<sub></sub> norm of the vector.
     *
     * The L<sub></sub> norm is the max of the absolute
     * values of the elements.
     *
     * @return the norm.
     * @see .getNorm
     * @see .getL1Norm
     * @see .getLInfDistance
     */
    open fun getLInfNorm(): Double {
        var norm = 0.0
        val it = iterator()
        while (it.hasNext()) {
            val e = it.next()
            norm = max(norm, abs(e.value))
        }
        return norm
    }

    /**
     * Distance between two vectors.
     *
     * This method computes the distance consistent with
     * L<sub>1</sub> norm, i.e. the sum of the absolute values of
     * the elements differences.
     *
     * @param v Vector to which distance is requested.
     * @return the distance between two vectors.
     * @throws DimensionMismatchException if `v` is not the same size as
     * `this` vector.
     */
    @Throws(DimensionMismatchException::class)
    open fun getL1Distance(v: RealVector): Double {
        checkVectorDimensions(v)
        var d = 0.0
        val it = iterator()
        while (it.hasNext()) {
            val e = it.next()
            d += abs(e.value - v.getEntry(e.index))
        }
        return d
    }

    /**
     * Distance between two vectors.
     *
     * This method computes the distance consistent with
     * L<sub></sub> norm, i.e. the max of the absolute values of
     * element differences.
     *
     * @param v Vector to which distance is requested.
     * @return the distance between two vectors.
     * @throws DimensionMismatchException if `v` is not the same size as
     * `this` vector.
     * @see .getDistance
     * @see .getL1Distance
     * @see .getLInfNorm
     */
    @Throws(DimensionMismatchException::class)
    open fun getLInfDistance(v: RealVector): Double {
        checkVectorDimensions(v)
        var d = 0.0
        val it = iterator()
        while (it.hasNext()) {
            val e = it.next()
            d = max(abs(e.value - v.getEntry(e.index)), d)
        }
        return d
    }

    /**
     * Get the index of the minimum entry.
     *
     * @return the index of the minimum entry or -1 if vector length is 0
     * or all entries are `NaN`.
     */
    fun getMinIndex(): Int {
        var minIndex = -1
        var minValue = Double.POSITIVE_INFINITY
        val iterator = iterator()
        while (iterator.hasNext()) {
            val entry = iterator.next()
            if (entry.value <= minValue) {
                minIndex = entry.index
                minValue = entry.value
            }
        }
        return minIndex
    }

    /**
     * Get the value of the minimum entry.
     *
     * @return the value of the minimum entry or `NaN` if all
     * entries are `NaN`.
     */
    fun getMinValue(): Double {
        val minIndex = getMinIndex()
        return if (minIndex < 0) Double.NaN else getEntry(minIndex)
    }

    /**
     * Get the index of the maximum entry.
     *
     * @return the index of the maximum entry or -1 if vector length is 0
     * or all entries are `NaN`
     */
    fun getMaxIndex(): Int {
        var maxIndex = -1
        var maxValue = Double.NEGATIVE_INFINITY
        val iterator = iterator()
        while (iterator.hasNext()) {
            val entry = iterator.next()
            if (entry.value >= maxValue) {
                maxIndex = entry.index
                maxValue = entry.value
            }
        }
        return maxIndex
    }

    /**
     * Get the value of the maximum entry.
     *
     * @return the value of the maximum entry or `NaN` if all
     * entries are `NaN`.
     */
    fun getMaxValue(): Double {
        val maxIndex = getMaxIndex()
        return if (maxIndex < 0) Double.NaN else getEntry(maxIndex)
    }

    /**
     * Multiply each entry by the argument. Returns a new vector.
     * Does not change instance data.
     *
     * @param d Multiplication factor.
     * @return `this` * `d`.
     */
    open fun mapMultiply(d: Double): RealVector {
        return copy().mapMultiplyToSelf(d)
    }

    /**
     * Multiply each entry.
     * The instance is changed in-place.
     *
     * @param d Multiplication factor.
     * @return `this`.
     */
    open fun mapMultiplyToSelf(d: Double): RealVector {
        return mapToSelf(fix2ndArgument(Multiply(), d))
    }

    /**
     * Subtract a value from each entry. Returns a new vector.
     * Does not change instance data.
     *
     * @param d Value to be subtracted.
     * @return `this` - `d`.
     */
    open fun mapSubtract(d: Double): RealVector? {
        return copy().mapSubtractToSelf(d)
    }

    /**
     * Subtract a value from each entry.
     * The instance is changed in-place.
     *
     * @param d Value to be subtracted.
     * @return `this`.
     */
    open fun mapSubtractToSelf(d: Double): RealVector? {
        return mapAddToSelf(-d)
    }

    /**
     * Divide each entry by the argument. Returns a new vector.
     * Does not change instance data.
     *
     * @param d Value to divide by.
     * @return `this` / `d`.
     */
    open fun mapDivide(d: Double): RealVector? {
        return copy().mapDivideToSelf(d)
    }

    /**
     * Divide each entry by the argument.
     * The instance is changed in-place.
     *
     * @param d Value to divide by.
     * @return `this`.
     */
    open fun mapDivideToSelf(d: Double): RealVector? {
        return mapToSelf(fix2ndArgument(Divide(), d))
    }

    /**
     * Find the orthogonal projection of this vector onto another vector.
     *
     * @param v vector onto which instance must be projected.
     * @return projection of the instance onto `v`.
     * @throws DimensionMismatchException if `v` is not the same size as
     * `this` vector.
     * @throws MathArithmeticException    if `this` or `v` is the null
     * vector
     */
    @Throws(DimensionMismatchException::class, MathArithmeticException::class)
    fun projection(v: RealVector): RealVector {
        val norm2 = v.dotProduct(v)
        if (norm2 == 0.0) {
            throw MathArithmeticException(LocalizedFormats.ZERO_NORM)
        }
        return v.mapMultiply(dotProduct(v) / v.dotProduct(v))
    }

    /**
     * Set all elements to a single value.
     *
     * @param value Single value to set for all elements.
     */
    open fun set(value: Double) {
        val it = iterator()
        while (it.hasNext()) {
            val e = it.next()
            e.value = value
        }
    }

    /**
     * Convert the vector to an array of `double`s.
     * The array is independent from this vector data: the elements
     * are copied.
     *
     * @return an array containing a copy of the vector elements.
     */
    open fun toArray(): DoubleArray {
        val dim = getDimension()
        val values = DoubleArray(dim)
        for (i in 0 until dim) {
            values[i] = getEntry(i)
        }
        return values
    }

    /**
     * Creates a unit vector pointing in the direction of this vector.
     * The instance is not changed by this method.
     *
     * @return a unit vector pointing in direction of this vector.
     * @throws MathArithmeticException if the norm is zero.
     */
    @Throws(MathArithmeticException::class)
    open fun unitVector(): RealVector? {
        val norm = getNorm()
        if (norm == 0.0) {
            throw MathArithmeticException(LocalizedFormats.ZERO_NORM)
        }
        return mapDivide(norm)
    }

    /**
     * Converts this vector into a unit vector.
     * The instance itself is changed by this method.
     *
     * @throws MathArithmeticException if the norm is zero.
     */
    @Throws(MathArithmeticException::class)
    open fun unitize() {
        val norm = getNorm()
        if (norm == 0.0) {
            throw MathArithmeticException(LocalizedFormats.ZERO_NORM)
        }
        mapDivideToSelf(norm)
    }


    /**
     * Generic dense iterator. Iteration is in increasing order
     * of the vector index.
     *
     *
     * Note: derived classes are required to return an [Iterator] that
     * returns non-null [Entry] objects as long as [Iterator.hasNext]
     * returns `true`.
     *
     * @return a dense iterator.
     */
    open operator fun iterator(): Iterator<Entry> {
        val dim = getDimension()
        return object : MutableIterator<Entry> {
            /** Current entry.  */
            private val e: Entry = Entry()

            /** Current index.  */
            private var i = 0

            /** {@inheritDoc}  */
            override fun hasNext(): Boolean {
                return i < dim
            }

            /** {@inheritDoc}  */
            override fun next(): Entry {
                return if (i < dim) {
                    e.index = i++
                    e
                } else {
                    throw NoDataException()
                }
            }

            /**
             * {@inheritDoc}
             *
             * @throws MathUnsupportedOperationException in all circumstances.
             */
            override fun remove() {
                throw MathUnsupportedOperationException()
            }
        }
    }

    /**
     * Acts as if implemented as:
     * <pre>
     * return copy().mapToSelf(function);
    </pre> *
     * Returns a new vector. Does not change instance data.
     *
     * @param function Function to apply to each entry.
     * @return a new vector.
     */
    open fun map(function: UnivariateFunction): RealVector {
        return copy().mapToSelf(function)
    }

    /**
     * Acts as if it is implemented as:
     * <pre>
     * Entry e = null;
     * for(Iterator<Entry> it = iterator(); it.hasNext(); e = it.next()) {
     * e.setValue(function.value(e.getValue()));
     * }
    </Entry></pre> *
     * Entries of this vector are modified in-place by this method.
     *
     * @param function Function to apply to each entry.
     * @return a reference to this vector.
     */
    open fun mapToSelf(function: UnivariateFunction): RealVector {
        val it = iterator()
        while (it.hasNext()) {
            val e = it.next()
            e.value = function.value(e.value)
        }
        return this
    }

    /**
     * Returns a new vector representing `a * this + b * y`, the linear
     * combination of `this` and `y`.
     * Returns a new vector. Does not change instance data.
     *
     * @param a Coefficient of `this`.
     * @param b Coefficient of `y`.
     * @param y Vector with which `this` is linearly combined.
     * @return a vector containing `a * this[i] + b * y[i]` for all
     * `i`.
     * @throws DimensionMismatchException if `y` is not the same size as
     * `this` vector.
     */
    @Throws(DimensionMismatchException::class)
    open fun combine(a: Double, b: Double, y: RealVector): RealVector? {
        return copy().combineToSelf(a, b, y)
    }

    /**
     * Updates `this` with the linear combination of `this` and
     * `y`.
     *
     * @param a Weight of `this`.
     * @param b Weight of `y`.
     * @param y Vector with which `this` is linearly combined.
     * @return `this`, with components equal to
     * `a * this[i] + b * y[i]` for all `i`.
     * @throws DimensionMismatchException if `y` is not the same size as
     * `this` vector.
     */
    @Throws(DimensionMismatchException::class)
    open fun combineToSelf(a: Double, b: Double, y: RealVector): RealVector? {
        checkVectorDimensions(y)
        for (i in 0 until getDimension()) {
            val xi = getEntry(i)
            val yi = y.getEntry(i)
            setEntry(i, a * xi + b * yi)
        }
        return this
    }

    /**
     * Visits (but does not alter) all entries of this vector in default order
     * (increasing index).
     *
     * @param visitor the visitor to be used to process the entries of this
     * vector
     * @return the value returned by [RealVectorPreservingVisitor.end]
     * at the end of the walk
     * @since 3.1
     */
    open fun walkInDefaultOrder(visitor: RealVectorPreservingVisitor): Double {
        val dim = getDimension()
        visitor.start(dim, 0, dim - 1)
        for (i in 0 until dim) {
            visitor.visit(i, getEntry(i))
        }
        return visitor.end()
    }

    /**
     * Visits (but does not alter) some entries of this vector in default order
     * (increasing index).
     *
     * @param visitor visitor to be used to process the entries of this vector
     * @param start   the index of the first entry to be visited
     * @param end     the index of the last entry to be visited (inclusive)
     * @return the value returned by [RealVectorPreservingVisitor.end]
     * at the end of the walk
     * @throws NumberIsTooSmallException if `end < start`.
     * @throws OutOfRangeException       if the indices are not valid.
     * @since 3.1
     */
    @Throws(NumberIsTooSmallException::class, OutOfRangeException::class)
    open fun walkInDefaultOrder(
        visitor: RealVectorPreservingVisitor,
        start: Int, end: Int
    ): Double {
        checkIndices(start, end)
        visitor.start(getDimension(), start, end)
        for (i in start..end) {
            visitor.visit(i, getEntry(i))
        }
        return visitor.end()
    }

    /**
     * Visits (but does not alter) all entries of this vector in optimized
     * order. The order in which the entries are visited is selected so as to
     * lead to the most efficient implementation; it might depend on the
     * concrete implementation of this abstract class.
     *
     * @param visitor the visitor to be used to process the entries of this
     * vector
     * @return the value returned by [RealVectorPreservingVisitor.end]
     * at the end of the walk
     * @since 3.1
     */
    open fun walkInOptimizedOrder(visitor: RealVectorPreservingVisitor): Double {
        return walkInDefaultOrder(visitor)
    }

    /**
     * Visits (but does not alter) some entries of this vector in optimized
     * order. The order in which the entries are visited is selected so as to
     * lead to the most efficient implementation; it might depend on the
     * concrete implementation of this abstract class.
     *
     * @param visitor visitor to be used to process the entries of this vector
     * @param start   the index of the first entry to be visited
     * @param end     the index of the last entry to be visited (inclusive)
     * @return the value returned by [RealVectorPreservingVisitor.end]
     * at the end of the walk
     * @throws NumberIsTooSmallException if `end < start`.
     * @throws OutOfRangeException       if the indices are not valid.
     * @since 3.1
     */
    @Throws(NumberIsTooSmallException::class, OutOfRangeException::class)
    open fun walkInOptimizedOrder(
        visitor: RealVectorPreservingVisitor,
        start: Int, end: Int
    ): Double {
        return walkInDefaultOrder(visitor, start, end)
    }

    /**
     * Visits (and possibly alters) all entries of this vector in default order
     * (increasing index).
     *
     * @param visitor the visitor to be used to process and modify the entries
     * of this vector
     * @return the value returned by [RealVectorChangingVisitor.end]
     * at the end of the walk
     * @since 3.1
     */
    open fun walkInDefaultOrder(visitor: RealVectorChangingVisitor): Double {
        val dim = getDimension()
        visitor.start(dim, 0, dim - 1)
        for (i in 0 until dim) {
            setEntry(i, visitor.visit(i, getEntry(i)))
        }
        return visitor.end()
    }

    /**
     * Visits (and possibly alters) some entries of this vector in default order
     * (increasing index).
     *
     * @param visitor visitor to be used to process the entries of this vector
     * @param start   the index of the first entry to be visited
     * @param end     the index of the last entry to be visited (inclusive)
     * @return the value returned by [RealVectorChangingVisitor.end]
     * at the end of the walk
     * @throws NumberIsTooSmallException if `end < start`.
     * @throws OutOfRangeException       if the indices are not valid.
     * @since 3.1
     */
    @Throws(NumberIsTooSmallException::class, OutOfRangeException::class)
    open fun walkInDefaultOrder(
        visitor: RealVectorChangingVisitor,
        start: Int, end: Int
    ): Double {
        checkIndices(start, end)
        visitor.start(getDimension(), start, end)
        for (i in start..end) {
            setEntry(i, visitor.visit(i, getEntry(i)))
        }
        return visitor.end()
    }

    /**
     * Visits (and possibly alters) all entries of this vector in optimized
     * order. The order in which the entries are visited is selected so as to
     * lead to the most efficient implementation; it might depend on the
     * concrete implementation of this abstract class.
     *
     * @param visitor the visitor to be used to process the entries of this
     * vector
     * @return the value returned by [RealVectorChangingVisitor.end]
     * at the end of the walk
     * @since 3.1
     */
    open fun walkInOptimizedOrder(visitor: RealVectorChangingVisitor): Double {
        return walkInDefaultOrder(visitor)
    }

    /**
     * Visits (and possibly change) some entries of this vector in optimized
     * order. The order in which the entries are visited is selected so as to
     * lead to the most efficient implementation; it might depend on the
     * concrete implementation of this abstract class.
     *
     * @param visitor visitor to be used to process the entries of this vector
     * @param start   the index of the first entry to be visited
     * @param end     the index of the last entry to be visited (inclusive)
     * @return the value returned by [RealVectorChangingVisitor.end]
     * at the end of the walk
     * @throws NumberIsTooSmallException if `end < start`.
     * @throws OutOfRangeException       if the indices are not valid.
     * @since 3.1
     */
    @Throws(NumberIsTooSmallException::class, OutOfRangeException::class)
    open fun walkInOptimizedOrder(
        visitor: RealVectorChangingVisitor,
        start: Int, end: Int
    ): Double {
        return walkInDefaultOrder(visitor, start, end)
    }

    /**
     *
     *
     * Test for the equality of two real vectors. If all coordinates of two real
     * vectors are exactly the same, and none are `NaN`, the two real
     * vectors are considered to be equal. `NaN` coordinates are
     * considered to affect globally the vector and be equals to each other -
     * i.e, if either (or all) coordinates of the real vector are equal to
     * `NaN`, the real vector is equal to a vector with all `NaN`
     * coordinates.
     *
     *
     *
     * This method *must* be overriden by concrete subclasses of
     * [RealVector] (the current implementation throws an exception).
     *
     *
     * @param other Object to test for equality.
     * @return `true` if two vector objects are equal, `false` if
     * `other` is null, not an instance of `RealVector`, or
     * not equal to this `RealVector` instance.
     * @throws MathUnsupportedOperationException if this method is not
     * overridden.
     */
    override fun equals(other: Any?): Boolean {
        throw MathUnsupportedOperationException()
    }

    /**
     * {@inheritDoc}. This method *must* be overriden by concrete
     * subclasses of [RealVector] (current implementation throws an
     * exception).
     *
     * @throws MathUnsupportedOperationException if this method is not
     * overridden.
     */
    override fun hashCode(): Int {
        throw MathUnsupportedOperationException()
    }

    /**
     * An entry in the vector.
     */
    open inner class Entry {
        /**
         * Get the index of the entry.
         *
         * @return the index of the entry.
         */
        /**
         * Set the index of the entry.
         *
         * @param index New index for the entry.
         */
        /**
         * Index of this entry.
         */
        var index = 0
        /**
         * Get the value of the entry.
         *
         * @return the value of the entry.
         */
        /**
         * Set the value of the entry.
         *
         * @param value New value for the entry.
         */
        open var value: Double
            get() = getEntry(index)
            set(value) {
                setEntry(index, value)
            }

        /**
         * Simple constructor.
         */
        init {
            index = 0
        }
    }

    /**
     * This class should rarely be used, but is here to provide
     * a default implementation of sparseIterator(), which is implemented
     * by walking over the entries, skipping those that are zero.
     *
     *
     * Concrete subclasses which are SparseVector implementations should
     * make their own sparse iterator, rather than using this one.
     *
     *
     * This implementation might be useful for ArrayRealVector, when expensive
     * operations which preserve the default value are to be done on the entries,
     * and the fraction of non-default values is small (i.e. someone took a
     * SparseVector, and passed it into the copy-constructor of ArrayRealVector)
     */
    protected inner class SparseEntryIterator : MutableIterator<Entry> {
        /**
         * Dimension of the vector.
         */
        private val dim: Int

        /**
         * Last entry returned by [.next].
         */
        private val current: Entry

        /**
         * Next entry for [.next] to return.
         */
        private val next: Entry

        /**
         * Advance an entry up to the next nonzero one.
         *
         * @param e entry to advance.
         */
        protected fun advance(e: Entry?) {
            if (e == null) {
                return
            }
            do {
                e.index = e.index + 1
            } while (e.index < dim && e.value == 0.0)
            if (e.index >= dim) {
                e.index = -1
            }
        }

        /**
         * {@inheritDoc}
         */
        override fun hasNext(): Boolean {
            return next.index >= 0
        }

        /**
         * {@inheritDoc}
         */
        override fun next(): Entry {
            val index = next.index
            if (index < 0) {
                throw NoDataException()
            }
            current.index = index
            advance(next)
            return current
        }

        /**
         * {@inheritDoc}
         *
         * @throws MathUnsupportedOperationException in all circumstances.
         */
        override fun remove() {
            throw MathUnsupportedOperationException()
        }

        /**
         * Simple constructor.
         */
        init {
            dim = getDimension()
            current = Entry()
            next = Entry()
            if (next.value == 0.0) {
                advance(next)
            }
        }
    }

}