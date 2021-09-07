package linear.visitor

/**
 * This interface defines a visitor for the entries of a vector. Visitors
 * implementing this interface may alter the entries of the vector being
 * visited.
 *
 * @author haokangkang
 * @since 3.1
 */
interface RealVectorChangingVisitor {
    /**
     * Start visiting a vector. This method is called once, before any entry
     * of the vector is visited.
     *
     * @param dimension the size of the vector
     * @param start     the index of the first entry to be visited
     * @param end       the index of the last entry to be visited (inclusive)
     */
    fun start(dimension: Int, start: Int, end: Int)

    /**
     * Visit one entry of the vector.
     *
     * @param index the index of the entry being visited
     * @param value the value of the entry being visited
     * @return the new value of the entry being visited
     */
    fun visit(index: Int, value: Double): Double

    /**
     * End visiting a vector. This method is called once, after all entries of
     * the vector have been visited.
     *
     * @return the value returned by
     * [RealVector.walkInDefaultOrder],
     * [RealVector.walkInDefaultOrder],
     * [RealVector.walkInOptimizedOrder]
     * or
     * [RealVector.walkInOptimizedOrder]
     */
    fun end(): Double
}