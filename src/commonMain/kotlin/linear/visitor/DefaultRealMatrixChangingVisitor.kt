package linear.visitor

/**
 * Default implementation of the [RealMatrixChangingVisitor] interface.
 *
 *
 * This class is a convenience to create custom visitors without defining all
 * methods. This class provides default implementations that do nothing.
 *
 *
 * @author haokangkang
 * @since 2.0
 */
class DefaultRealMatrixChangingVisitor : RealMatrixChangingVisitor {
    /**
     * {@inheritDoc}
     */
    override fun start(
        rows: Int, columns: Int,
        startRow: Int, endRow: Int, startColumn: Int, endColumn: Int
    ) {
    }

    /**
     * {@inheritDoc}
     */
    override fun visit(row: Int, column: Int, value: Double): Double {
        return value
    }

    /**
     * {@inheritDoc}
     */
    override fun end(): Double {
        return 0.0
    }
}