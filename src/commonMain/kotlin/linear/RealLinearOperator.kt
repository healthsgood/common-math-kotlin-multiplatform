package linear

import exception.DimensionMismatchException

/**
 * This class defines a linear operator operating on real (`double`)
 * vector spaces. No direct access to the coefficients of the underlying matrix
 * is provided.
 *
 *
 * The motivation for such an interface is well stated by
 * [Barrett et al. (1994)](#BARR1994):
 * <blockquote>
 * We restrict ourselves to iterative methods, which work by repeatedly
 * improving an approximate solution until it is accurate enough. These
 * methods access the coefficient matrix A of the linear system only via the
 * matrix-vector product y = A  x
 * (and perhaps z = A<sup>T</sup>  x). Thus the user need only
 * supply a subroutine for computing y (and perhaps z) given x, which permits
 * full exploitation of the sparsity or other special structure of A.
</blockquote> *
 * <br></br>
 *
 * <dl>
 * <dt><a name="BARR1994">Barret et al. (1994)</a></dt>
 * <dd>
 * R. Barrett, M. Berry, T. F. Chan, J. Demmel, J. M. Donato, J. Dongarra,
 * V. Eijkhout, R. Pozo, C. Romine and H. Van der Vorst,
 * *Templates for the Solution of Linear Systems: Building Blocks for
 * Iterative Methods*, SIAM
</dd> *
</dl> *
 *
 * @author haokangkang
 * @since 3.0
 */
abstract class RealLinearOperator {
    /**
     * Returns the dimension of the codomain of this operator.
     *
     * @return the number of rows of the underlying matrix
     */
    abstract fun getRowDimension(): Int

    /**
     * Returns the dimension of the domain of this operator.
     *
     * @return the number of columns of the underlying matrix
     */
    abstract fun getColumnDimension(): Int

    /**
     * Returns the result of multiplying `this` by the vector `x`.
     *
     * @param x the vector to operate on
     * @return the product of `this` instance with `x`
     * @throws DimensionMismatchException if the column dimension does not match
     * the size of `x`
     */
    @Throws(DimensionMismatchException::class)
    abstract fun operate(x: RealVector?): RealVector?
}