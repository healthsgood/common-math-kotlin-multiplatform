package leastsquares

import leastsquares.LeastSquaresProblem.Evaluation
import linear.RealMatrix
import linear.RealVector

/**
 * Applies a dense weight matrix to an evaluation.
 *
 * @since 3.3
 */
internal class DenseWeightedEvaluation
/**
 * Create a weighted evaluation from an unweighted one.
 * // weight square root is square, nR=nC=number of observations
 * @param unweighted the evalutation before weights are applied
 * @param weightSqrt the matrix square root of the weight matrix
 */(
    /**
     * the unweighted evaluation
     */
    private val unweighted: Evaluation,
    /**
     * reference to the weight square root matrix
     */
    private val weightSqrt: RealMatrix
) : AbstractEvaluation(weightSqrt.getColumnDimension()) {
    /* apply weights */
    /**
     * {@inheritDoc}
     */
    override fun getJacobian(): RealMatrix? {
        return weightSqrt.multiply(unweighted.getJacobian())
    }

    /**
     * {@inheritDoc}
     */
    override fun getResiduals(): RealVector? {
        return weightSqrt.operate(unweighted.getResiduals())
    }
    /* delegate */
    /**
     * {@inheritDoc}
     */
    override fun getPoint(): RealVector? {
        return unweighted.getPoint()
    }
}