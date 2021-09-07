package leastsquares

import leastsquares.LeastSquaresProblem.Evaluation
import linear.ArrayRealVector
import linear.QRDecomposition
import linear.RealMatrix
import linear.RealVector
import util.FastMath.sqrt

/**
 * An implementation of [Evaluation] that is designed for extension. All of the
 * methods implemented here use the methods that are left unimplemented.
 *
 *
 * TODO cache results?
 *
 * @author haokangkang
 * @since 3.3
 */
abstract class AbstractEvaluation
/**
 * Constructor.
 *
 * @param observationSize the number of observation. Needed for [                        ][.].
 */
internal constructor(observationSize: Int) : Evaluation {
    /**
     * {@inheritDoc}
     */
    override fun getCovariances(threshold: Double): RealMatrix? {
        // Set up the Jacobian.
        val j = getJacobian()

        // Compute transpose(J)J.
        val jTj = j!!.transpose()!!.multiply(j)

        // Compute the covariances matrix.
        val solver = QRDecomposition(jTj!!, threshold).solver
        return solver!!.getInverse()
    }

    /**
     * {@inheritDoc}
     */
    override fun getSigma(covarianceSingularityThreshold: Double): RealVector? {
        val cov = getCovariances(covarianceSingularityThreshold)
        val nC = cov!!.getColumnDimension()
        val sig: RealVector = ArrayRealVector(nC)
        for (i in 0 until nC) {
            sig.setEntry(i, sqrt(cov.getEntry(i, i)))
        }
        return sig
    }

    /**
     * {@inheritDoc}
     */
    override fun getCost(): Double {
        val r = ArrayRealVector(getResiduals())
        return sqrt(r.dotProduct(r))
    }
}