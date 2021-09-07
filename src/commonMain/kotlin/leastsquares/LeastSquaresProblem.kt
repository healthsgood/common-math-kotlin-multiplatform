package leastsquares

import exception.TooManyEvaluationsException
import leastsquares.LeastSquaresProblem.Evaluation
import linear.RealMatrix
import linear.RealVector
import linear.exception.SingularMatrixException
import optim.OptimizationProblem

/**
 * The data necessary to define a non-linear least squares problem.
 *
 *
 * Includes the observed values, computed model function, and
 * convergence/divergence criteria. Weights are implicit in [ ][Evaluation.getResiduals] and [Evaluation.getJacobian].
 *
 *
 *
 * Instances are typically either created progressively using a [ //LeastSquaresBuilder builder] or created at once using a [ factory][LeastSquaresFactory].
 *
 * //@see LeastSquaresBuilder
 *
 * @author haokangkang
 * @see LeastSquaresFactory
 *
 * @see LeastSquaresAdapter
 *
 * @since 3.3
 */
interface LeastSquaresProblem : OptimizationProblem<Evaluation?> {
    /**
     * Gets the initial guess.
     *
     * @return the initial guess values.
     */
    fun getStart(): RealVector?

    /**
     * Get the number of observations (rows in the Jacobian) in this problem.
     *
     * @return the number of scalar observations
     */
    fun getObservationSize(): Int

    /**
     * Get the number of parameters (columns in the Jacobian) in this problem.
     *
     * @return the number of scalar parameters
     */
    fun getParameterSize(): Int

    /**
     * Evaluate the model at the specified point.
     *
     * @param point the parameter values.
     * @return the model's value and derivative at the given point.
     * @throws exception.TooManyEvaluationsException if the maximal number of evaluations (of the model vector function) is
     * exceeded.
     */
    @Throws(TooManyEvaluationsException::class)
    fun evaluate(point: RealVector?): Evaluation?

    /**
     * An evaluation of a [LeastSquaresProblem] at a particular point. This class
     * also computes several quantities derived from the value and its Jacobian.
     */
    interface Evaluation {
        /**
         * Get the covariance matrix of the optimized parameters. <br></br> Note that this
         * operation involves the inversion of the `J<sup>T</sup>J` matrix,
         * where `J` is the Jacobian matrix. The `threshold` parameter is a
         * way for the caller to specify that the result of this computation should be
         * considered meaningless, and thus trigger an exception.
         *
         * @param threshold Singularity threshold.
         * @return the covariance matrix.
         * @throws SingularMatrixException if the covariance matrix cannot be computed (singular problem).
         */
        @Throws(SingularMatrixException::class)
        fun getCovariances(threshold: Double): RealMatrix?

        /**
         * Get an estimate of the standard deviation of the parameters. The returned
         * values are the square root of the diagonal coefficients of the covariance
         * matrix, `sd(a[i]) ~= sqrt(C[i][i])`, where `a[i]` is the optimized
         * value of the `i`-th parameter, and `C` is the covariance matrix.
         *
         * @param covarianceSingularityThreshold Singularity threshold (see [                                       ][.getCovariances]).
         * @return an estimate of the standard deviation of the optimized parameters
         * @throws SingularMatrixException if the covariance matrix cannot be computed.
         */
        fun getSigma(covarianceSingularityThreshold: Double): RealVector?

        /**
         * Get the weighted Jacobian matrix.
         *
         * @return the weighted Jacobian: W<sup>1/2</sup> J.
         * @throws exception.DimensionMismatchException if the Jacobian dimension does not match problem dimension.
         */
        fun getJacobian(): RealMatrix?

        /**
         * Get the cost.
         *
         * @return the cost.
         * @see .getResiduals
         */
        fun getCost(): Double

        /**
         * Get the weighted residuals. The residual is the difference between the
         * observed (target) values and the model (objective function) value. There is one
         * residual for each element of the vector-valued function. The raw residuals are
         * then multiplied by the square root of the weight matrix.
         *
         * @return the weighted residuals: W<sup>1/2</sup> K.
         * @throws exception.DimensionMismatchException if the residuals have the wrong length.
         */
        fun getResiduals(): RealVector?

        /**
         * Get the abscissa (independent variables) of this evaluation.
         *
         * @return the point provided to [.evaluate].
         */
        fun getPoint(): RealVector?
    }
}