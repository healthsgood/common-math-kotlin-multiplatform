package leastsquares

import exception.MathIllegalStateException
import exception.util.LocalizedFormats
import leastsquares.LeastSquaresProblem
import leastsquares.LeastSquaresProblem.Evaluation
import linear.DiagonalMatrix
import linear.EigenDecomposition
import linear.RealMatrix
import linear.RealVector
import optim.AbstractOptimizationProblem
import optim.ConvergenceChecker
import util.FastMath.sqrt

/**
 * A Factory for creating [LeastSquaresProblem]s.
 *
 * @author haokangkang
 * @since 3.3
 */
object LeastSquaresFactory {
    /**
     * Create a [org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem]
     * from the given elements. There will be no weights applied (unit weights).
     *
     * @param model          the model function. Produces the computed values.
     * @param observed       the observed (target) values
     * @param start          the initial guess.
     * @param weight         the weight matrix
     * @param checker        convergence checker
     * @param maxEvaluations the maximum number of times to evaluate the model
     * @param maxIterations  the maximum number to times to iterate in the algorithm
     * @param lazyEvaluation Whether the call to [Evaluation.evaluate]
     * will defer the evaluation until access to the value is requested.
     * @param paramValidator Model parameters validator.
     * @return the specified General Least Squares problem.
     * @since 3.4
     */
    fun create(
        model: MultivariateJacobianFunction,
        observed: RealVector,
        start: RealVector?,
        weight: RealMatrix?,
        checker: ConvergenceChecker<Evaluation?>?,
        maxEvaluations: Int,
        maxIterations: Int,
        lazyEvaluation: Boolean,
        paramValidator: ParameterValidator?
    ): LeastSquaresProblem {
        val p: LeastSquaresProblem = LocalLeastSquaresProblem(
            model,
            observed,
            start,
            checker,
            maxEvaluations,
            maxIterations,
            lazyEvaluation,
            paramValidator
        )
        return if (weight != null) {
            weightMatrix(p, weight)
        } else {
            p
        }
    }

    /**
     * Create a [org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem]
     * from the given elements. There will be no weights applied (unit weights).
     *
     * @param model          the model function. Produces the computed values.
     * @param observed       the observed (target) values
     * @param start          the initial guess.
     * @param checker        convergence checker
     * @param maxEvaluations the maximum number of times to evaluate the model
     * @param maxIterations  the maximum number to times to iterate in the algorithm
     * @return the specified General Least Squares problem.
     */
    fun create(
        model: MultivariateJacobianFunction,
        observed: RealVector,
        start: RealVector?,
        checker: ConvergenceChecker<Evaluation?>?,
        maxEvaluations: Int,
        maxIterations: Int
    ): LeastSquaresProblem {
        return create(
            model,
            observed,
            start,
            null,
            checker,
            maxEvaluations,
            maxIterations,
            false,
            null
        )
    }

    /**
     * Create a [org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem]
     * from the given elements.
     *
     * @param model          the model function. Produces the computed values.
     * @param observed       the observed (target) values
     * @param start          the initial guess.
     * @param weight         the weight matrix
     * @param checker        convergence checker
     * @param maxEvaluations the maximum number of times to evaluate the model
     * @param maxIterations  the maximum number to times to iterate in the algorithm
     * @return the specified General Least Squares problem.
     */

    fun create(
        model: MultivariateJacobianFunction,
        observed: RealVector,
        start: RealVector?,
        weight: RealMatrix,
        checker: ConvergenceChecker<Evaluation?>?,
        maxEvaluations: Int,
        maxIterations: Int
    ): LeastSquaresProblem {
        return weightMatrix(
            create(
                model,
                observed,
                start,
                checker,
                maxEvaluations,
                maxIterations
            ),
            weight
        )
    }

    /**
     * Apply a dense weight matrix to the [LeastSquaresProblem].
     *
     * @param problem the unweighted problem
     * @param weights the matrix of weights
     * @return a new [LeastSquaresProblem] with the weights applied. The original
     * `problem` is not modified.
     */
    fun weightMatrix(
        problem: LeastSquaresProblem?,
        weights: RealMatrix
    ): LeastSquaresProblem {
        val weightSquareRoot = squareRoot(weights)
        return object : LeastSquaresAdapter(problem!!) {
            /** {@inheritDoc}  */
            override fun evaluate(point: RealVector?): Evaluation? {
                return DenseWeightedEvaluation(super.evaluate(point)!!, weightSquareRoot)
            }
        }
    }

    /**
     * Computes the square-root of the weight matrix.
     *
     * @param m Symmetric, positive-definite (weight) matrix.
     * @return the square-root of the weight matrix.
     */
    private fun squareRoot(m: RealMatrix): RealMatrix {
        return if (m is DiagonalMatrix) {
            val dim = m.getRowDimension()
            val sqrtM: RealMatrix = DiagonalMatrix(dim)
            for (i in 0 until dim) {
                sqrtM.setEntry(i, i, sqrt(m.getEntry(i, i)))
            }
            sqrtM
        } else {
            val dec = EigenDecomposition(m)
            dec.squareRoot as RealMatrix
        }
    }

    /**
     * A private, "field" immutable (not "real" immutable) implementation of [ ].
     *
     * @since 3.3
     */
    private class LocalLeastSquaresProblem(
        /**
         * Model function.
         */
        private val model: MultivariateJacobianFunction,
        /**
         * Target values for the model function at optimum.
         */
        private val target: RealVector,
        /**
         * Initial guess.
         */
        private val start: RealVector?,
        checker: ConvergenceChecker<Evaluation?>?,
        maxEvaluations: Int,
        maxIterations: Int,
        /**
         * Whether to use lazy evaluation.
         */
        private val lazyEvaluation: Boolean,
        /**
         * Model parameters validator.
         */
        private val paramValidator: ParameterValidator?
    ) : AbstractOptimizationProblem<Evaluation?>(maxEvaluations, maxIterations, checker), LeastSquaresProblem {
        /**
         * {@inheritDoc}
         */
        override fun getObservationSize(): Int {
            return target.getDimension()
        }

        /**
         * {@inheritDoc}
         */
        override fun getParameterSize(): Int {
            return start!!.getDimension()
        }

        /**
         * {@inheritDoc}
         */
        override fun getStart(): RealVector? {
            return start?.copy()
        }

        /**
         * {@inheritDoc}
         */
        override fun evaluate(point: RealVector?): Evaluation? {
            // Copy so optimizer can change point without changing our instance.
            val p = if (paramValidator == null) point!!.copy() else paramValidator.validate(point!!.copy())!!
            return if (lazyEvaluation) {
                LazyUnweightedEvaluation(
                    model as ValueAndJacobianFunction,
                    target,
                    p
                )
            } else {
                // Evaluate value and jacobian in one function call.
                val value = model.value(p)
                UnweightedEvaluation(
                    value!!.first,
                    value.second,
                    target,
                    p
                )
            }
        }

        /**
         * Container with the model evaluation at a particular point.
         */
        private class UnweightedEvaluation(
            values: RealVector?,
            /**
             * Derivative at point.
             */
            private val jacobian: RealMatrix?,
            target: RealVector,
            /**
             * Point of evaluation.
             */
            private val point: RealVector
        ) : AbstractEvaluation(target.getDimension()) {
            /**
             * Computed residuals.
             */
            private val residuals: RealVector

            /**
             * {@inheritDoc}
             */
            override fun getJacobian(): RealMatrix? {
                return jacobian
            }

            /**
             * {@inheritDoc}
             */
            override fun getPoint(): RealVector? {
                return point
            }

            /**
             * {@inheritDoc}
             */
            override fun getResiduals(): RealVector? {
                return residuals
            }

            /**
             * Create an [Evaluation] with no weights.
             *
             * @param values   the computed function values
             * @param jacobian the computed function Jacobian
             * @param target   the observed values
             * @param point    the abscissa
             */
            init {
                residuals = target.subtract(values!!)
            }
        }

        /**
         * Container with the model *lazy* evaluation at a particular point.
         */
        private class LazyUnweightedEvaluation
        /**
         * Create an [Evaluation] with no weights.
         * // Safe to cast as long as we control usage of this class.
         * @param model  the model function
         * @param target the observed values
         * @param point  the abscissa
         */(
            /**
             * Model and Jacobian functions.
             */
            private val model: ValueAndJacobianFunction,
            /**
             * Target values for the model function at optimum.
             */
            private val target: RealVector,
            /**
             * Point of evaluation.
             */
            private val point: RealVector
        ) : AbstractEvaluation(target.getDimension()) {
            /**
             * {@inheritDoc}
             */
            override fun getJacobian(): RealMatrix? {
                return model.computeJacobian(point.toArray())
            }

            /**
             * {@inheritDoc}
             */
            override fun getPoint(): RealVector? {
                return point
            }

            /**
             * {@inheritDoc}
             */
            override fun getResiduals(): RealVector? {
                return target.subtract(model.computeValue(point.toArray())!!)
            }
        }

        /**
         * Create a [LeastSquaresProblem] from the given data.
         *
         * @param model          the model function
         * @param target         the observed data
         * @param start          the initial guess
         * @param checker        the convergence checker
         * @param maxEvaluations the allowed evaluations
         * @param maxIterations  the allowed iterations
         * @param lazyEvaluation Whether the call to [Evaluation.evaluate]
         * will defer the evaluation until access to the value is requested.
         * @param paramValidator Model parameters validator.
         */
        init {
            if (lazyEvaluation &&
                model !is ValueAndJacobianFunction
            ) {
                // Lazy evaluation requires that value and Jacobian
                // can be computed separately.
                throw MathIllegalStateException(
                    LocalizedFormats.INVALID_IMPLEMENTATION,
                    model::class.simpleName
                )
            }
        }
    }
}