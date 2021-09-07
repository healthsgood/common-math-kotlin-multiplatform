package leastsquares

import crossjvm.ArrayUtil
import exception.ConvergenceException
import exception.util.LocalizedFormats
import leastsquares.LeastSquaresOptimizer.Optimum
import linear.ArrayRealVector
import linear.RealMatrix
import util.FastMath.abs
import util.FastMath.max
import util.FastMath.min
import util.FastMath.sqrt
import util.Precision

/**
 * This class solves a least-squares problem using the Levenberg-Marquardt
 * algorithm.
 *
 *
 * This implementation *should* work even for over-determined systems
 * (i.e. systems having more point than equations). Over-determined systems
 * are solved by ignoring the point which have the smallest impact according
 * to their jacobian column norm. Only the rank of the matrix and some loop bounds
 * are changed to implement this.
 *
 *
 * The resolution engine is a simple translation of the MINPACK [lmder](http://www.netlib.org/minpack/lmder.f) routine with minor
 * changes. The changes include the over-determined resolution, the use of
 * inherited convergence checker and the Q.R. decomposition which has been
 * rewritten following the algorithm described in the
 * P. Lascaux and R. Theodor book *Analyse numrique matricielle
 * applique  l'art de l'ingnieur*, Masson 1986.
 *
 * The authors of the original fortran version are:
 *
 *  * Argonne National Laboratory. MINPACK project. March 1980
 *  * Burton S. Garbow
 *  * Kenneth E. Hillstrom
 *  * Jorge J. More
 *
 * The redistribution policy for MINPACK is available [here](http://www.netlib.org/minpack/disclaimer), for convenience, it
 * is reproduced below.
 *
 * <table border="0" width="80%" cellpadding="10" align="center" bgcolor="#E0E0E0">
 * <tr><td>
 * Minpack Copyright Notice (1999) University of Chicago.
 * All rights reserved
</td></tr> *
 * <tr><td>
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *  1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *  1. Redistributions in binary form must reproduce the above
 * copyright notice, this list of conditions and the following
 * disclaimer in the documentation and/or other materials provided
 * with the distribution.
 *  1. The end-user documentation included with the redistribution, if any,
 * must include the following acknowledgment:
 * `This product includes software developed by the University of
 * Chicago, as Operator of Argonne National Laboratory.`
 * Alternately, this acknowledgment may appear in the software itself,
 * if and wherever such third-party acknowledgments normally appear.
 *  1. **WARRANTY DISCLAIMER. THE SOFTWARE IS SUPPLIED "AS IS"
 * WITHOUT WARRANTY OF ANY KIND. THE COPYRIGHT HOLDER, THE
 * UNITED STATES, THE UNITED STATES DEPARTMENT OF ENERGY, AND
 * THEIR EMPLOYEES: (1) DISCLAIM ANY WARRANTIES, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE
 * OR NON-INFRINGEMENT, (2) DO NOT ASSUME ANY LEGAL LIABILITY
 * OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR
 * USEFULNESS OF THE SOFTWARE, (3) DO NOT REPRESENT THAT USE OF
 * THE SOFTWARE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS, (4)
 * DO NOT WARRANT THAT THE SOFTWARE WILL FUNCTION
 * UNINTERRUPTED, THAT IT IS ERROR-FREE OR THAT ANY ERRORS WILL
 * BE CORRECTED.**
 *  1. **LIMITATION OF LIABILITY. IN NO EVENT WILL THE COPYRIGHT
 * HOLDER, THE UNITED STATES, THE UNITED STATES DEPARTMENT OF
 * ENERGY, OR THEIR EMPLOYEES: BE LIABLE FOR ANY INDIRECT,
 * INCIDENTAL, CONSEQUENTIAL, SPECIAL OR PUNITIVE DAMAGES OF
 * ANY KIND OR NATURE, INCLUDING BUT NOT LIMITED TO LOSS OF
 * PROFITS OR LOSS OF DATA, FOR ANY REASON WHATSOEVER, WHETHER
 * SUCH LIABILITY IS ASSERTED ON THE BASIS OF CONTRACT, TORT
 * (INCLUDING NEGLIGENCE OR STRICT LIABILITY), OR OTHERWISE,
 * EVEN IF ANY OF SAID PARTIES HAS BEEN WARNED OF THE
 * POSSIBILITY OF SUCH LOSS OR DAMAGES.**
 * </td></tr>
</table> *
 *
 * @author haokangkang
 * @since 3.3
 */
class LevenbergMarquardtOptimizer
/**
 * Default constructor.
 *
 *
 * The default values for the algorithm settings are:
 *
 *  * Initial step bound factor: 100
 *  * Cost relative tolerance: 1e-10
 *  * Parameters relative tolerance: 1e-10
 *  * Orthogonality tolerance: 1e-10
 *  * QR ranking threshold: [Precision.SAFE_MIN]
 *
 */ constructor(
    /**
     * Positive input variable used in determining the initial step bound.
     */
    val initialStepBoundFactor: Double = 100.0,
    /**
     * Desired relative error in the sum of squares.
     */
    val costRelativeTolerance: Double = 1e-10,
    /**
     * Desired relative error in the approximate solution parameters.
     */
    val parameterRelativeTolerance: Double = 1e-10,
    /**
     * Desired max cosine on the orthogonality between the function vector
     * and the columns of the jacobian.
     */
    val orthoTolerance: Double = 1e-10,
    /**
     * Threshold for QR ranking.
     */
    val rankingThreshold: Double = Precision.SAFE_MIN
) : LeastSquaresOptimizer {
    /* configuration parameters */
    /**
     * Gets the value of a tuning parameter.
     *
     * @return the parameter's value.
     * @see .withInitialStepBoundFactor
     */
    /**
     * Gets the value of a tuning parameter.
     *
     * @return the parameter's value.
     * @see .withCostRelativeTolerance
     */
    /**
     * Gets the value of a tuning parameter.
     *
     * @return the parameter's value.
     * @see .withParameterRelativeTolerance
     */
    /**
     * Gets the value of a tuning parameter.
     *
     * @return the parameter's value.
     * @see .withOrthoTolerance
     */
    /**
     * Gets the value of a tuning parameter.
     *
     * @return the parameter's value.
     * @see .withRankingThreshold
     */

    /**
     * @param newInitialStepBoundFactor Positive input variable used in
     * determining the initial step bound. This bound is set to the
     * product of initialStepBoundFactor and the euclidean norm of
     * `diag * x` if non-zero, or else to `newInitialStepBoundFactor`
     * itself. In most cases factor should lie in the interval
     * `(0.1, 100.0)`. `100` is a generally recommended value.
     * of the matrix is reduced.
     * @return a new instance.
     */
    fun withInitialStepBoundFactor(newInitialStepBoundFactor: Double): LevenbergMarquardtOptimizer {
        return LevenbergMarquardtOptimizer(
            newInitialStepBoundFactor,
            costRelativeTolerance,
            parameterRelativeTolerance,
            orthoTolerance,
            rankingThreshold
        )
    }

    /**
     * @param newCostRelativeTolerance Desired relative error in the sum of squares.
     * @return a new instance.
     */
    fun withCostRelativeTolerance(newCostRelativeTolerance: Double): LevenbergMarquardtOptimizer {
        return LevenbergMarquardtOptimizer(
            initialStepBoundFactor,
            newCostRelativeTolerance,
            parameterRelativeTolerance,
            orthoTolerance,
            rankingThreshold
        )
    }

    /**
     * @param newParRelativeTolerance Desired relative error in the approximate solution
     * parameters.
     * @return a new instance.
     */
    fun withParameterRelativeTolerance(newParRelativeTolerance: Double): LevenbergMarquardtOptimizer {
        return LevenbergMarquardtOptimizer(
            initialStepBoundFactor,
            costRelativeTolerance,
            newParRelativeTolerance,
            orthoTolerance,
            rankingThreshold
        )
    }

    /**
     * Modifies the given parameter.
     *
     * @param newOrthoTolerance Desired max cosine on the orthogonality between
     * the function vector and the columns of the Jacobian.
     * @return a new instance.
     */
    fun withOrthoTolerance(newOrthoTolerance: Double): LevenbergMarquardtOptimizer {
        return LevenbergMarquardtOptimizer(
            initialStepBoundFactor,
            costRelativeTolerance,
            parameterRelativeTolerance,
            newOrthoTolerance,
            rankingThreshold
        )
    }

    /**
     * @param newQRRankingThreshold Desired threshold for QR ranking.
     * If the squared norm of a column vector is smaller or equal to this
     * threshold during QR decomposition, it is considered to be a zero vector
     * and hence the rank of the matrix is reduced.
     * @return a new instance.
     */
    fun withRankingThreshold(newQRRankingThreshold: Double): LevenbergMarquardtOptimizer {
        return LevenbergMarquardtOptimizer(
            initialStepBoundFactor,
            costRelativeTolerance,
            parameterRelativeTolerance,
            orthoTolerance,
            newQRRankingThreshold
        )
    }

    /**
     * {@inheritDoc}
     */
    override fun optimize(problem: LeastSquaresProblem?): Optimum? {
        // Pull in relevant data from the problem as locals.
        val nR = problem!!.getObservationSize() // Number of observed data.
        val nC = problem.getParameterSize() // Number of parameters.
        // Counters.
        val iterationCounter = problem.iterationCounter
        val evaluationCounter = problem.evaluationCounter
        // Convergence criterion.
        val checker = problem.convergenceChecker

        // arrays shared with the other private methods
        val solvedCols = min(nR, nC)
        /* Parameters evolution direction associated with lmPar. */
        val lmDir = DoubleArray(nC)
        /* Levenberg-Marquardt parameter. */
        var lmPar = 0.0

        // local point
        var delta = 0.0
        var xNorm = 0.0
        val diag = DoubleArray(nC)
        val oldX = DoubleArray(nC)
        var oldRes = DoubleArray(nR)
        val qtf = DoubleArray(nR)
        val work1 = DoubleArray(nC)
        val work2 = DoubleArray(nC)
        val work3 = DoubleArray(nC)


        // Evaluate the function at the starting point and calculate its norm.
        evaluationCounter.incrementCount()
        //value will be reassigned in the loop
        var current = problem.evaluate(problem.getStart())
        var currentResiduals = current!!.getResiduals()!!.toArray()
        var currentCost = current.getCost()
        var currentPoint = current.getPoint()!!.toArray()

        // Outer loop.
        var firstIteration = true
        while (true) {
            iterationCounter.incrementCount()
            val previous = current

            // QR decomposition of the jacobian matrix
            val internalData = qrDecomposition(current!!.getJacobian(), solvedCols)
            val weightedJacobian = internalData.weightedJacobian
            val permutation = internalData.permutation
            val diagR = internalData.diagR
            val jacNorm = internalData.jacNorm

            //residuals already have weights applied
            var weightedResidual = currentResiduals
            for (i in 0 until nR) {
                qtf[i] = weightedResidual[i]
            }

            // compute Qt.res
            qTy(qtf, internalData)

            // now we don't need Q anymore,
            // so let jacobian contain the R matrix with its diagonal elements
            for (k in 0 until solvedCols) {
                val pk = permutation[k]
                weightedJacobian[k][pk] = diagR[pk]
            }
            if (firstIteration) {
                // scale the point according to the norms of the columns
                // of the initial jacobian
                xNorm = 0.0
                for (k in 0 until nC) {
                    var dk = jacNorm[k]
                    if (dk == 0.0) {
                        dk = 1.0
                    }
                    val xk = dk * currentPoint[k]
                    xNorm += xk * xk
                    diag[k] = dk
                }
                xNorm = sqrt(xNorm)

                // initialize the step bound delta
                delta = if (xNorm == 0.0) initialStepBoundFactor else initialStepBoundFactor * xNorm
            }

            // check orthogonality between function vector and jacobian columns
            var maxCosine = 0.0
            if (currentCost != 0.0) {
                for (j in 0 until solvedCols) {
                    val pj = permutation[j]
                    val s = jacNorm[pj]
                    if (s != 0.0) {
                        var sum = 0.0
                        for (i in 0..j) {
                            sum += weightedJacobian[i][pj] * qtf[i]
                        }
                        maxCosine = max(maxCosine, abs(sum) / (s * currentCost))
                    }
                }
            }
            if (maxCosine <= orthoTolerance) {
                // Convergence has been reached.
                return OptimumImpl(
                    current,
                    evaluationCounter.count,
                    iterationCounter.count
                )
            }

            // rescale if necessary
            for (j in 0 until nC) {
                diag[j] = max(diag[j], jacNorm[j])
            }

            // Inner loop.
            var ratio = 0.0
            while (ratio < 1.0e-4) {


                // save the state
                for (j in 0 until solvedCols) {
                    val pj = permutation[j]
                    oldX[pj] = currentPoint[pj]
                }
                val previousCost = currentCost
                var tmpVec = weightedResidual
                weightedResidual = oldRes
                oldRes = tmpVec

                // determine the Levenberg-Marquardt parameter
                lmPar = determineLMParameter(
                    qtf, delta, diag,
                    internalData, solvedCols,
                    work1, work2, work3, lmDir, lmPar
                )

                // compute the new point and the norm of the evolution direction
                var lmNorm = 0.0
                for (j in 0 until solvedCols) {
                    val pj = permutation[j]
                    lmDir[pj] = -lmDir[pj]
                    currentPoint[pj] = oldX[pj] + lmDir[pj]
                    val s = diag[pj] * lmDir[pj]
                    lmNorm += s * s
                }
                lmNorm = sqrt(lmNorm)
                // on the first iteration, adjust the initial step bound.
                if (firstIteration) {
                    delta = min(delta, lmNorm)
                }

                // Evaluate the function at x + p and calculate its norm.
                evaluationCounter.incrementCount()
                current = problem.evaluate(ArrayRealVector(currentPoint))
                currentResiduals = current!!.getResiduals()!!.toArray()
                currentCost = current.getCost()
                currentPoint = current.getPoint()!!.toArray()

                // compute the scaled actual reduction
                var actRed = -1.0
                if (0.1 * currentCost < previousCost) {
                    val r = currentCost / previousCost
                    actRed = 1.0 - r * r
                }

                // compute the scaled predicted reduction
                // and the scaled directional derivative
                for (j in 0 until solvedCols) {
                    val pj = permutation[j]
                    val dirJ = lmDir[pj]
                    work1[j] = 0.0
                    for (i in 0..j) {
                        work1[i] += weightedJacobian[i][pj] * dirJ
                    }
                }
                var coeff1 = 0.0
                for (j in 0 until solvedCols) {
                    coeff1 += work1[j] * work1[j]
                }
                val pc2 = previousCost * previousCost
                coeff1 /= pc2
                val coeff2 = lmPar * lmNorm * lmNorm / pc2
                val preRed = coeff1 + 2 * coeff2
                val dirDer = -(coeff1 + coeff2)

                // ratio of the actual to the predicted reduction
                ratio = if (preRed == 0.0) 0.0 else actRed / preRed

                // update the step bound
                if (ratio <= 0.25) {
                    var tmp = if (actRed < 0) 0.5 * dirDer / (dirDer + 0.5 * actRed) else 0.5
                    if (0.1 * currentCost >= previousCost || tmp < 0.1) {
                        tmp = 0.1
                    }
                    delta = tmp * min(delta, 10.0 * lmNorm)
                    lmPar /= tmp
                } else if (lmPar == 0.0 || ratio >= 0.75) {
                    delta = 2 * lmNorm
                    lmPar *= 0.5
                }

                // test for successful iteration.
                if (ratio >= 1.0e-4) {
                    // successful iteration, update the norm
                    firstIteration = false
                    xNorm = 0.0
                    for (k in 0 until nC) {
                        val xK = diag[k] * currentPoint[k]
                        xNorm += xK * xK
                    }
                    xNorm = sqrt(xNorm)

                    // tests for convergence.
                    if (checker != null && checker.converged(iterationCounter.count, previous, current)) {
                        return OptimumImpl(current, evaluationCounter.count, iterationCounter.count)
                    }
                } else {
                    // failed iteration, reset the previous values
                    currentCost = previousCost
                    for (j in 0 until solvedCols) {
                        val pj = permutation[j]
                        currentPoint[pj] = oldX[pj]
                    }
                    tmpVec = weightedResidual
                    weightedResidual = oldRes
                    oldRes = tmpVec
                    // Reset "current" to previous values.
                    current = previous
                }

                // Default convergence criteria.
                if (abs(actRed) <= costRelativeTolerance && preRed <= costRelativeTolerance && ratio <= 2.0 ||
                    delta <= parameterRelativeTolerance * xNorm
                ) {
                    return OptimumImpl(current!!, evaluationCounter.count, iterationCounter.count)
                }

                // tests for termination and stringent tolerances
                if (abs(actRed) <= TWO_EPS && preRed <= TWO_EPS && ratio <= 2.0) {
                    throw ConvergenceException(
                        LocalizedFormats.TOO_SMALL_COST_RELATIVE_TOLERANCE,
                        costRelativeTolerance
                    )
                } else if (delta <= TWO_EPS * xNorm) {
                    throw ConvergenceException(
                        LocalizedFormats.TOO_SMALL_PARAMETERS_RELATIVE_TOLERANCE,
                        parameterRelativeTolerance
                    )
                } else if (maxCosine <= TWO_EPS) {
                    throw ConvergenceException(
                        LocalizedFormats.TOO_SMALL_ORTHOGONALITY_TOLERANCE,
                        orthoTolerance
                    )
                }
            }
        }
    }

    /**
     * Determines the Levenberg-Marquardt parameter.
     *
     *
     * This implementation is a translation in Java of the MINPACK
     * [lmpar](http://www.netlib.org/minpack/lmpar.f)
     * routine.
     *
     * This method sets the lmPar and lmDir attributes.
     *
     * The authors of the original fortran function are:
     *
     *  * Argonne National Laboratory. MINPACK project. March 1980
     *  * Burton  S. Garbow
     *  * Kenneth E. Hillstrom
     *  * Jorge   J. More
     *
     *
     * Luc Maisonobe did the Java translation.
     *
     * @param qy           Array containing qTy.
     * @param delta        Upper bound on the euclidean norm of diagR * lmDir.
     * @param diag         Diagonal matrix.
     * @param internalData Data (modified in-place in this method).
     * @param solvedCols   Number of solved point.
     * @param work1        work array
     * @param work2        work array
     * @param work3        work array
     * @param lmDir        the "returned" LM direction will be stored in this array.
     * @param lmPar        the value of the LM parameter from the previous iteration.
     * @return the new LM parameter
     */
    private fun determineLMParameter(
        qy: DoubleArray, delta: Double, diag: DoubleArray,
        internalData: InternalData, solvedCols: Int,
        work1: DoubleArray, work2: DoubleArray, work3: DoubleArray,
        lmDir: DoubleArray, lmPar: Double
    ): Double {
        var lmPar = lmPar
        val weightedJacobian = internalData.weightedJacobian
        val permutation = internalData.permutation
        val rank = internalData.rank
        val diagR = internalData.diagR
        val nC: Int = weightedJacobian[0].size

        // compute and store in x the gauss-newton direction, if the
        // jacobian is rank-deficient, obtain a least squares solution
        for (j in 0 until rank) {
            lmDir[permutation[j]] = qy[j]
        }
        for (j in rank until nC) {
            lmDir[permutation[j]] = 0.0
        }
        for (k in rank - 1 downTo 0) {
            val pk = permutation[k]
            val ypk = lmDir[pk] / diagR[pk]
            for (i in 0 until k) {
                lmDir[permutation[i]] -= ypk * weightedJacobian[i][pk]
            }
            lmDir[pk] = ypk
        }

        // evaluate the function at the origin, and test
        // for acceptance of the Gauss-Newton direction
        var dxNorm = 0.0
        for (j in 0 until solvedCols) {
            val pj = permutation[j]
            val s = diag[pj] * lmDir[pj]
            work1[pj] = s
            dxNorm += s * s
        }
        dxNorm = sqrt(dxNorm)
        var fp = dxNorm - delta
        if (fp <= 0.1 * delta) {
            lmPar = 0.0
            return lmPar
        }

        // if the jacobian is not rank deficient, the Newton step provides
        // a lower bound, parl, for the zero of the function,
        // otherwise set this bound to zero
        var sum2: Double
        var parl = 0.0
        if (rank == solvedCols) {
            for (j in 0 until solvedCols) {
                val pj = permutation[j]
                work1[pj] *= diag[pj] / dxNorm
            }
            sum2 = 0.0
            for (j in 0 until solvedCols) {
                val pj = permutation[j]
                var sum = 0.0
                for (i in 0 until j) {
                    sum += weightedJacobian[i][pj] * work1[permutation[i]]
                }
                val s = (work1[pj] - sum) / diagR[pj]
                work1[pj] = s
                sum2 += s * s
            }
            parl = fp / (delta * sum2)
        }

        // calculate an upper bound, paru, for the zero of the function
        sum2 = 0.0
        for (j in 0 until solvedCols) {
            val pj = permutation[j]
            var sum = 0.0
            for (i in 0..j) {
                sum += weightedJacobian[i][pj] * qy[i]
            }
            sum /= diag[pj]
            sum2 += sum * sum
        }
        val gNorm = sqrt(sum2)
        var paru = gNorm / delta
        if (paru == 0.0) {
            paru = Precision.SAFE_MIN / min(delta, 0.1)
        }

        // if the input par lies outside of the interval (parl,paru),
        // set par to the closer endpoint
        lmPar = min(paru, max(lmPar, parl))
        if (lmPar == 0.0) {
            lmPar = gNorm / dxNorm
        }
        for (countdown in 10 downTo 0) {

            // evaluate the function at the current value of lmPar
            if (lmPar == 0.0) {
                lmPar = max(Precision.SAFE_MIN, 0.001 * paru)
            }
            val sPar = sqrt(lmPar)
            for (j in 0 until solvedCols) {
                val pj = permutation[j]
                work1[pj] = sPar * diag[pj]
            }
            determineLMDirection(qy, work1, work2, internalData, solvedCols, work3, lmDir)
            dxNorm = 0.0
            for (j in 0 until solvedCols) {
                val pj = permutation[j]
                val s = diag[pj] * lmDir[pj]
                work3[pj] = s
                dxNorm += s * s
            }
            dxNorm = sqrt(dxNorm)
            val previousFP = fp
            fp = dxNorm - delta

            // if the function is small enough, accept the current value
            // of lmPar, also test for the exceptional cases where parl is zero
            if (abs(fp) <= 0.1 * delta ||
                parl == 0.0 && fp <= previousFP && previousFP < 0
            ) {
                return lmPar
            }

            // compute the Newton correction
            for (j in 0 until solvedCols) {
                val pj = permutation[j]
                work1[pj] = work3[pj] * diag[pj] / dxNorm
            }
            for (j in 0 until solvedCols) {
                val pj = permutation[j]
                work1[pj] /= work2[j]
                val tmp = work1[pj]
                for (i in j + 1 until solvedCols) {
                    work1[permutation[i]] -= weightedJacobian[i][pj] * tmp
                }
            }
            sum2 = 0.0
            for (j in 0 until solvedCols) {
                val s = work1[permutation[j]]
                sum2 += s * s
            }
            val correction = fp / (delta * sum2)

            // depending on the sign of the function, update parl or paru.
            if (fp > 0) {
                parl = max(parl, lmPar)
            } else if (fp < 0) {
                paru = min(paru, lmPar)
            }

            // compute an improved estimate for lmPar
            lmPar = max(parl, lmPar + correction)
        }
        return lmPar
    }

    /**
     * Solve a*x = b and d*x = 0 in the least squares sense.
     *
     * This implementation is a translation in Java of the MINPACK
     * [qrsolv](http://www.netlib.org/minpack/qrsolv.f)
     * routine.
     *
     * This method sets the lmDir and lmDiag attributes.
     *
     * The authors of the original fortran function are:
     *
     *  * Argonne National Laboratory. MINPACK project. March 1980
     *  * Burton  S. Garbow
     *  * Kenneth E. Hillstrom
     *  * Jorge   J. More
     *
     *
     * Luc Maisonobe did the Java translation.
     *
     * @param qy           array containing qTy
     * @param diag         diagonal matrix
     * @param lmDiag       diagonal elements associated with lmDir
     * @param internalData Data (modified in-place in this method).
     * @param solvedCols   Number of sloved point.
     * @param work         work array
     * @param lmDir        the "returned" LM direction is stored in this array
     */
    private fun determineLMDirection(
        qy: DoubleArray, diag: DoubleArray,
        lmDiag: DoubleArray,
        internalData: InternalData,
        solvedCols: Int,
        work: DoubleArray,
        lmDir: DoubleArray
    ) {
        val permutation = internalData.permutation
        val weightedJacobian = internalData.weightedJacobian
        val diagR = internalData.diagR

        // copy R and Qty to preserve input and initialize s
        //  in particular, save the diagonal elements of R in lmDir
        for (j in 0 until solvedCols) {
            val pj = permutation[j]
            for (i in j + 1 until solvedCols) {
                weightedJacobian[i][pj] = weightedJacobian[j][permutation[i]]
            }
            lmDir[j] = diagR[pj]
            work[j] = qy[j]
        }

        // eliminate the diagonal matrix d using a Givens rotation
        for (j in 0 until solvedCols) {

            // prepare the row of d to be eliminated, locating the
            // diagonal element using p from the Q.R. factorization
            val pj = permutation[j]
            val dpj = diag[pj]
            if (dpj != 0.0) {
                ArrayUtil.fill(lmDiag, j + 1, lmDiag.size, 0.0)
            }
            lmDiag[j] = dpj

            //  the transformations to eliminate the row of d
            // modify only a single element of Qty
            // beyond the first n, which is initially zero.
            var qtbpj = 0.0
            for (k in j until solvedCols) {
                val pk = permutation[k]

                // determine a Givens rotation which eliminates the
                // appropriate element in the current row of d
                if (lmDiag[k] != 0.0) {
                    val sin: Double
                    val cos: Double
                    val rkk = weightedJacobian[k][pk]
                    if (abs(rkk) < abs(lmDiag[k])) {
                        val cotan = rkk / lmDiag[k]
                        sin = 1.0 / sqrt(1.0 + cotan * cotan)
                        cos = sin * cotan
                    } else {
                        val tan = lmDiag[k] / rkk
                        cos = 1.0 / sqrt(1.0 + tan * tan)
                        sin = cos * tan
                    }

                    // compute the modified diagonal element of R and
                    // the modified element of (Qty,0)
                    weightedJacobian[k][pk] = cos * rkk + sin * lmDiag[k]
                    val temp = cos * work[k] + sin * qtbpj
                    qtbpj = -sin * work[k] + cos * qtbpj
                    work[k] = temp

                    // accumulate the tranformation in the row of s
                    for (i in k + 1 until solvedCols) {
                        val rik = weightedJacobian[i][pk]
                        val temp2 = cos * rik + sin * lmDiag[i]
                        lmDiag[i] = -sin * rik + cos * lmDiag[i]
                        weightedJacobian[i][pk] = temp2
                    }
                }
            }

            // store the diagonal element of s and restore
            // the corresponding diagonal element of R
            lmDiag[j] = weightedJacobian[j][permutation[j]]
            weightedJacobian[j][permutation[j]] = lmDir[j]
        }

        // solve the triangular system for z, if the system is
        // singular, then obtain a least squares solution
        var nSing = solvedCols
        for (j in 0 until solvedCols) {
            if (lmDiag[j] == 0.0 && nSing == solvedCols) {
                nSing = j
            }
            if (nSing < solvedCols) {
                work[j] = 0.0
            }
        }
        if (nSing > 0) {
            for (j in nSing - 1 downTo 0) {
                val pj = permutation[j]
                var sum = 0.0
                for (i in j + 1 until nSing) {
                    sum += weightedJacobian[i][pj] * work[i]
                }
                work[j] = (work[j] - sum) / lmDiag[j]
            }
        }

        // permute the components of z back to components of lmDir
        for (j in lmDir.indices) {
            lmDir[permutation[j]] = work[j]
        }
    }

    /**
     * Decompose a matrix A as A.P = Q.R using Householder transforms.
     *
     * As suggested in the P. Lascaux and R. Theodor book
     * *Analyse numrique matricielle applique
     * l'art de l'ingnieur* (Masson, 1986), instead of representing
     * the Householder transforms with u<sub>k</sub> unit vectors such that:
     * <pre>
     * H<sub>k</sub> = I - 2u<sub>k</sub>.u<sub>k</sub><sup>t</sup>
    </pre> *
     * we use <sub>k</sub> non-unit vectors such that:
     * <pre>
     * H<sub>k</sub> = I - beta<sub>k</sub>v<sub>k</sub>.v<sub>k</sub><sup>t</sup>
    </pre> *
     * where v<sub>k</sub> = a<sub>k</sub> - alpha<sub>k</sub> e<sub>k</sub>.
     * The beta<sub>k</sub> coefficients are provided upon exit as recomputing
     * them from the v<sub>k</sub> vectors would be costly.
     *
     * This decomposition handles rank deficient cases since the tranformations
     * are performed in non-increasing columns norms order thanks to columns
     * pivoting. The diagonal elements of the R matrix are therefore also in
     * non-increasing absolute values order.
     *
     * @param jacobian   Weighted Jacobian matrix at the current point.
     * @param solvedCols Number of solved point.
     * @return data used in other methods of this class.
     * @throws ConvergenceException if the decomposition cannot be performed.
     */
    @Throws(ConvergenceException::class)
    private fun qrDecomposition(
        jacobian: RealMatrix?,
        solvedCols: Int
    ): InternalData {
        // Code in this class assumes that the weighted Jacobian is -(W^(1/2) J),
        // hence the multiplication by -1.
        val weightedJacobian = jacobian!!.scalarMultiply(-1.0)!!.getData()
        val nR = weightedJacobian!!.size
        val nC: Int = weightedJacobian[0]!!.size
        val permutation = IntArray(nC)
        val diagR = DoubleArray(nC)
        val jacNorm = DoubleArray(nC)
        val beta = DoubleArray(nC)

        // initializations
        for (k in 0 until nC) {
            permutation[k] = k
            var norm2 = 0.0
            for (i in 0 until nR) {
                val akk = weightedJacobian[i]!![k]
                norm2 += akk * akk
            }
            jacNorm[k] = sqrt(norm2)
        }

        // transform the matrix column after column
        for (k in 0 until nC) {

            // select the column with the greatest norm on active components
            var nextColumn = -1
            var ak2 = Double.NEGATIVE_INFINITY
            for (i in k until nC) {
                var norm2 = 0.0
                for (j in k until nR) {
                    val aki = weightedJacobian[j]!![permutation[i]]
                    norm2 += aki * aki
                }
                if (norm2.isInfinite() || norm2.isNaN()) {
                    throw ConvergenceException(
                        LocalizedFormats.UNABLE_TO_PERFORM_QR_DECOMPOSITION_ON_JACOBIAN,
                        nR, nC
                    )
                }
                if (norm2 > ak2) {
                    nextColumn = i
                    ak2 = norm2
                }
            }
            if (ak2 <= rankingThreshold) {
                return InternalData(weightedJacobian as Array<DoubleArray>, permutation, k, diagR, jacNorm, beta)
            }
            val pk = permutation[nextColumn]
            permutation[nextColumn] = permutation[k]
            permutation[k] = pk

            // choose alpha such that Hk.u = alpha ek
            val akk = weightedJacobian[k]!![pk]
            val alpha = if (akk > 0) -sqrt(ak2) else sqrt(ak2)
            val betak = 1.0 / (ak2 - akk * alpha)
            beta[pk] = betak

            // transform the current column
            diagR[pk] = alpha
            weightedJacobian[k]!![pk] -= alpha

            // transform the remaining columns
            for (dk in nC - 1 - k downTo 1) {
                var gamma = 0.0
                for (j in k until nR) {
                    gamma += weightedJacobian[j]!![pk] * weightedJacobian[j]!![permutation[k + dk]]
                }
                gamma *= betak
                for (j in k until nR) {
                    weightedJacobian[j]!![permutation[k + dk]] -= gamma * weightedJacobian[j]!![pk]
                }
            }
        }
        return InternalData(weightedJacobian as Array<DoubleArray>, permutation, solvedCols, diagR, jacNorm, beta)
    }

    /**
     * Compute the product Qt.y for some Q.R. decomposition.
     *
     * @param y            vector to multiply (will be overwritten with the result)
     * @param internalData Data.
     */
    private fun qTy(
        y: DoubleArray,
        internalData: InternalData
    ) {
        val weightedJacobian = internalData.weightedJacobian
        val permutation = internalData.permutation
        val beta = internalData.beta
        val nR = weightedJacobian.size
        val nC: Int = weightedJacobian[0].size
        for (k in 0 until nC) {
            val pk = permutation[k]
            var gamma = 0.0
            for (i in k until nR) {
                gamma += weightedJacobian[i][pk] * y[i]
            }
            gamma *= beta[pk]
            for (i in k until nR) {
                y[i] -= gamma * weightedJacobian[i][pk]
            }
        }
    }

    /**
     * Holds internal data.
     * This structure was created so that all optimizer fields can be "final".
     * Code should be further refactored in order to not pass around arguments
     * that will modified in-place (cf. "work" arrays).
     */
    private class InternalData
    /**
     * @param weightedJacobian Weighted Jacobian.
     * @param permutation      Columns permutation array.
     * @param rank             Rank of the Jacobian matrix.
     * @param diagR            Diagonal elements of the R matrix in the QR decomposition.
     * @param jacNorm          Norms of the columns of the jacobian matrix.
     * @param beta             Coefficients of the Householder transforms vectors.
     */(
        /**
         * Weighted Jacobian.
         */
        val weightedJacobian: Array<DoubleArray>,
        /**
         * Columns permutation array.
         */
        val permutation: IntArray,
        /**
         * Rank of the Jacobian matrix.
         */
        val rank: Int,
        /**
         * Diagonal elements of the R matrix in the QR decomposition.
         */
        val diagR: DoubleArray,
        /**
         * Norms of the columns of the jacobian matrix.
         */
        val jacNorm: DoubleArray,
        /**
         * Coefficients of the Householder transforms vectors.
         */
        val beta: DoubleArray
    )

    companion object {
        /**
         * Twice the "epsilon machine".
         */
        private const val TWO_EPS = 2 * Precision.EPSILON
    }
    /**
     * Construct an instance with all parameters specified.
     *
     * @param initialStepBoundFactor initial step bound factor
     * @param costRelativeTolerance  cost relative tolerance
     * @param parameterRelativeTolerance   parameters relative tolerance
     * @param orthoTolerance         orthogonality tolerance
     * @param rankingThreshold     threshold in the QR decomposition. Columns with a 2
     * norm less than this threshold are considered to be
     * all 0s.
     */
}