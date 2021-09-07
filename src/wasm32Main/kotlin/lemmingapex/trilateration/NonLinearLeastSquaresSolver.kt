package lemmingapex.trilateration

import leastsquares.LeastSquaresFactory.create
import leastsquares.LeastSquaresOptimizer
import leastsquares.LeastSquaresOptimizer.Optimum
import linear.ArrayRealVector
import linear.DiagonalMatrix

/**
 * Solves a Trilateration problem with an instance of a
 * [LeastSquaresOptimizer]
 *
 * @author scott
 */
open class NonLinearLeastSquaresSolver(
    protected val function: TrilaterationFunction,
    private val leastSquaresOptimizer: LeastSquaresOptimizer
) {
    fun solve(
        target: DoubleArray?,
        weights: DoubleArray?,
        initialPoint: DoubleArray?,
        debugInfo: Boolean = false
    ): Optimum? {
        if (debugInfo) {
            println("Max Number of Iterations : " + MAXNUMBEROFITERATIONS)
        }
        val leastSquaresProblem = create( // function to be optimized
            function,  // target values at optimal point in least square equation
            // (x0+xi)^2 + (y0+yi)^2 + ri^2 = target[i]
            ArrayRealVector(target, false),
            ArrayRealVector(initialPoint, false),
            DiagonalMatrix(weights, true),
            null,
            MAXNUMBEROFITERATIONS,
            MAXNUMBEROFITERATIONS
        )
        return leastSquaresOptimizer.optimize(leastSquaresProblem)
    }

    fun solve(debugInfo: Boolean = false): Optimum? {
        val numberOfPositions = function.positions.size
        val positionDimension: Int = function.positions[0].size
        val initialPoint = DoubleArray(positionDimension)
        // initial point, use average of the vertices
        for (i in function.positions.indices) {
            val vertex = function.positions[i]
            for (j in vertex.indices) {
                initialPoint[j] += vertex[j]
            }
        }
        for (j in initialPoint.indices) {
            initialPoint[j] = initialPoint[j] / numberOfPositions
        }
        if (debugInfo) {
            val output = StringBuilder("initialPoint: ")
            for (i in initialPoint.indices) {
                output.append(initialPoint[i]).append(" ")
            }
            println(output.toString())
        }
        val target = DoubleArray(numberOfPositions)
        val distances = function.distances
        val weights = DoubleArray(target.size)
        for (i in target.indices) {
            target[i] = 0.0
            weights[i] = inverseSquareLaw(distances[i])
        }
        return solve(target, weights, initialPoint, debugInfo)
    }

    private fun inverseSquareLaw(distance: Double): Double {
        return 1 / (distance * distance)
    }

    companion object {
        protected const val MAXNUMBEROFITERATIONS = 1000
    }
}