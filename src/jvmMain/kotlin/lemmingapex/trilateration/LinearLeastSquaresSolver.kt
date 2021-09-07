package com.lemmingapex.trilateration

import linear.*

/**
 *
 * For testing only. A linear approach to solve the Trilateration problem.
 * see http://inside.mines.edu/~whereman/talks/TurgutOzal-11-Trilateration.pdf
 *
 * @author scott
 */
open class LinearLeastSquaresSolver(protected val function: TrilaterationFunction) {
    @JvmOverloads
    fun solve(debugInfo: Boolean = false): RealVector {
        val numberOfPositions = function.positions.size
        val positionDimension: Int = function.positions[0].size
        val Ad = Array<DoubleArray?>(numberOfPositions - 1) { DoubleArray(positionDimension) }

        // TODO: which reference position should be used?  currently using postion and distance in index 0.
        for (i in 1 until numberOfPositions) {
            val Adi = DoubleArray(positionDimension)
            for (j in 0 until positionDimension) {
                Adi[j] = function.positions[i][j] - function.positions[0][j]
            }
            Ad[i - 1] = Adi
        }
        if (debugInfo) {
            println(Array2DRowRealMatrix(Ad))
        }

        // reference point is function.getPositions()[0], with distance function.getDistances()[0]
        val referenceDistance = function.distances[0]
        val r0squared = referenceDistance * referenceDistance
        val bd = DoubleArray(numberOfPositions - 1)
        for (i in 1 until numberOfPositions) {
            val ri = function.distances[i]
            val risquared = ri * ri

            // find distance between ri and r0
            var di0squared = 0.0
            for (j in 0 until positionDimension) {
                val dij0j = function.positions[i][j] - function.positions[0][j]
                di0squared += dij0j * dij0j
            }
            bd[i - 1] = 0.5 * (r0squared - risquared + di0squared)
        }
        if (debugInfo) {
            println(ArrayRealVector(bd))
        }
        val A: RealMatrix = Array2DRowRealMatrix(Ad, false)
        val b: RealVector = ArrayRealVector(bd, false)
        val solver = QRDecomposition(A, 0.0).solver
        val x: RealVector?
        x = if (!solver!!.isNonSingular) {
            // bummer...
            ArrayRealVector(DoubleArray(positionDimension))
        } else {
            solver.solve(b)
        }
        return x!!.add(ArrayRealVector(function.positions[0]))
    }
}