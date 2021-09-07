package lemmingapex.trilateration

import leastsquares.MultivariateJacobianFunction
import linear.Array2DRowRealMatrix
import linear.ArrayRealVector
import linear.RealMatrix
import linear.RealVector
import util.Pair
import kotlin.math.max

/**
 * Models the Trilateration problem. This is a formulation for a nonlinear least
 * squares optimizer.
 *
 * @author scott
 */
class TrilaterationFunction(positions: Array<DoubleArray>, distances: DoubleArray) : MultivariateJacobianFunction {
    /**
     * Known positions of static nodes
     */
    val positions: Array<DoubleArray>

    /**
     * Euclidean distances from static nodes to mobile node
     */
    val distances: DoubleArray

    /**
     * Calculate and return Jacobian function Actually return initialized function
     *
     * Jacobian matrix, [i][j] at
     * J[i][0] = delta_[(x0-xi)^2 + (y0-yi)^2 - ri^2]/delta_[x0] at
     * J[i][1] = delta_[(x0-xi)^2 + (y0-yi)^2 - ri^2]/delta_[y0] partial derivative with respect to the parameters passed to value() method
     *
     * @param point for which to calculate the slope
     * @return Jacobian matrix for point
     */
    fun jacobian(point: RealVector?): RealMatrix {
        val pointArray = point!!.toArray()
        val jacobian = Array<DoubleArray?>(distances.size) { DoubleArray(pointArray.size) }
        for (i in jacobian.indices) {
            for (j in pointArray.indices) {
                jacobian[i]!![j] = 2 * pointArray[j] - 2 * positions[i][j]
            }
        }
        return Array2DRowRealMatrix(jacobian)
    }

    override fun value(point: RealVector?): Pair<RealVector?, RealMatrix?>? {

        // input
        val pointArray = point!!.toArray()

        // output
        val resultPoint = DoubleArray(distances.size)

        // compute least squares
        for (i in resultPoint.indices) {
            resultPoint[i] = 0.0
            // calculate sum, add to overall
            for (j in pointArray.indices) {
                resultPoint[i] += (pointArray[j] - positions[i][j]) * (pointArray[j] - positions[i][j])
            }
            resultPoint[i] -= distances[i] * distances[i]
        }
        val jacobian = jacobian(point)
        return Pair(ArrayRealVector(resultPoint), jacobian)
    }

    companion object {
        protected const val epsilon = 1E-7
    }

    init {
        require(positions.size >= 2) { "Need at least two positions." }
        require(positions.size == distances.size) { "The number of positions you provided, " + positions.size + ", does not match the number of distances, " + distances.size + "." }

        // bound distances to strictly positive domain
        for (i in distances.indices) {
            distances[i] = max(distances[i], epsilon)
        }
        val positionDimension: Int = positions[0].size
        for (i in 1 until positions.size) {
            require(positionDimension == positions[i].size) { "The dimension of all positions should be the same." }
        }
        this.positions = positions
        this.distances = distances
    }
}